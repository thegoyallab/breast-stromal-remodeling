#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/StageM22_ClinicalTrialStromalValidation_v0_1')
MAP=ROOT/'derived/GSE211434_Reed_state_mapping.tsv.gz'
OUT=ROOT/'derived'

def bh(pvals):
    p=np.asarray(pvals,float); q=np.full(len(p),np.nan)
    good=np.where(np.isfinite(p))[0]
    if not len(good):return q
    pg=p[good]; order=np.argsort(pg); adj=np.empty(len(good)); prev=1.
    for k in range(len(good)-1,-1,-1):
        i=order[k]; rank=k+1
        val=min(prev,pg[i]*len(good)/rank); adj[i]=val; prev=val
    q[good]=adj
    return q

m=pd.read_csv(MAP,sep='\t')
states={
 'Fibroblast':['FB1','FB2','FB3','FB4'],
 'Perivascular':['PV1','PV2','PV3','PV4','PV5'],
 'Endothelial':['LE1','LE2','VEA','VEAT','VEC','VEV']
}
rows=[]
for (pt,day,lin),g in m.groupby(['Patient','Day','mapper_lineage'],sort=False):
    cnt=g.nearest_normal_state.value_counts()
    d={'Patient':pt,'Day':int(day),'mapper_lineage':lin,'n_cells':len(g),
       'Treatment':g.Treatment.iloc[0],'TreatmentArm':g.TreatmentArm.iloc[0],
       'Response':g.Response.iloc[0],'ResponseYN':g.ResponseYN.iloc[0],
       'median_max_corr':g.max_corr.median(),'median_margin':g.margin.median()}
    for st in states[lin]: d[st]=int(cnt.get(st,0))
    n=len(g)
    for st in states[lin]:d[st+'_fraction']=d[st]/n
    if lin=='Fibroblast':
        a=d['FB1']+d['FB3']; b=d['FB2']+d['FB4']
        d['endpoint']='FB13_vs_FB24'; d['balance']=np.log((a+0.5)/(b+0.5))
        d['support_fraction']=a/n
    elif lin=='Perivascular':
        a=d['PV3']; b=d['PV1']+d['PV2']+d['PV4']+d['PV5']
        d['endpoint']='PV3_vs_rest'; d['balance']=np.log((a+0.5)/(b+0.5))
        d['support_fraction']=a/n
    else:
        a=d['LE2']; b=d['LE1']
        d['endpoint']='LE2_vs_LE1'; d['balance']=np.log((a+0.5)/(b+0.5))
        d['support_fraction']=(a/(a+b)) if (a+b)>0 else np.nan
    rows.append(d)
pb=pd.DataFrame(rows)
pb.to_csv(OUT/'GSE211434_Reed_state_patient_values.tsv',sep='\t',index=False)

tests=[]; loo=[]
for thr in [20,10]:
  for lin in ['Fibroblast','Perivascular','Endothelial']:
    d=pb[(pb.mapper_lineage==lin)&(pb.n_cells>=thr)].copy()
    endpoint=d.endpoint.iloc[0]
    for day in [14,180]:
      w=d[d.Day.isin([0,day])].pivot(index='Patient',columns='Day',values='balance').dropna()
      pts=w.index.tolist(); a=w[0].to_numpy(); b=w[day].to_numpy(); delta=b-a
      p=wilcoxon(b,a,alternative='two-sided',zero_method='wilcox').pvalue if len(delta) else np.nan
      med=float(np.median(delta)) if len(delta) else np.nan
      # paired raw supporting fraction medians
      sf=d[d.Day.isin([0,day])].pivot(index='Patient',columns='Day',values='support_fraction').reindex(pts)
      tests.append(dict(threshold=thr,lineage=lin,endpoint=endpoint,contrast=f'D0_D{day}',
                        test_family='paired_change',n=len(delta),median_balance_change=med,
                        median_support_fraction_D0=float(np.nanmedian(sf[0])),
                        median_support_fraction_Dx=float(np.nanmedian(sf[day])),
                        p=float(p)))
      if len(delta)>=3:
        for i,pt in enumerate(pts):
          lm=float(np.median(np.delete(delta,i)))
          loo.append(dict(threshold=thr,lineage=lin,endpoint=endpoint,contrast=f'D0_D{day}',
                          omitted_patient=pt,full_median_change=med,loo_median_change=lm,
                          same_direction=bool(np.sign(lm)==np.sign(med) or lm==0 or med==0)))
      mm=d.drop_duplicates('Patient').set_index('Patient')
      change=pd.DataFrame({'Patient':pts,'delta':delta})
      for family,field,A,B in [
          ('response_change','ResponseYN','Response','Non-response'),
          ('treatment_change','Treatment','Combination ribociclib','Letrozole alone')]:
          x=[r.delta for r in change.itertuples() if mm.loc[r.Patient,field]==A]
          y=[r.delta for r in change.itertuples() if mm.loc[r.Patient,field]==B]
          pu=mannwhitneyu(x,y,alternative='two-sided').pvalue if x and y else np.nan
          tests.append(dict(threshold=thr,lineage=lin,endpoint=endpoint,contrast=f'D0_D{day}',
                            test_family=family,n=len(delta),n_A=len(x),n_B=len(y),
                            group_A=A,group_B=B,median_change_A=float(np.median(x)) if x else np.nan,
                            median_change_B=float(np.median(y)) if y else np.nan,
                            difference_medians=float(np.median(x)-np.median(y)) if x and y else np.nan,
                            p=float(pu) if np.isfinite(pu) else np.nan))

res=pd.DataFrame(tests); res['FDR']=np.nan
for (thr,fam),idx in res.groupby(['threshold','test_family']).groups.items():
    ix=list(idx); res.loc[ix,'FDR']=bh(res.loc[ix,'p'].to_numpy())
loo=pd.DataFrame(loo)
rob=(loo.groupby(['threshold','lineage','endpoint','contrast']).same_direction.all()
     .rename('LOO_direction_robust').reset_index())
res=res.merge(rob,on=['threshold','lineage','endpoint','contrast'],how='left')
res.to_csv(OUT/'GSE211434_Reed_state_longitudinal_stats.tsv',sep='\t',index=False)
loo.to_csv(OUT/'GSE211434_Reed_state_LOO.tsv',sep='\t',index=False)

# state-fraction descriptive table
statecols=[c for c in pb.columns if c.endswith('_fraction')]
desc=pb.groupby(['mapper_lineage','Day'])[statecols+['median_max_corr','median_margin']].median(numeric_only=True).reset_index()
desc.to_csv(OUT/'GSE211434_Reed_state_descriptive_medians.tsv',sep='\t',index=False)

primary=res[res.threshold.eq(20)].sort_values(['test_family','lineage','contrast'])
print('PRIMARY_TESTS')
print(primary.to_string(index=False))
print('\nPRIMARY_FDR_LT_0.05')
print(primary[primary.FDR<0.05].to_string(index=False))
print('\nSENSITIVITY_PAIRED')
print(res[(res.threshold==10)&(res.test_family=='paired_change')].sort_values(['lineage','contrast']).to_string(index=False))
print('\nMAPPING_CONFIDENCE')
print(pb.groupby(['mapper_lineage','Day']).agg(n_patients=('Patient','nunique'),median_patient_max_corr=('median_max_corr','median'),median_patient_margin=('median_margin','median')).reset_index().to_string(index=False))
