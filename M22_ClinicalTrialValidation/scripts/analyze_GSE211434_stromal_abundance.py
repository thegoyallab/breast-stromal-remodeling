import gzip,csv,math,json
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu

ROOT='/mnt/stk35data/BreastCancer_Census_v1/StageM22_ClinicalTrialStromalValidation_v0_1'
ann=ROOT+'/raw/GSE211434_cell_types.txt.gz'
rows=[]
with gzip.open(ann,'rt',newline='') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
meta={}
counts=Counter()
totals=Counter()
for r in rows:
    p=r['Patient']; d=int(r['Day']); ct=r['CellType']
    totals[(p,d)] += 1
    counts[(p,d,ct)] += 1
    meta[p]=(r['ResponseYN'],r['Treatment'],r['TreatmentArm'],r['Response'])

def bh(pvals):
    p=np.asarray(pvals,float); n=len(p); order=np.argsort(p); q=np.empty(n,float); prev=1.0
    for rank_idx in range(n-1,-1,-1):
        idx=order[rank_idx]; rank=rank_idx+1
        val=min(prev,p[idx]*n/rank); q[idx]=val; prev=val
    return q

lineages=['Fibroblasts','Pericytes','Endothelial cells']
records=[]
values=[]
for thr in [20,10]:
  for lin in lineages:
    for p in sorted(meta):
      for d in [0,14,180]:
        n=counts[(p,d,lin)]; total=totals[(p,d)]
        frac=n/total if total else np.nan
        bal=math.log((n+0.5)/(total-n+0.5)) if total else np.nan
        values.append(dict(threshold=thr,Patient=p,Day=d,lineage=lin,n_cells=n,total_cells=total,fraction=frac,balance=bal,
                           ResponseYN=meta[p][0],Treatment=meta[p][1],TreatmentArm=meta[p][2],ResponseRaw=meta[p][3],
                           eligible=n>=thr))
    for day in [14,180]:
      eligible=[p for p in sorted(meta) if counts[(p,0,lin)]>=thr and counts[(p,day,lin)]>=thr]
      a=np.array([math.log((counts[(p,0,lin)]+0.5)/(totals[(p,0)]-counts[(p,0,lin)]+0.5)) for p in eligible])
      b=np.array([math.log((counts[(p,day,lin)]+0.5)/(totals[(p,day)]-counts[(p,day,lin)]+0.5)) for p in eligible])
      af=np.array([counts[(p,0,lin)]/totals[(p,0)] for p in eligible])
      bf=np.array([counts[(p,day,lin)]/totals[(p,day)] for p in eligible])
      delta=b-a
      try: ps=wilcoxon(b,a,alternative='two-sided',zero_method='wilcox').pvalue
      except ValueError: ps=1.0
      records.append(dict(threshold=thr,lineage=lin,contrast=f'D0_D{day}',test_family='paired_change',
                          n=len(eligible),effect_median_balance_change=float(np.median(delta)),
                          median_fraction_D0=float(np.median(af)),median_fraction_Dx=float(np.median(bf)),
                          median_fraction_change=float(np.median(bf-af)),p=float(ps)))
      for groupvar,labelA,labelB in [('ResponseYN','Response','Non-response'),('Treatment','Combination ribociclib','Letrozole alone')]:
        da=[]; db=[]
        for i,p in enumerate(eligible):
          val=delta[i]
          g=meta[p][0] if groupvar=='ResponseYN' else meta[p][1]
          if g==labelA: da.append(val)
          elif g==labelB: db.append(val)
        pu=mannwhitneyu(da,db,alternative='two-sided').pvalue if da and db else np.nan
        records.append(dict(threshold=thr,lineage=lin,contrast=f'D0_D{day}',test_family='response_change' if groupvar=='ResponseYN' else 'treatment_change',
                            n=len(eligible),n_A=len(da),n_B=len(db),group_A=labelA,group_B=labelB,
                            median_change_A=float(np.median(da)) if da else np.nan,
                            median_change_B=float(np.median(db)) if db else np.nan,
                            effect_difference_medians=float(np.median(da)-np.median(db)) if da and db else np.nan,p=float(pu) if np.isfinite(pu) else np.nan))
res=pd.DataFrame(records)
# FDR per threshold + family across six lineage/timepoint tests
res['FDR']=np.nan
for (thr,fam),idx in res.groupby(['threshold','test_family']).groups.items():
    ix=list(idx); good=[i for i in ix if np.isfinite(res.loc[i,'p'])]
    if good:
        res.loc[good,'FDR']=bh(res.loc[good,'p'].values)
res.to_csv(ROOT+'/derived/GSE211434_stromal_abundance_stats.tsv',sep='\t',index=False)
pd.DataFrame(values).to_csv(ROOT+'/derived/GSE211434_stromal_abundance_patient_values.tsv.gz',sep='\t',index=False,compression='gzip')
primary=res[res.threshold.eq(20)].sort_values(['test_family','lineage','contrast'])
print(primary.to_string(index=False))
print('\nPRIMARY FDR<0.05')
print(primary[primary.FDR<0.05].to_string(index=False))
