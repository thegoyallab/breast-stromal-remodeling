#!/home/pankaj/venvs/bc_scanpy/bin/python
import os
os.environ['OMP_NUM_THREADS']='8'
os.environ['OPENBLAS_NUM_THREADS']='8'
os.environ['MKL_NUM_THREADS']='8'
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
CCC=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage07e'
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage09'
RES=OUT

z=np.load(CCC/'patient_compartment_pseudobulk_CPM_v0_1.npz',allow_pickle=True)
X=z['X'].astype(np.float64)
genes=z['genes'].astype(str)
gidx={g:i for i,g in enumerate(genes)}
md=pd.read_csv(CCC/'patient_compartment_pseudobulk_metadata_v0_1.tsv',sep='\t')
md=md.reset_index().rename(columns={'index':'row_idx'})
res=pd.read_csv(CCC/'LIANA_consensus_resource_v1_10_0.tsv',sep='\t')
res=res.drop_duplicates(['ligand','receptor']).reset_index(drop=True)

cts=['Epithelial','Fibroblast','Endothelial','Perivascular']
roles=['reference_normal','core_primary']
min_cells=10

def bh(p):
    p=np.asarray(p,float)
    out=np.full(len(p),np.nan)
    ok=np.isfinite(p)
    if not ok.any(): return out
    vals=p[ok]; n=len(vals)
    order=np.argsort(vals)
    ranked=vals[order]
    q=ranked*n/np.arange(1,n+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    q=np.minimum(q,1)
    tmp=np.empty(n); tmp[order]=q
    out[ok]=tmp
    return out

def comp_expr(row_idx, complex_name):
    parts=complex_name.split('_')
    inds=[gidx.get(x,-1) for x in parts]
    if any(i<0 for i in inds): return np.nan
    return float(np.min(X[row_idx,inds]))

# patient-role-compartment lookup
eligible=md[md.n_cells>=min_cells].copy()
lookup={}
for role in roles:
    for ct in cts:
        sub=eligible[(eligible.role==role)&(eligible.cell_type==ct)]
        lookup[(role,ct)]={r.patient_uid:(int(r.row_idx),r.gse) for _,r in sub.iterrows()}

rows=[]
for source in cts:
  for target in cts:
    normal_pat=set(lookup[('reference_normal',source)]) & set(lookup[('reference_normal',target)])
    primary_pat=set(lookup[('core_primary',source)]) & set(lookup[('core_primary',target)])
    if len(normal_pat)<3 or len(primary_pat)<3:
        continue
    for ligand,receptor in res[['ligand','receptor']].itertuples(index=False,name=None):
        # Skip if any subunit unavailable
        if any(gidx.get(x,-1)<0 for x in ligand.split('_')+receptor.split('_')):
            continue
        def scores(role,pats):
            vals=[]; cohort=[]
            for pid in pats:
                srow,sgse=lookup[(role,source)][pid]
                trow,tgse=lookup[(role,target)][pid]
                lv=comp_expr(srow,ligand); rv=comp_expr(trow,receptor)
                vals.append(np.sqrt(max(lv,0)*max(rv,0)))
                cohort.append(sgse)
            return np.asarray(vals),np.asarray(cohort)
        nv,nc=scores('reference_normal',sorted(normal_pat))
        pv,pc=scores('core_primary',sorted(primary_pat))
        nl=np.log2(nv+1); pl=np.log2(pv+1)
        try:
            p=float(mannwhitneyu(pl,nl,alternative='two-sided').pvalue)
        except Exception:
            p=np.nan
        delta=float(np.median(pl)-np.median(nl))
        # descriptive primary-cohort direction consistency vs normal median
        ref=float(np.median(nl))
        signs=[]
        cohort_medians=[]
        for gse in sorted(set(pc)):
            vv=pl[pc==gse]
            if len(vv)>=3:
                d=float(np.median(vv)-ref)
                cohort_medians.append((gse,len(vv),d))
                if delta!=0 and d!=0:
                    signs.append(np.sign(d)==np.sign(delta))
        consistency=(float(np.mean(signs)) if signs else np.nan)
        rows.append([source,target,ligand,receptor,len(nl),len(pl),
                     float(np.median(nl)),float(np.median(pl)),delta,p,
                     consistency,len(signs),
                     ';'.join(f'{g}:{n}:{d:.4f}' for g,n,d in cohort_medians)])
    print('route',source,'->',target,'done',flush=True)

df=pd.DataFrame(rows,columns=['source','target','ligand_complex','receptor_complex',
    'n_normal','n_primary','median_log2score_normal','median_log2score_primary',
    'delta_median_log2score','p','primary_cohort_direction_consistency',
    'primary_cohorts_evaluable','primary_cohort_deltas'])
df['FDR_global']=bh(df.p)
df['FDR_route']=np.nan
for (s,t),idx in df.groupby(['source','target']).groups.items():
    df.loc[idx,'FDR_route']=bh(df.loc[idx,'p'].to_numpy())
df=df.sort_values(['FDR_route','FDR_global','p'])
df.to_csv(RES/'LR_pseudobulk_primary_vs_normal_validation_v0_1.tsv.gz',
          sep='\t',index=False,compression='gzip')

# Merge with LIANA primary / normal consensus evidence
keys=['source','target','ligand_complex','receptor_complex']
P=pd.read_csv(CCC/'liana_consensus_core_primary_v0_1.tsv.gz',sep='\t',
              usecols=keys+['lr_means','expr_prod','lrscore','specificity_rank','magnitude_rank'])
N=pd.read_csv(CCC/'liana_consensus_reference_normal_v0_1.tsv.gz',sep='\t',
              usecols=keys+['lr_means','expr_prod','lrscore','specificity_rank','magnitude_rank'])
P=P.rename(columns={c:c+'_primary' for c in P.columns if c not in keys})
N=N.rename(columns={c:c+'_normal' for c in N.columns if c not in keys})
m=df.merge(P,on=keys,how='left').merge(N,on=keys,how='left')
m['delta_liana_lr_means']=m['lr_means_primary']-m['lr_means_normal']
m['delta_liana_expr_prod']=m['expr_prod_primary']-m['expr_prod_normal']
m['direction_concordant']=np.sign(m['delta_median_log2score'])==np.sign(m['delta_liana_lr_means'])
# high-confidence: donor route FDR, global FDR, LIANA primary present,
# concordant direction, and cohort directional consistency >=0.75 when evaluable
m['high_confidence']=(m.FDR_route<0.05)&(m.FDR_global<0.05)&m.lr_means_primary.notna()&m.direction_concordant&(
    m.primary_cohort_direction_consistency.isna() | (m.primary_cohort_direction_consistency>=0.75))
m=m.sort_values(['high_confidence','FDR_global','FDR_route','p'],
                ascending=[False,True,True,True])
m.to_csv(RES/'LR_integrated_LIANA_pseudobulk_primary_vs_normal_v0_1.tsv.gz',
         sep='\t',index=False,compression='gzip')
summary=(m.groupby(['source','target'],observed=True)
         .agg(tested=('p','size'),route_FDR05=('FDR_route',lambda x:int((x<0.05).sum())),
              global_FDR05=('FDR_global',lambda x:int((x<0.05).sum())),
              high_confidence=('high_confidence','sum'))
         .reset_index())
summary.to_csv(RES/'LR_validation_route_summary_v0_1.tsv',sep='\t',index=False)
print(summary.to_string(index=False))
print('TOTAL high confidence',int(m.high_confidence.sum()))
