#!/usr/bin/env python3
import os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT=Path(os.environ.get('BC_PROJECT_ROOT', Path.cwd())).resolve()
S10=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage10a'
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage10c'

core=pd.read_csv(S10/'Stage10A_mechanistic_core_34_v0_1.tsv',sep='\t')
md=pd.read_csv(S10/'state_pseudobulk_metadata_v0_1.tsv',sep='\t')
z=np.load(S10/'state_pseudobulk_counts_v0_1.npz',allow_pickle=True)
X=z['X'].astype(float); genes=z['genes'].astype(str)
gidx={g:i for i,g in enumerate(genes)}
assert len(md)==X.shape[0]
md=md.reset_index().rename(columns={'index':'row_idx'})

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    v=p[ok]; n=len(v); o=np.argsort(v); r=v[o]
    q=np.minimum.accumulate((r*n/np.arange(1,n+1))[::-1])[::-1]
    tmp=np.empty(n); tmp[o]=np.minimum(q,1); out[ok]=tmp
    return out

def comp_cpm(row, comp):
    inds=[gidx.get(g,-1) for g in comp.split('_')]
    if any(i<0 for i in inds): return np.nan
    lib=float(md.loc[row,'library_total_counts'])
    if lib<=0: return np.nan
    return float(np.min(X[row,inds]/lib*1e6))

def one_score(pid,role,source,sstate,target,tstate,lig,rec,floor,gse_filter=None):
    q=md[(md.patient_uid==pid)&(md.role==role)]
    if gse_filter is not None: q=q[q.gse==gse_filter]
    s=q[(q.cell_type==source)&(q.state_affinity==sstate)&(q.n_cells>=floor)]
    t=q[(q.cell_type==target)&(q.state_affinity==tstate)&(q.n_cells>=floor)]
    if len(s)!=1 or len(t)!=1: return None
    sr=int(s.iloc[0].row_idx); tr=int(t.iloc[0].row_idx)
    lv=comp_cpm(sr,lig); rv=comp_cpm(tr,rec)
    if not (np.isfinite(lv) and np.isfinite(rv)): return None
    score=np.log2(np.sqrt(max(lv,0)*max(rv,0))+1)
    return score,int(s.iloc[0].n_cells),int(t.iloc[0].n_cells),str(s.iloc[0].gse)

rows=[]; samples=[]
for floor in [5,10,20,30]:
  for scope,gse_filter in [('GSE161529_only','GSE161529'),('all_cohorts',None)]:
    for r in core.itertuples(index=False):
      vals={}
      for role in ['reference_normal','core_primary']:
        pids=md.loc[md.role==role,'patient_uid'].drop_duplicates().tolist()
        arr=[]
        for pid in pids:
            o=one_score(pid,role,r.source,r.sender_state_top,r.target,r.receiver_state_top,
                        r.ligand_complex,r.receptor_complex,floor,gse_filter)
            if o is not None:
                arr.append((pid,*o))
                samples.append([r.core_rank,floor,scope,role,pid,o[3],o[0],o[1],o[2]])
        vals[role]=arr
      n=np.array([x[1] for x in vals['reference_normal']],float)
      p=np.array([x[1] for x in vals['core_primary']],float)
      if len(n)>=3 and len(p)>=3:
        stat,pv=mannwhitneyu(p,n,alternative='two-sided',method='auto')
        d=float(np.median(p)-np.median(n))
        expected=1 if r.rewiring=='gain' else -1
        concord=(np.sign(d)==expected and d!=0)
      else:
        pv=np.nan; d=np.nan; concord=np.nan
      rows.append([r.core_rank,r.rewiring,r.source,r.sender_state_top,r.target,r.receiver_state_top,
                   r.ligand_complex,r.receptor_complex,floor,scope,
                   len(n),len(p),float(np.median(n)) if len(n) else np.nan,
                   float(np.median(p)) if len(p) else np.nan,d,pv,concord,
                   r.delta_median_log2score,r.FDR_global])
R=pd.DataFrame(rows,columns=['core_rank','stage10A_rewiring','source','sender_state','target','receiver_state',
    'ligand_complex','receptor_complex','cell_floor','scope','n_normal','n_primary',
    'median_stateLR_normal','median_stateLR_primary','delta_state_matched_primary_minus_normal',
    'p_mannwhitney','direction_concordant_with_stage09','stage09_broad_delta','stage09_broad_FDR'])
R['FDR_within_floor_scope']=np.nan
for key,idx in R.groupby(['cell_floor','scope']).groups.items():
    R.loc[idx,'FDR_within_floor_scope']=bh(R.loc[idx,'p_mannwhitney'].to_numpy())
R['evaluable']=R.n_normal.ge(3)&R.n_primary.ge(3)
R['supported_FDR05']=R.evaluable & R.direction_concordant_with_stage09.fillna(False) & R.FDR_within_floor_scope.lt(0.05)
R.to_csv(OUT/'StateMatched_LR_primary_vs_normal_v0_1.tsv',sep='\t',index=False)

S=pd.DataFrame(samples,columns=['core_rank','cell_floor','scope','role','patient_uid','gse',
                                'log2_stateLR_score','sender_n_cells','receiver_n_cells'])
S.to_csv(OUT/'StateMatched_LR_patient_scores_v0_1.tsv.gz',sep='\t',index=False,compression='gzip')

summ=(R.groupby(['scope','cell_floor'])
      .agg(core_total=('core_rank','size'),evaluable=('evaluable','sum'),
           direction_concordant=('direction_concordant_with_stage09',lambda x:int((x==True).sum())),
           supported_FDR05=('supported_FDR05','sum'))
      .reset_index())
summ['concordance_among_evaluable']=summ.direction_concordant/summ.evaluable.replace(0,np.nan)
summ.to_csv(OUT/'StateMatched_LR_summary_v0_1.tsv',sep='\t',index=False)
print(summ.to_string(index=False))

print('\nMAIN FLOOR20 GSE161529')
q=R[(R.cell_floor==20)&(R.scope=='GSE161529_only')].copy()
print(q[['core_rank','stage10A_rewiring','sender_state','receiver_state','ligand_complex','receptor_complex',
         'n_normal','n_primary','delta_state_matched_primary_minus_normal','p_mannwhitney',
         'FDR_within_floor_scope','direction_concordant_with_stage09','supported_FDR05']].to_string(index=False))

print('\nDiscordant/evaluable floor20 main:')
print(q[q.evaluable & (q.direction_concordant_with_stage09==False)]
      [['core_rank','ligand_complex','receptor_complex','n_normal','n_primary',
        'delta_state_matched_primary_minus_normal','stage09_broad_delta']].to_string(index=False))
