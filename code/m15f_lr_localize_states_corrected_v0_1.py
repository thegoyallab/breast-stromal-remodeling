#!/home/pankaj/venvs/bc_scanpy/bin/python
import os
os.environ['OMP_NUM_THREADS']='8'
os.environ['OPENBLAS_NUM_THREADS']='8'
os.environ['MKL_NUM_THREADS']='8'
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage10a'
S9=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage09'
RES=OUT

cand=pd.read_csv(S9/'LR_mechanistic_shortlist_v0_1.tsv',sep='\t')
md=pd.read_csv(RES/'state_pseudobulk_metadata_v0_1.tsv',sep='\t')
z=np.load(RES/'state_pseudobulk_counts_v0_1.npz',allow_pickle=True)
X=z['X'].astype(np.float64)
genes=z['genes'].astype(str)
gidx={g:i for i,g in enumerate(genes)}

def bh(p):
    p=np.asarray(p,float)
    out=np.full(len(p),np.nan)
    ok=np.isfinite(p)
    if not ok.any(): return out
    vals=p[ok]; n=len(vals); o=np.argsort(vals); r=vals[o]
    q=r*n/np.arange(1,n+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; q=np.minimum(q,1.0)
    tmp=np.empty(n); tmp[o]=q; out[ok]=tmp
    return out

# fast group lookup
md=md.reset_index().rename(columns={'index':'row_idx'})
states={ct:sorted(md.loc[md.cell_type==ct,'state_affinity'].dropna().unique())
        for ct in md.cell_type.unique()}

def parts(c): return c.split('_')
def comp_expr(row, complex_name, lib):
    inds=[gidx.get(g,-1) for g in parts(complex_name)]
    if any(i<0 for i in inds) or lib<=0: return np.nan
    cpm=X[row,inds]/lib*1e6
    return float(np.min(cpm))

# Build total counts/lib per patient-role-celltype from state rows
keycols=['role','patient_uid','gse','cell_type']
totals={}
for key,sub in md.groupby(keycols,observed=True):
    rows=sub.row_idx.to_numpy(int)
    totals[key]=(X[rows,:].sum(axis=0), float(sub.library_total_counts.sum()),
                 int(sub.n_cells.sum()))

# complexes needed per role/lineage
tests=set()
for r in cand.itertuples(index=False):
    role='core_primary' if r.delta_median_log2score>0 else 'reference_normal'
    tests.add((role,r.source,r.ligand_complex))
    tests.add((role,r.target,r.receptor_complex))

rows_out=[]
for role,ct,comp in sorted(tests):
    subct=md[(md.role==role)&(md.cell_type==ct)].copy()
    for st in states[ct]:
        diffs=[]; se=[]; re=[]; det=[]; pids=[]; cohorts=[]
        for pid,subp in subct.groupby('patient_uid',observed=True):
            ss=subp[subp.state_affinity==st]
            if len(ss)!=1: continue
            rr=int(ss.iloc[0].row_idx); nstate=int(ss.iloc[0].n_cells)
            if nstate<5: continue
            key=(role,pid,str(ss.iloc[0].gse),ct)
            if key not in totals: continue
            total_counts,total_lib,total_cells=totals[key]
            rest_cells=total_cells-nstate
            rest_lib=total_lib-float(ss.iloc[0].library_total_counts)
            if rest_cells<10 or rest_lib<=0: continue
            inds=[gidx.get(g,-1) for g in parts(comp)]
            if any(i<0 for i in inds): continue
            state_cpm=X[rr,inds]/float(ss.iloc[0].library_total_counts)*1e6
            rest_counts=total_counts[inds]-X[rr,inds]
            rest_cpm=rest_counts/rest_lib*1e6
            sv=float(np.min(state_cpm)); rv=float(np.min(rest_cpm))
            diffs.append(np.log2(sv+1)-np.log2(rv+1))
            se.append(sv); re.append(rv); det.append(sv>0); pids.append(pid); cohorts.append(key[2])
        n=len(diffs)
        if n>=5 and np.any(np.asarray(diffs)!=0):
            try:
                p=float(wilcoxon(diffs,alternative='greater',zero_method='wilcox',
                                 correction=False,method='auto').pvalue)
            except Exception:
                p=np.nan
        else:
            p=np.nan
        rows_out.append([role,ct,comp,st,n,
                         float(np.median(diffs)) if n else np.nan,
                         float(np.median(se)) if n else np.nan,
                         float(np.median(re)) if n else np.nan,
                         float(np.mean(det)) if n else np.nan,
                         p,
                         ';'.join(sorted(set(cohorts))) if cohorts else ''])
loc=pd.DataFrame(rows_out,columns=['role','cell_type','complex','state','n_paired_patients',
    'median_log2_state_vs_rest','median_state_CPM','median_rest_CPM',
    'state_detection_fraction','p_enrichment','cohorts'])

# BH within each complex (across states) and lineage/role global
loc['FDR_complex']=np.nan
for key,idx in loc.groupby(['role','cell_type','complex']).groups.items():
    loc.loc[idx,'FDR_complex']=bh(loc.loc[idx,'p_enrichment'].to_numpy())
loc['FDR_lineage_role']=np.nan
for key,idx in loc.groupby(['role','cell_type']).groups.items():
    loc.loc[idx,'FDR_lineage_role']=bh(loc.loc[idx,'p_enrichment'].to_numpy())

loc['specific_supported']=(
    (loc.n_paired_patients>=5) &
    (loc.median_log2_state_vs_rest>=0.5) &
    (loc.state_detection_fraction>=0.5) &
    (loc.FDR_complex<0.05)
)
loc['localization_score']=(
    loc.median_log2_state_vs_rest.clip(lower=0).fillna(0) *
    np.sqrt(loc.n_paired_patients.clip(lower=0)) *
    loc.state_detection_fraction.fillna(0)
)
loc.to_csv(RES/'complex_state_localization_tests_v0_1.tsv.gz',sep='\t',
           index=False,compression='gzip')

# choose supported states, else best evaluable descriptive state
best={}
for key,sub in loc.groupby(['role','cell_type','complex'],observed=True):
    sp=sub[sub.specific_supported].sort_values(['localization_score','median_state_CPM'],
                                               ascending=False)
    if len(sp):
        best[key]=[(r.state,'supported',r.localization_score,r.FDR_complex,
                    r.median_log2_state_vs_rest,r.n_paired_patients,
                    r.state_detection_fraction) for r in sp.itertuples(index=False)]
    else:
        ev=sub[sub.n_paired_patients>=5].sort_values(['median_state_CPM','localization_score'],
                                                     ascending=False)
        if len(ev):
            r=next(ev.itertuples(index=False))
            best[key]=[(r.state,'descriptive',r.localization_score,r.FDR_complex,
                        r.median_log2_state_vs_rest,r.n_paired_patients,
                        r.state_detection_fraction)]
        else:
            best[key]=[]

# interaction-level state assignments
assign=[]
pair_rows=[]
for i,r in cand.reset_index(drop=True).iterrows():
    role='core_primary' if r.delta_median_log2score>0 else 'reference_normal'
    s=best.get((role,r.source,r.ligand_complex),[])
    t=best.get((role,r.target,r.receptor_complex),[])
    s_txt=';'.join([f'{x[0]}|{x[1]}|score={x[2]:.3f}|FDR={x[3]:.3g}|d={x[4]:.3f}|n={x[5]}' for x in s]) or 'not_evaluable'
    t_txt=';'.join([f'{x[0]}|{x[1]}|score={x[2]:.3f}|FDR={x[3]:.3g}|d={x[4]:.3f}|n={x[5]}' for x in t]) or 'not_evaluable'
    status=('both_supported' if s and t and s[0][1]=='supported' and t[0][1]=='supported'
            else 'one_or_both_descriptive' if s and t
            else 'not_fully_evaluable')
    assign.append([i,role,status,s_txt,t_txt,
                   s[0][0] if s else np.nan,t[0][0] if t else np.nan,
                   s[0][1] if s else 'not_evaluable',t[0][1] if t else 'not_evaluable',
                   s[0][2] if s else np.nan,t[0][2] if t else np.nan])
    for sx in (s if s else [(None,'not_evaluable',0,np.nan,np.nan,0,np.nan)]):
      for tx in (t if t else [(None,'not_evaluable',0,np.nan,np.nan,0,np.nan)]):
        if sx[0] is None or tx[0] is None: continue
        pair_score=float(r.priority_score)*(1+sx[2])*(1+tx[2])
        pair_rows.append([i,role,r.source,sx[0],r.target,tx[0],
                          r.ligand_complex,r.receptor_complex,
                          r.delta_median_log2score,r.FDR_global,r.priority_score,
                          sx[1],tx[1],sx[2],tx[2],pair_score])

A=pd.DataFrame(assign,columns=['candidate_id','localization_role','assignment_status',
    'sender_states_all','receiver_states_all','sender_state_top','receiver_state_top',
    'sender_assignment','receiver_assignment','sender_localization_score','receiver_localization_score'])
C=cand.reset_index(drop=True).reset_index().rename(columns={'index':'candidate_id'})
M=C.merge(A,on='candidate_id',how='left')
M.to_csv(RES/'LR_state_assignments_v0_1.tsv.gz',sep='\t',index=False,compression='gzip')

P=pd.DataFrame(pair_rows,columns=['candidate_id','role','source_lineage','sender_state',
    'target_lineage','receiver_state','ligand_complex','receptor_complex',
    'delta_median_log2score','FDR_global','stage09_priority_score',
    'sender_assignment','receiver_assignment','sender_localization_score',
    'receiver_localization_score','state_pair_priority_score'])
P=P.sort_values('state_pair_priority_score',ascending=False)
P.to_csv(RES/'LR_state_pair_expansion_v0_1.tsv.gz',sep='\t',index=False,compression='gzip')

# compact summary
summ=(M.groupby(['assignment_status'],dropna=False)
      .size().reset_index(name='interactions'))
summ.to_csv(RES/'Stage10A_assignment_summary_v0_1.tsv',sep='\t',index=False)
print(summ.to_string(index=False))
print('\nTop state-resolved pairs:')
print(P.head(50).to_string(index=False))
