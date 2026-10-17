#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd,re

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage10a'
RES=OUT
M=pd.read_csv(RES/'LR_state_assignments_v0_1.tsv.gz',sep='\t')
L=pd.read_csv(RES/'complex_state_localization_tests_v0_1.tsv.gz',sep='\t')

def family_lig(x):
    gs=x.split('_')
    if all(g.startswith('COL') for g in gs): return 'Collagen'
    if any(g.startswith('LAM') for g in gs): return 'Laminin'
    if any(g in {'JAG1','JAG2','DLL1','DLL3','DLL4'} for g in gs): return 'Notch_ligand'
    if any(g in {'FN1','COMP','HSPG2','AGRN','MMRN2','MFGE8'} for g in gs): return 'ECM_other'
    if any(g.startswith('S100') for g in gs): return 'S100'
    if any(g in {'VEGFA','VEGFB','VEGFC','PGF','ANGPT1','ANGPT2','ANGPTL4','EDN1'} for g in gs): return 'Vascular_GF'
    if any(g in {'SELE','SELL','SELP','ICAM1','VCAM1'} for g in gs): return 'Adhesion'
    if any(g.startswith('THBS') for g in gs): return 'Thrombospondin'
    if any(g.startswith('MMP') for g in gs): return 'MMP'
    if any(g.startswith('CALM') for g in gs): return 'Calmodulin'
    return 'Other'

def family_rec(x):
    gs=x.split('_')
    if any(g.startswith('ITGA') or g.startswith('ITGB') for g in gs): return 'Integrin'
    if any(g.startswith('NOTCH') for g in gs): return 'NOTCH'
    if any(g.startswith('SDC') for g in gs): return 'Syndecan'
    if any(g.startswith('ERBB') for g in gs): return 'ERBB'
    if any(g.startswith('PDGFR') for g in gs): return 'PDGFR'
    if any(g in {'CD93','CD248','PECAM1','ENG','ROBO4','NRP1','DYSF','AQP1'} for g in gs): return 'Vascular_receptor'
    if any(g in {'CD44','CD47','CD46','CD36'} for g in gs): return 'CD_family'
    return 'Other'

# enrich assignment table with top-state cohort count / FDR / effect
def loc_lookup(role,ct,comp,state):
    q=L[(L.role==role)&(L.cell_type==ct)&(L['complex']==comp)&(L.state==state)]
    if len(q)!=1: return {}
    r=q.iloc[0]
    return {'n_loc':r.n_paired_patients,'loc_fdr':r.FDR_complex,
            'loc_delta':r.median_log2_state_vs_rest,
            'loc_cohorts':r.cohorts,'n_cohorts':len(str(r.cohorts).split(';')) if pd.notna(r.cohorts) and str(r.cohorts) else 0}

rows=[]
for _,r in M[M.assignment_status=='both_supported'].iterrows():
    s=loc_lookup(r.localization_role,r.source,r.ligand_complex,r.sender_state_top)
    t=loc_lookup(r.localization_role,r.target,r.receptor_complex,r.receiver_state_top)
    d=r.to_dict()
    for k,v in s.items(): d['sender_'+k]=v
    for k,v in t.items(): d['receiver_'+k]=v
    d['ligand_family']=family_lig(r.ligand_complex)
    d['receptor_family']=family_rec(r.receptor_complex)
    d['circuit_family']=d['ligand_family']+'→'+d['receptor_family']
    d['state_pair']=str(r.sender_state_top)+'→'+str(r.receiver_state_top)
    rows.append(d)
B=pd.DataFrame(rows)
B['rewiring']='gain'
B.loc[B.delta_median_log2score<0,'rewiring']='loss'
# primary state specificity ideally represented in >=2 primary cohorts; label, do not delete
B['cross_cohort_state_support']=np.where(
    B.rewiring.eq('gain') & (B.sender_n_cohorts>=2) & (B.receiver_n_cohorts>=2),
    'two_or_more_primary_cohorts',
    np.where(B.rewiring.eq('gain'),'single_primary_cohort_state_localization','reference_normal_localization')
)
B.to_csv(RES/'LR_both_supported_state_resolved_v0_1.tsv',sep='\t',index=False)

# Collapse redundant interactions: retain max 2 within each state-pair x circuit family x rewiring.
B=B.sort_values('priority_score',ascending=False)
B['family_rank']=B.groupby(['rewiring','state_pair','circuit_family']).cumcount()+1
collapsed=B[B.family_rank<=2].copy()

# Final core target: 24 gained + 10 lost; prioritize cross-cohort state support in gains.
g=collapsed[collapsed.rewiring=='gain'].copy()
g['cc_bonus']=(g.cross_cohort_state_support=='two_or_more_primary_cohorts').astype(int)
g=g.sort_values(['cc_bonus','priority_score'],ascending=[False,False])
l=collapsed[collapsed.rewiring=='loss'].sort_values('priority_score',ascending=False)
core=pd.concat([g.head(24),l.head(10)],ignore_index=True)
core['core_rank']=np.arange(1,len(core)+1)
cols=['core_rank','rewiring','source','sender_state_top','target','receiver_state_top',
      'ligand_complex','receptor_complex','ligand_family','receptor_family','circuit_family',
      'delta_median_log2score','FDR_global','priority_score',
      'sender_n_loc','receiver_n_loc','sender_loc_fdr','receiver_loc_fdr',
      'sender_loc_delta','receiver_loc_delta','sender_n_cohorts','receiver_n_cohorts',
      'cross_cohort_state_support','de_details']
core=core[cols]
core.to_csv(RES/'Stage10A_mechanistic_core_34_v0_1.tsv',sep='\t',index=False)

# Core family/state summaries
fam=(core.groupby(['rewiring','circuit_family']).size().reset_index(name='n_core')
     .sort_values(['rewiring','n_core'],ascending=[True,False]))
fam.to_csv(RES/'Stage10A_core_family_summary_v0_1.tsv',sep='\t',index=False)
sp=(core.groupby(['rewiring','sender_state_top','receiver_state_top']).size().reset_index(name='n_core')
    .sort_values(['rewiring','n_core'],ascending=[True,False]))
sp.to_csv(RES/'Stage10A_core_statepair_summary_v0_1.tsv',sep='\t',index=False)

print('BOTH SUPPORTED:',len(B),'gained',sum(B.rewiring=='gain'),'lost',sum(B.rewiring=='loss'))
print('GAINED cross-cohort state-supported:',sum((B.rewiring=='gain')&(B.cross_cohort_state_support=='two_or_more_primary_cohorts')))
print('COLLAPSED DISTINCT:',len(collapsed))
print('CORE:',len(core))
print('\nCORE TABLE')
print(core.to_string(index=False))
print('\nFAMILIES')
print(fam.to_string(index=False))
