from pathlib import Path
import pandas as pd, numpy as np, json
R=Path('/mnt/stk35data/BreastCancer_Census_v1/StageM22_ClinicalTrialStromalValidation_v0_1')
D=R/'derived'

# Table S20A: cohort design + applicability.
fel_map=pd.read_csv(D/'GSE211434_Reed_state_mapping.tsv.gz',sep='\t')
br_app=pd.read_csv(D/'BREAKFAST_Reed_applicability_gate.tsv',sep='\t')
rows=[]
wave_gates={'Fibroblast':(0.3168938385,0.00269698504),'Endothelial':(0.348450157,0.0020778835),'Perivascular':(0.3765755535,0.003566771765)}
for L,g in fel_map.groupby('mapper_lineage'):
    cg,mg=wave_gates[L]
    cm=float(g.max_corr.median()); mm=float(g.margin.median())
    rows.append(dict(cohort='FELINE GSE211434',trial='NCT02712723',setting='ER+ endocrine/CDK4/6',
                     timepoints='D0,D14,D180',paired_design='27 complete source cohort; lineage-specific >=20-cell paired n varies',
                     lineage=L,n_mapped_cells=len(g),median_max_corr=cm,max_corr_applicability_gate=cg,
                     median_margin=mm,margin_applicability_gate=mg,
                     fine_state_applicability='FAIL' if not(cm>=cg and mm>=mg) else 'PASS',
                     response_linkage='Exact source sensitive/resistant labels',
                     permitted_inference='Broad abundance/program + response-change; fine state audit only' if not(cm>=cg and mm>=mg) else 'All prespecified'))
for _,r in br_app.iterrows():
    rows.append(dict(cohort='BREAKFAST GSE279219',trial='NCT04248998',setting='Early TNBC preoperative therapy',
                     timepoints='T0,T1 (~14-20 d)',paired_design='10 complete public T0/T1 pairs',
                     lineage=r.lineage,n_mapped_cells=int(r.n_cells),median_max_corr=float(r.median_max_corr),
                     max_corr_applicability_gate=float(r.max_corr_gate),median_margin=float(r.median_margin),
                     margin_applicability_gate=float(r.margin_gate),
                     fine_state_applicability='PASS' if bool(r.applicability_pass) else 'FAIL',
                     response_linkage='Unresolved at molecular sample level',
                     permitted_inference='Paired treatment-state only; no response/arm claim'))
pd.DataFrame(rows).to_csv(R/'TableS20A_ClinicalTrial_Cohorts_and_Applicability.tsv',sep='\t',index=False)

# Table S20B: standardized primary/secondary endpoint tests.
out=[]
fa=pd.read_csv(D/'GSE211434_stromal_abundance_stats.tsv',sep='\t')
for _,r in fa[fa.threshold.eq(20)].iterrows():
    out.append(dict(cohort='FELINE GSE211434',endpoint_family='broad_lineage_abundance',lineage=r.lineage,
                    contrast=r.contrast,test_family=r.test_family,n=int(r.n),effect=r.get('effect_median_balance_change',np.nan),
                    group_A=r.get('group_A',np.nan),group_B=r.get('group_B',np.nan),
                    median_change_A=r.get('median_change_A',np.nan),median_change_B=r.get('median_change_B',np.nan),
                    p=r.p,FDR=r.FDR,reference_applicability='NOT_REQUIRED',interpretation='ACCEPTED'))
fp=pd.read_csv(D/'GSE211434_stromal_program_stats.tsv',sep='\t')
for _,r in fp[fp.threshold.eq(20)].iterrows():
    out.append(dict(cohort='FELINE GSE211434',endpoint_family='fixed_expression_program',lineage=r.lineage,
                    contrast=r.contrast,test_family=r.test_family,n=int(r.n),effect=r.get('median_change',np.nan),
                    group_A=r.get('group_A',np.nan),group_B=r.get('group_B',np.nan),
                    median_change_A=r.get('median_change_A',np.nan),median_change_B=r.get('median_change_B',np.nan),
                    p=r.p,FDR=r.FDR,reference_applicability='NOT_REQUIRED',interpretation='ACCEPTED'))
fs=pd.read_csv(D/'GSE211434_Reed_state_longitudinal_stats.tsv',sep='\t')
for _,r in fs[fs.threshold.eq(20)].iterrows():
    L=r.lineage
    cm=float(fel_map[fel_map.mapper_lineage.eq(L)].max_corr.median()); mm=float(fel_map[fel_map.mapper_lineage.eq(L)].margin.median())
    cg,mg=wave_gates[L];app=cm>=cg and mm>=mg
    out.append(dict(cohort='FELINE GSE211434',endpoint_family='frozen_Reed_state_balance',lineage=L,
                    contrast=r.contrast,test_family=r.test_family,n=int(r.n),effect=r.get('median_balance_change',np.nan),
                    group_A=r.get('group_A',np.nan),group_B=r.get('group_B',np.nan),
                    median_change_A=r.get('median_change_A',np.nan),median_change_B=r.get('median_change_B',np.nan),
                    p=r.p,FDR=r.FDR,reference_applicability='PASS' if app else 'FAIL',
                    interpretation='ACCEPTED' if app else 'AUDIT_ONLY_NOT_INTERPRETABLE'))
ba=pd.read_csv(D/'BREAKFAST_broad_abundance_stats.tsv',sep='\t')
for _,r in ba.iterrows():
    out.append(dict(cohort='BREAKFAST GSE279219',endpoint_family='broad_lineage_abundance',lineage=r.lineage,
                    contrast='T0_T1',test_family='paired_change',n=int(r.n),effect=r.median_balance_change,
                    p=r.p,FDR=r.FDR,reference_applicability='NOT_REQUIRED',interpretation='ACCEPTED'))
bp=pd.read_csv(D/'BREAKFAST_program_stats.tsv',sep='\t')
for _,r in bp.iterrows():
    out.append(dict(cohort='BREAKFAST GSE279219',endpoint_family='fixed_expression_program',lineage=r.lineage,
                    contrast='T0_T1',test_family='paired_change',n=int(r.n),effect=r.median_score_change,
                    p=r.p,FDR=r.FDR,reference_applicability='NOT_REQUIRED',interpretation='ACCEPTED'))
bs=pd.read_csv(D/'BREAKFAST_Reed_state_stats.tsv',sep='\t')
for _,r in bs.iterrows():
    out.append(dict(cohort='BREAKFAST GSE279219',endpoint_family='frozen_Reed_state_balance',lineage=r.lineage,
                    contrast='T0_T1',test_family='paired_change',n=int(r.n),effect=r.median_balance_change,
                    p=r.p,FDR=r.FDR,reference_applicability='PASS' if bool(r.applicability_pass) else 'FAIL',
                    interpretation='ACCEPTED' if bool(r.applicability_pass) else 'AUDIT_ONLY_NOT_INTERPRETABLE'))
tab=pd.DataFrame(out)
cols=['cohort','endpoint_family','lineage','contrast','test_family','n','effect','group_A','group_B','median_change_A','median_change_B','p','FDR','reference_applicability','interpretation']
for c in cols:
    if c not in tab:tab[c]=np.nan
tab[cols].to_csv(R/'TableS20B_ClinicalTrial_Stromal_Endpoint_Tests.tsv',sep='\t',index=False)

gate={
 'stage':'M22','date':'2026-10-16','status':'PASS_WITH_MIXED_BOUNDARY_RESULT',
 'parent_submission_package':'AtlasExpansion_PreSubmissionIntegrity_20261016_v0_97',
 'parent_package_modified':False,
 'cohorts_analyzed':['GSE211434/FELINE','GSE279219/BREAKFAST'],
 'cohort_selection_based_on_outcomes':False,
 'headline_positive_findings':['FELINE broad endothelial D0-D180 abundance increase: FDR 0.0107'],
 'negative_or_boundary_findings':[
   'No FELINE fixed stromal program passed FDR',
   'No FELINE response-associated stromal change passed FDR',
   'FELINE Reed fine-state mapping failed reference-domain applicability and is not interpreted',
   'BREAKFAST passed Reed applicability in all three lineages but no prespecified abundance/program/state endpoint passed FDR'
 ],
 'recommended_manuscript_role':'bounded supplementary longitudinal stress test',
 'title_change_recommended':False,'abstract_change_recommended':False,
 'main_claim_changed':False
}
(R/'M22_FINAL_GATE.json').write_text(json.dumps(gate,indent=2))
print(json.dumps(gate,indent=2))
