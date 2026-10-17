#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import gzip,json
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
M14=ROOT/"StageM14_FinalSubmissionAssembly_v0_1"
M15A=ROOT/"StageM15_DCISEligibilityAudit_v0_1"
M15C=ROOT/"StageM15C_SubtypeSensitivity_v0_1"
M15D=ROOT/"StageM15D_PanTissueSpecificityDecision_v0_1"
M15E=ROOT/"StageM15E_AmbientRNAFeasibility_v0_1"
states={"Fibroblast":["FB1","FB2","FB3","FB4"],"Endothelial":["VEV","VEA","VEAT","VEC","LE1","LE2"],"Perivascular":["PV1","PV2","PV3","PV4","PV5"]}
floors=[1,5,10,20,30]

def bh(p):
 p=np.asarray(p,float); o=np.full(len(p),np.nan); ok=np.isfinite(p)
 if not ok.any(): return o
 v=p[ok]; idx=np.argsort(v); q=v[idx]*len(v)/np.arange(1,len(v)+1)
 q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(v));z[idx]=np.minimum(q,1);o[ok]=z
 return o

# 1) authoritative corrected stromal patient table
ps=pd.read_csv(M15A/"results/M15A_patient_state_roles_corrected_v0_1.tsv",sep="\t")
ps=ps.drop(columns=["role_original"],errors="ignore")
ps.to_csv(ST/"source_data/main/Figure2/patient_state_fraction_summary.tsv",sep="\t",index=False)

# Also relabel Supplementary Fig S4 patient confidence rows where patient IDs are present.
s4=ST/"source_data/supplementary/FigureS4_patient_confidence_v0_1.tsv"
if s4.exists():
 d=pd.read_csv(s4,sep="\t")
 dcis={"GSE195861:NCCBC2","GSE195861:NCCBC3","GSE195861:NCCBC5","GSE195861:NCCBC6","GSE195861:NCCBC7","GSE195861:NCCBC11","GSE195861:NCCBC14"}
 if "patient_uid" in d:
  mask=d.patient_uid.isin(dcis)
  if "role" in d:d.loc[mask,"role"]="reference_dcis"
 d.to_csv(s4,sep="\t",index=False)

# 2) corrected analytical denominators (Table S2)
tab=(ps.groupby(["lineage","role","gse"],as_index=False)
       .agg(n_patients=("patient_uid","nunique"),n_cells=("n_cells","sum"))
       .sort_values(["lineage","role","gse"]))
tab.to_csv(ST/"tables/TableS2_analytical_cohort_denominators.tsv",sep="\t",index=False)

# 3) complete floor-sensitivity raw state fraction table, corrected at low floors
rows=[]
for scope in ["all_cohorts","GSE161529_only"]:
 for floor in floors:
  for lin,ss in states.items():
   q=ps[(ps.lineage==lin)&ps.role.isin(["reference_normal","core_primary"])&(ps.n_cells>=floor)].copy()
   if scope=="GSE161529_only":q=q[q.gse=="GSE161529"]
   for st in ss:
    a=q.loc[q.role=="reference_normal",st].dropna().to_numpy(float)
    b=q.loc[q.role=="core_primary",st].dropna().to_numpy(float)
    U,p=mannwhitneyu(a,b,alternative="two-sided") if len(a) and len(b) else (np.nan,np.nan)
    rows.append([scope,floor,lin,st,len(a),len(b),np.median(a) if len(a) else np.nan,np.median(b) if len(b) else np.nan,
                 (np.median(b)-np.median(a)) if len(a) and len(b) else np.nan,U,p])
FS=pd.DataFrame(rows,columns=["scope","min_cells","lineage","state","n_normal","n_primary","median_normal","median_primary","delta_primary_minus_normal","U","p"])
FS["FDR"]=FS.groupby(["scope","min_cells","lineage"])["p"].transform(lambda x:bh(x))
FS.to_csv(ST/"source_data/main/Figure2/Figure2CDE_complete_state_composition_floor_sensitivity_v0_1.tsv",sep="\t",index=False)
FS.to_csv(ST/"tables/TableS4_composition_sensitivity_full.tsv",sep="\t",index=False)

# 4) Figure 1 scale: remove unused pan-tissue bar from main source
f1=ST/"source_data/main/Figure1/Figure1B_normal_reference_scale_v0_1.tsv"
q=pd.read_csv(f1,sep="\t")
q=q[q.layer!="CELLxGENE_pan_tissue_normal"].copy()
q.to_csv(f1,sep="\t",index=False)

# 5) IMC: replace two-contrast table/source with explicit three-contrast family
imc=pd.read_csv(M15A/"results/M15A_IMC_threeway_marker_tests_v0_1.tsv",sep="\t")
imc.to_csv(ST/"tables/TableS9_IMC_gate_sensitivity.tsv",sep="\t",index=False)
imc.to_csv(ST/"source_data/main/Figure5/Figure5C_IMC_contractile_perivascular_validation_v0_1.tsv",sep="\t",index=False)
# add patient-level values needed for revised Figure 5C
pd.read_csv(M15A/"results/M15A_IMC_patient_marker_medians_three_gates_v0_1.tsv",sep="\t").to_csv(
 ST/"source_data/main/Figure5/Figure5C_IMC_patient_marker_medians_threeway_v0_1.tsv",sep="\t",index=False)

# 6) subtype sensitivity as new Supplementary Table S16
g=pd.read_csv(M15C/"results/M15C_global_composition_subtype_tests_FDR_v0_1.tsv",sep="\t")
a=pd.read_csv(M15C/"results/M15C_mapper_robust_axis_tests_FDR_v0_1.tsv",sep="\t")
g["result_type"]="global_CLR_PERMANOVA"; a["result_type"]="aggregate_axis"
# harmonized publication table via common columns
cols=["result_type","mapper","lineage","comparison","group_a","group_b","n_a","n_b"]
G=g[cols+["PERMANOVA_p","PERMANOVA_FDR_across_lineages","PERMDISP_p"]].rename(columns={"PERMANOVA_p":"p","PERMANOVA_FDR_across_lineages":"FDR","PERMDISP_p":"dispersion_p"})
G["effect_or_axis"]=np.nan;G["delta_b_minus_a"]=np.nan
A=a[cols+["axis","delta_b_minus_a","p","axis_FDR_across_lineages"]].rename(columns={"axis_FDR_across_lineages":"FDR","axis":"effect_or_axis"})
A["dispersion_p"]=np.nan
out=pd.concat([G,A],ignore_index=True,sort=False)
out.to_csv(ST/"tables/TableS16_subtype_sensitivity_GSE161529.tsv",sep="\t",index=False)
out.to_csv(ST/"source_data/supplementary/TableS16_subtype_sensitivity_GSE161529.tsv",sep="\t",index=False)

# 7) Correct scANVI cross-mapper role concordance by relabeling 7 DCIS patients.
dcis={"GSE195861:NCCBC2","GSE195861:NCCBC3","GSE195861:NCCBC5","GSE195861:NCCBC6","GSE195861:NCCBC7","GSE195861:NCCBC11","GSE195861:NCCBC14"}
bench=ROOT/"StageM13_PostIndependentReview_v0_1/results/scanvi_benchmark"
cr=[]
for lin in states:
 p=bench/f"{lin}_WaveA_scanvi_predictions_v0_1.tsv.gz"
 d=pd.read_csv(p,sep="\t")
 d.loc[d.patient_uid.isin(dcis),"role"]="reference_dcis"
 for role in ["reference_normal","core_primary","reference_dcis"]:
  z=d[d.role==role]
  if len(z):
   cr.append([lin,role,len(z),float(z.state_concordant_historical.mean()),float(z.scanvi_maxprob.median()),float(z.scanvi_entropy.median())])
CR=pd.DataFrame(cr,columns=["lineage","role","n_cells","exact_state_concordance","median_scanvi_maxprob","median_scanvi_entropy"])
CR.to_csv(ST/"source_data/supplementary/FigureS12B_crossmapper_role_concordance_v0_1.tsv",sep="\t",index=False)
CR.to_csv(ST/"tables/TableS15F_crossmapper_role_concordance_DCIscorrected_v0_1.tsv",sep="\t",index=False)

# 8) retain explicit M15E provenance in a supplementary audit table
pd.read_csv(M15E/"results/M15E_public_data_provenance_v0_1.tsv",sep="\t").to_csv(
 ST/"tables/TableS17_ambientRNA_data_availability_v0_1.tsv",sep="\t",index=False)

summary={
 "status":"M15F_CORRECTED_SOURCE_BUILD_PASS",
 "corrected_TableS2_rows":int(len(tab)),
 "reference_dcis_rows_TableS2":int((tab.role=="reference_dcis").sum()),
 "floor20_principal_counts_unchanged":{
  lin:[int(((ps.lineage==lin)&(ps.role=="reference_normal")&(ps.n_cells>=20)).sum()),
       int(((ps.lineage==lin)&(ps.role=="core_primary")&(ps.n_cells>=20)).sum())] for lin in states},
 "pan_tissue_removed_from_main_Figure1_source":bool("CELLxGENE_pan_tissue_normal" not in q.layer.tolist()),
 "IMC_three_contrast_rows":int(len(imc)),
 "new_subtype_table_rows":int(len(out)),
 "scanvi_corrected_core_primary":{
   r.lineage:{"n_cells":int(r.n_cells),"concordance":float(r.exact_state_concordance)}
   for r in CR[CR.role=="core_primary"].itertuples()},
 "M14_untouched":True
}
(ST/"audit/M15F_source_build_gate_v0_1.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
print("\nTABLE S2 GSE195861\n",tab[tab.gse=="GSE195861"].to_string(index=False))
print("\nSCANVI CORE PRIMARY CORRECTED\n",CR[CR.role=="core_primary"].to_string(index=False))
