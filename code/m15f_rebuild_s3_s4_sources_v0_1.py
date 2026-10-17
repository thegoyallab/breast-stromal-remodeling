#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.stats import mannwhitneyu

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
SRC=ST/"source_data/supplementary"
D=pd.read_csv(ST/"source_data/main/Figure2/patient_state_fraction_summary.tsv",sep="\t")
metrics={
 "Fibroblast":["median_max_corr","FB1","FB2","FB3","FB4"],
 "Endothelial":["median_max_corr","VEV","VEA","VEAT","VEC","LE1","LE2"],
 "Perivascular":["median_max_corr","PV1","PV2","PV3","PV4","PV5"]}
contrasts=[
 ("reference_normal","reference_preneoplastic","Normal_vs_Preneoplastic"),
 ("reference_normal","core_primary","Normal_vs_Primary"),
 ("reference_normal","linked_progression_LN","Normal_vs_LN"),
 ("core_primary","linked_progression_LN","Primary_vs_LN")]

def bh(x):
 x=np.asarray(x,float);out=np.full(len(x),np.nan);ok=np.isfinite(x)
 if not ok.any():return out
 v=x[ok];o=np.argsort(v);q=v[o]*len(v)/np.arange(1,len(v)+1)
 q=np.minimum.accumulate(q[::-1])[::-1];z=np.empty(len(v));z[o]=np.minimum(q,1);out[ok]=z
 return out

# Rebuild S3a threshold/cohort sensitivity exactly as historical script but from corrected roles.
rows=[]
for scope in ["all_cohorts","GSE161529_only"]:
 ds=D.copy()
 if scope=="GSE161529_only": ds=ds[ds.gse=="GSE161529"]
 for floor in [1,5,10,20,30]:
  qf=ds[ds.n_cells>=floor]
  for lin in metrics:
   for ra,rb,cname in contrasts:
    for m in metrics[lin]:
     a=qf[(qf.lineage==lin)&(qf.role==ra)][m].dropna().to_numpy(float)
     b=qf[(qf.lineage==lin)&(qf.role==rb)][m].dropna().to_numpy(float)
     if len(a)>=3 and len(b)>=3:
      _,p=mannwhitneyu(b,a,alternative="two-sided")
      delta=float(np.median(b)-np.median(a))
     else:p=np.nan;delta=np.nan
     rows.append([scope,floor,lin,cname,m,len(a),len(b),
       float(np.median(a)) if len(a) else np.nan,float(np.median(b)) if len(b) else np.nan,delta,p])
R=pd.DataFrame(rows,columns=["scope","min_cells","lineage","contrast","metric","n_A","n_B","median_A","median_B","delta_B_minus_A","p"])
R["FDR_within_lineage_contrast"]=np.nan
for _,idx in R.groupby(["scope","min_cells","lineage","contrast"]).groups.items():
 R.loc[list(idx),"FDR_within_lineage_contrast"]=bh(R.loc[list(idx),"p"])
R.to_csv(SRC/"FigureS3_cell_threshold_and_cohort_sensitivity.tsv",sep="\t",index=False)

# Rebuild S3b mapping-quality exclusion using historical QC exclusion list and corrected roles.
qc=pd.read_csv(ROOT/"Reference_Arms/Normal_Reference/Mapping_Benchmark/stroma_broad_v0_1/patient_concordance_v0_1.tsv",sep="\t")
bad=qc.loc[qc.accuracy<0.90,"patient_uid"].drop_duplicates().tolist()
headline={"Fibroblast":["FB2","FB3","FB4"],"Endothelial":["VEV","VEA","VEC","LE1","LE2"],"Perivascular":["PV1","PV2","PV3","PV5"]}
rows=[]
for floor in [1,20,30]:
 for excl in [False,True]:
  d=D[D.n_cells>=floor].copy()
  if excl:d=d[~d.patient_uid.isin(bad)]
  for lin,mm in headline.items():
   for m in mm:
    a=d[(d.lineage==lin)&(d.role=="reference_normal")][m].dropna().to_numpy()
    b=d[(d.lineage==lin)&(d.role=="core_primary")][m].dropna().to_numpy()
    if len(a)>=3 and len(b)>=3:
     _,p=mannwhitneyu(b,a,alternative="two-sided");delta=float(np.median(b)-np.median(a))
    else:p=np.nan;delta=np.nan
    rows.append([floor,excl,lin,m,len(a),len(b),delta,p])
Q=pd.DataFrame(rows,columns=["min_cells","exclude_accuracy_lt0p90","lineage","metric","n_normal","n_primary","delta","p"])
Q["FDR"]=np.nan
for _,idx in Q.groupby(["min_cells","exclude_accuracy_lt0p90","lineage"]).groups.items():
 Q.loc[list(idx),"FDR"]=bh(Q.loc[list(idx),"p"])
Q.to_csv(SRC/"FigureS3_mapping_quality_exclusion_sensitivity.tsv",sep="\t",index=False)

# Correct S4 patient confidence role overlay.
s4=SRC/"FigureS4_patient_confidence_v0_1.tsv"
P=pd.read_csv(s4,sep="\t")
dcis={"GSE195861:NCCBC2","GSE195861:NCCBC3","GSE195861:NCCBC5","GSE195861:NCCBC6","GSE195861:NCCBC7","GSE195861:NCCBC11","GSE195861:NCCBC14"}
mask=P.patient_uid.isin(dcis)
assert mask.sum()==15,(mask.sum(),P.loc[mask,["lineage","patient_uid"]])
P.loc[mask,"stage06a_role"]="reference_dcis"
P.to_csv(s4,sep="\t",index=False)

# QA: floor20 headline normal-primary inputs identical to historical M14 for S3.
old=pd.read_csv(ROOT/"StageM14_FinalSubmissionAssembly_v0_1/source_data/supplementary/FigureS3_cell_threshold_and_cohort_sensitivity.tsv",sep="\t")
new=R[(R.scope=="all_cohorts")&(R.min_cells==20)&(R.contrast=="Normal_vs_Primary")]
old20=old[(old.scope=="all_cohorts")&(old.min_cells==20)&(old.contrast=="Normal_vs_Primary")]
key=["lineage","metric"]
j=new.set_index(key).join(old20.set_index(key),lsuffix="_new",rsuffix="_old")
assert np.allclose(j.delta_B_minus_A_new,j.delta_B_minus_A_old,equal_nan=True)
assert (j.n_A_new==j.n_A_old).all() and (j.n_B_new==j.n_B_old).all()

gate={"status":"M15F_S3_S4_SOURCE_REBUILD_PASS",
 "S3_rows":len(R),"S3_mapping_quality_rows":len(Q),
 "S4_DCIspatient_lineage_rows_relabelled":int(mask.sum()),
 "S3_floor20_normal_invasive_unchanged":True,
 "excluded_mapping_accuracy_lt0p90_n":len(bad)}
(ST/"audit/M15F_S3_S4_source_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
print("\nS3 floor1 corrected Normal-v-Primary max corr:")
print(R[(R.scope=="all_cohorts")&(R.min_cells==1)&(R.contrast=="Normal_vs_Primary")&(R.metric=="median_max_corr")][["lineage","n_B","median_B","p","FDR_within_lineage_contrast"]].to_string(index=False))
