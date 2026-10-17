#!/home/pankaj/venvs/bc_scanpy/bin/python
"""M17A4: re-run IMC protein validation in an independently defined mural-like phenotype.

Mural-like phenotype is based only on deposited cellPhenotype:
  Myofibroblasts OR Myofibroblasts PDPN^{+}
It does not use is_perivascular (vessel-contact) for lineage identity.
"""
from pathlib import Path
import os,zipfile,json
os.environ["OMP_NUM_THREADS"]="8"; os.environ["OPENBLAS_NUM_THREADS"]="8"; os.environ["MKL_NUM_THREADS"]="8"
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM17A_VascularSpatialClinicalPhospho_v0_1"
OUT=ST/"results"; AUD=ST/"audit"
ZIP=ROOT/"Parallel_Modalities/IMC_METABRIC_Danenberg/raw/METABRICBreastCancerIMCTME.zip"
MYO={"Myofibroblasts","Myofibroblasts PDPN^{+}"}
MARKERS=["SMA","PDGFRB","FSP1","Caveolin-1","CXCL12"]
MINCELLS=20

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    v=p[ok]; o=np.argsort(v); q=v[o]*len(v)/np.arange(1,len(v)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(v)); z[o]=np.minimum(q,1); out[ok]=z
    return out

pieces=[]
with zipfile.ZipFile(ZIP) as z:
  with z.open("toPublicRepository/SingleCells.csv") as f:
    for ch in pd.read_csv(f,usecols=["metabric_id","cellPhenotype","is_tumour","is_normal","is_dcis"]+MARKERS,
                          chunksize=250000):
      q=ch[ch.cellPhenotype.isin(MYO)].copy()
      if len(q):
        q["condition"]=np.select([q.is_normal.eq(1),q.is_dcis.eq(1),q.is_tumour.eq(1)],
                                 ["Normal","DCIS","Tumor"],default="Other")
        q=q[q.condition!="Other"]
        pieces.append(q)
d=pd.concat(pieces,ignore_index=True)
for m in MARKERS: d[m]=pd.to_numeric(d[m],errors="coerce")
g=d.groupby(["metabric_id","condition"])
p=g.size().rename("n_mural_like").reset_index()
p=p.merge(g[MARKERS].median().reset_index(),on=["metabric_id","condition"])
p["eligible_ge20"]=p.n_mural_like>=MINCELLS
p.to_csv(OUT/"M17A4_IMC_murallike_patient_marker_medians_v0_1.tsv",sep="\t",index=False)

rows=[]
for a,b,label in [("Tumor","Normal","Tumor_vs_Normal"),
                  ("DCIS","Normal","DCIS_vs_Normal"),
                  ("DCIS","Tumor","DCIS_vs_Tumor")]:
  for marker in MARKERS:
    q=p[p.eligible_ge20]
    x=q.loc[q.condition==a,marker].dropna().to_numpy(float)
    y=q.loc[q.condition==b,marker].dropna().to_numpy(float)
    if len(x)>=3 and len(y)>=3:
      U,pv=mannwhitneyu(x,y,alternative="two-sided")
      ma,mb=float(np.median(x)),float(np.median(y))
      rbc=float(2*U/(len(x)*len(y))-1)
    else: U=pv=ma=mb=rbc=np.nan
    rows.append([marker,label,a,b,len(x),len(y),ma,mb,ma-mb,U,pv,rbc])
R=pd.DataFrame(rows,columns=["marker","contrast","group_a","group_b","n_a","n_b","median_a","median_b",
                              "delta_a_minus_b","U","p","rank_biserial_a_over_b"])
R["FDR_all15"]=bh(R.p)
R["FDR_within_contrast5"]=np.nan
for c,ii in R.groupby("contrast").groups.items():
    R.loc[list(ii),"FDR_within_contrast5"]=bh(R.loc[list(ii),"p"])
R.to_csv(OUT/"M17A4_IMC_murallike_marker_tests_v0_1.tsv",sep="\t",index=False)

den=p[p.eligible_ge20].groupby("condition").agg(
    n_patients=("metabric_id","nunique"),n_mural_like_cells=("n_mural_like","sum")).reset_index()
den.to_csv(OUT/"M17A4_IMC_murallike_denominators_v0_1.tsv",sep="\t",index=False)

# Compare old semantics explicitly.
old=ROOT/"StageM15F_FinalScientificFreeze_v0_1/source_data/main/Figure5/Figure5C_IMC_contractile_perivascular_validation_v0_1.tsv"
oldr=pd.read_csv(old,sep="\t")
o=oldr[(oldr.gate=="exclude_is_epithelial") & oldr.marker.isin(["SMA","Caveolin-1"]) &
       oldr.contrast.isin(["Tumor_vs_Normal","DCIS_vs_Normal","DCIS_vs_Tumor"])].copy()
n=R[R.marker.isin(["SMA","Caveolin-1"])].copy()
cmp=o.merge(n,on=["marker","contrast"],suffixes=("_vesselcontactgate","_murallike"))
cmp.to_csv(OUT/"M17A4_old_vesselcontact_vs_murallike_comparison_v0_1.tsv",sep="\t",index=False)

def get(marker,contrast):
  return R[(R.marker==marker)&(R.contrast==contrast)].iloc[0]
s=get("SMA","Tumor_vs_Normal"); c=get("Caveolin-1","Tumor_vs_Normal")
d1=get("SMA","DCIS_vs_Normal"); d2=get("Caveolin-1","DCIS_vs_Normal")
gate={
 "status":"M17A4_MURALLIKE_IMC_REANALYSIS_COMPLETE",
 "semantic_correction":"Danenberg is_perivascular means vessel contact, not mural lineage identity.",
 "mural_like_definition":"cellPhenotype in {Myofibroblasts, Myofibroblasts PDPN+}",
 "min_cells_per_patient_condition":MINCELLS,
 "tumor_vs_normal_SMA":{"delta":float(s.delta_a_minus_b),"p":float(s.p),"FDR15":float(s.FDR_all15)},
 "tumor_vs_normal_Caveolin1":{"delta":float(c.delta_a_minus_b),"p":float(c.p),"FDR15":float(c.FDR_all15)},
 "DCIS_vs_normal_SMA":{"delta":float(d1.delta_a_minus_b),"p":float(d1.p),"FDR15":float(d1.FDR_all15)},
 "DCIS_vs_normal_Caveolin1":{"delta":float(d2.delta_a_minus_b),"p":float(d2.p),"FDR15":float(d2.FDR_all15)},
 "interpretation":"This analysis tests a myofibroblast-defined mural-like compartment; it is still not a direct pericyte lineage assay."
}
(AUD/"M17A4_murallike_imc_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
print("\nDENOMINATORS\n",den.to_string(index=False))
print("\nTESTS\n",R.to_string(index=False))
