#!/home/pankaj/venvs/bc_scanpy/bin/python
"""M17A1: lock an outcome-blind mural-like spatial phenotype in Danenberg IMC.

Important semantic correction:
- Danenberg 'is_perivascular' means *cell contacts a vessel mask*, not a cell lineage.
- Mural-like cells are therefore defined independently of spatial location as the deposited
  phenotypes 'Myofibroblasts' or 'Myofibroblasts PDPN^{+}'.
- Endothelial cells are the deposited phenotype 'Endothelial'.

Primary spatial endpoint:
  myofibroblast vessel-contact enrichment = log2((observed+0.5)/(expected+0.5)),
  with expected contact count from within-image random phenotype assignment among
  non-epithelial, non-endothelial cells while preserving the vessel-contact prevalence.

Secondary independent endpoint:
  myofibroblast-endothelial graph contact enrichment using the deposited undirected
  CellNeighbours graph, normalized to the exact random-mixing expectation within image.

No clinical outcome table is read in this script.
"""
from pathlib import Path
import os,zipfile,json,math
os.environ["OMP_NUM_THREADS"]="8"; os.environ["OPENBLAS_NUM_THREADS"]="8"; os.environ["MKL_NUM_THREADS"]="8"
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu,spearmanr

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM17A_VascularSpatialClinicalPhospho_v0_1"
OUT=ST/"results"; AUD=ST/"audit"
ZIP=ROOT/"Parallel_Modalities/IMC_METABRIC_Danenberg/raw/METABRICBreastCancerIMCTME.zip"
OUT.mkdir(parents=True,exist_ok=True); AUD.mkdir(parents=True,exist_ok=True)

MYO={"Myofibroblasts","Myofibroblasts PDPN^{+}"}
ENDO="Endothelial"
MIN_MYO=20
MIN_ENDO=20
PC=0.5

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    v=p[ok]; order=np.argsort(v); q=v[order]*len(v)/np.arange(1,len(v)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(v)); z[order]=np.minimum(q,1); out[ok]=z
    return out

# ---------- Single-cell-derived image counts and marker summaries ----------
rows=[]
with zipfile.ZipFile(ZIP) as z:
  with z.open("toPublicRepository/SingleCells.csv") as f:
    for ch in pd.read_csv(f,usecols=["ImageNumber","metabric_id","cellPhenotype","is_epithelial",
          "is_tumour","is_normal","is_dcis","is_perivascular","SMA","PDGFRB","Caveolin-1"],
          chunksize=250000):
      ch["condition"]=np.select([ch.is_normal.eq(1),ch.is_dcis.eq(1),ch.is_tumour.eq(1)],
                                ["Normal","DCIS","Tumor"],default="Other")
      for im,g in ch.groupby("ImageNumber",sort=False):
        # Dataset is image-sorted; nevertheless aggregation below will collapse chunk-split rows.
        pid=str(g.metabric_id.iloc[0]); cond=str(g.condition.iloc[0])
        myo=g.cellPhenotype.isin(MYO)
        endo=g.cellPhenotype.eq(ENDO)
        universe=g.is_epithelial.eq(0)&(~endo)
        contact=g.is_perivascular.eq(1)
        m=g.loc[myo,["SMA","PDGFRB","Caveolin-1"]].apply(pd.to_numeric,errors="coerce")
        rows.append(dict(ImageNumber=int(im),metabric_id=pid,condition=cond,
          n_cells=len(g),n_universe_non_epi_nonendo=int(universe.sum()),
          n_vesselcontact_universe=int((universe&contact).sum()),
          n_myo=int(myo.sum()),n_endo=int(endo.sum()),
          n_myo_vesselcontact=int((myo&contact).sum()),
          sum_SMA_myo=float(m.SMA.sum(skipna=True)),n_SMA_myo=int(m.SMA.notna().sum()),
          sum_PDGFRB_myo=float(m.PDGFRB.sum(skipna=True)),n_PDGFRB_myo=int(m.PDGFRB.notna().sum()),
          sum_CAV1_myo=float(m["Caveolin-1"].sum(skipna=True)),n_CAV1_myo=int(m["Caveolin-1"].notna().sum()),
          # medians cannot be safely combined across chunks; retained only as sums here.
        ))
img=pd.DataFrame(rows)
# Combine possible chunk-split images.
aggcols=[c for c in img.columns if c.startswith("n_") or c.startswith("sum_")]
img=img.groupby(["ImageNumber","metabric_id","condition"],as_index=False)[aggcols].sum()
# expected myofibroblast-vessel contact under within-image phenotype permutation
img["expected_myo_vesselcontact"]=np.where(img.n_universe_non_epi_nonendo>0,
    img.n_myo*img.n_vesselcontact_universe/img.n_universe_non_epi_nonendo,np.nan)

# ---------- independent graph endpoint ----------
# CellNeighbours is undirected but can contain repeated direction records; canonicalize.
edge_counts={}
seen=set()
with zipfile.ZipFile(ZIP) as z:
  with z.open("toPublicRepository/CellNeighbours.csv") as f:
    for ch in pd.read_csv(f,usecols=["ImageNumber","from","to","from_is_epithelial","to_is_epithelial",
                       "from_cellPhenotype","to_cellPhenotype"],chunksize=300000):
      for r in ch.itertuples(index=False):
        if int(r.from_is_epithelial)==1 or int(r.to_is_epithelial)==1: continue
        im=int(r.ImageNumber); a=int(r[1]); b=int(r[2])
        key=(im,min(a,b),max(a,b))
        if key in seen: continue
        seen.add(key)
        x=edge_counts.setdefault(im,[0,0]) # total non-epi edges, myo-endo edges
        x[0]+=1
        fa=str(r.from_cellPhenotype); ta=str(r.to_cellPhenotype)
        if ((fa in MYO and ta==ENDO) or (ta in MYO and fa==ENDO)): x[1]+=1

ed=pd.DataFrame([(im,v[0],v[1]) for im,v in edge_counts.items()],
                columns=["ImageNumber","n_non_epi_edges","n_myo_endo_edges"])
img=img.merge(ed,on="ImageNumber",how="left")
img[["n_non_epi_edges","n_myo_endo_edges"]]=img[["n_non_epi_edges","n_myo_endo_edges"]].fillna(0)
# non-epithelial graph universe includes endothelial
img["n_non_epi"]=img.n_universe_non_epi_nonendo+img.n_endo
img["expected_myo_endo_edges"]=np.where(img.n_non_epi>1,
    img.n_non_epi_edges*(2*img.n_myo*img.n_endo/(img.n_non_epi*(img.n_non_epi-1))),np.nan)

img.to_csv(OUT/"M17A1_image_level_spatial_components_v0_1.tsv.gz",sep="\t",index=False,compression="gzip")

# ---------- patient-condition aggregation ----------
sumcols=["n_cells","n_universe_non_epi_nonendo","n_vesselcontact_universe","n_myo","n_endo",
         "n_myo_vesselcontact","expected_myo_vesselcontact","n_non_epi_edges","n_myo_endo_edges",
         "expected_myo_endo_edges","sum_SMA_myo","n_SMA_myo","sum_PDGFRB_myo","n_PDGFRB_myo",
         "sum_CAV1_myo","n_CAV1_myo"]
p=img.groupby(["metabric_id","condition"],as_index=False)[sumcols].sum()
p["eligible_primary"]=p.n_myo.ge(MIN_MYO)&p.expected_myo_vesselcontact.gt(0)
p["eligible_graph"]=p.n_myo.ge(MIN_MYO)&p.n_endo.ge(MIN_ENDO)&p.expected_myo_endo_edges.gt(0)
p["myo_vesselcontact_fraction"]=p.n_myo_vesselcontact/np.where(p.n_myo>0,p.n_myo,np.nan)
p["vessel_contact_enrichment_log2"]=np.log2((p.n_myo_vesselcontact+PC)/(p.expected_myo_vesselcontact+PC))
p["spatial_decoupling_primary"]=-p.vessel_contact_enrichment_log2
p["myo_endo_edge_enrichment_log2"]=np.log2((p.n_myo_endo_edges+PC)/(p.expected_myo_endo_edges+PC))
p["spatial_decoupling_graph"]=-p.myo_endo_edge_enrichment_log2
p["mean_SMA_myo"]=p.sum_SMA_myo/p.n_SMA_myo.replace(0,np.nan)
p["mean_PDGFRB_myo"]=p.sum_PDGFRB_myo/p.n_PDGFRB_myo.replace(0,np.nan)
p["mean_CAV1_myo"]=p.sum_CAV1_myo/p.n_CAV1_myo.replace(0,np.nan)
p.to_csv(OUT/"M17A1_patient_spatial_locked_v0_1.tsv",sep="\t",index=False)

# ---------- outcome-blind condition tests ----------
metrics=[
 ("spatial_decoupling_primary","primary"),
 ("spatial_decoupling_graph","secondary"),
]
contrasts=[("Normal","DCIS"),("Normal","Tumor"),("DCIS","Tumor")]
tests=[]
for metric,level in metrics:
  elig="eligible_primary" if level=="primary" else "eligible_graph"
  q=p[p[elig]]
  for a,b in contrasts:
    x=q.loc[q.condition==a,metric].dropna().to_numpy(float)
    y=q.loc[q.condition==b,metric].dropna().to_numpy(float)
    if len(x)>=3 and len(y)>=3:
      U,pv=mannwhitneyu(x,y,alternative="two-sided")
      ma,mb=float(np.median(x)),float(np.median(y))
    else: U=pv=ma=mb=np.nan
    tests.append([metric,level,a,b,len(x),len(y),ma,mb,mb-ma,U,pv])
T=pd.DataFrame(tests,columns=["metric","level","group_a","group_b","n_a","n_b","median_a","median_b",
                              "delta_b_minus_a","U","p"])
T["FDR_all6"]=bh(T.p)
T.to_csv(OUT/"M17A1_condition_spatial_tests_v0_1.tsv",sep="\t",index=False)

# Concordance of two independent spatial definitions in invasive tumours
q=p[(p.condition=="Tumor")&p.eligible_primary&p.eligible_graph].dropna(
    subset=["spatial_decoupling_primary","spatial_decoupling_graph"])
rho,rhop=spearmanr(q.spatial_decoupling_primary,q.spatial_decoupling_graph) if len(q)>=5 else (np.nan,np.nan)

# Denominators
den=p.groupby("condition").agg(n_patient_condition=("metabric_id","size"),
  n_primary_eligible=("eligible_primary","sum"),n_graph_eligible=("eligible_graph","sum"),
  total_myo=("n_myo","sum"),total_endo=("n_endo","sum")).reset_index()
den.to_csv(OUT/"M17A1_spatial_denominators_v0_1.tsv",sep="\t",index=False)

# Semantic audit and gate
gate={
 "status":"M17A1_SPATIAL_PHENOTYPE_LOCKED",
 "no_clinical_outcome_table_read":True,
 "semantic_definition":{
   "is_perivascular":"cell is in contact with a vessel mask (Danenberg annotation), NOT a lineage",
   "mural_like_proxy":"deposited Myofibroblasts + Myofibroblasts PDPN+ phenotypes",
   "endothelial":"deposited Endothelial phenotype"
 },
 "primary_endpoint":"negative log2 enrichment of myofibroblast vessel-mask contact versus within-image expectation",
 "secondary_endpoint":"negative log2 enrichment of myofibroblast-endothelial graph edges versus within-image random-mixing expectation",
 "eligibility":{"min_myo":MIN_MYO,"min_endo_for_graph":MIN_ENDO},
 "tumour_crossmetric_n":int(len(q)),"tumour_crossmetric_spearman_rho":float(rho),
 "tumour_crossmetric_spearman_p":float(rhop),
 "condition_test_family":"six tests = 2 locked spatial metrics x 3 cross-sectional tissue contrasts",
 "interpretation_limit":"Mural-like spatial decoupling proxy; not direct pericyte lineage tracing or longitudinal detachment."
}
(AUD/"M17A1_spatial_lock_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
print("\nDENOMINATORS\n",den.to_string(index=False))
print("\nCONDITION TESTS\n",T.to_string(index=False))
