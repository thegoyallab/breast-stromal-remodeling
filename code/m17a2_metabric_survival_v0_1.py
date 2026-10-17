#!/home/pankaj/venvs/bc_scanpy/bin/python
"""M17A2: prespecified METABRIC survival analysis using locked M17A1 spatial scores.

Primary score: M17A1 vessel-mask spatial decoupling score.
Secondary scores: graph decoupling; predeclared contractile mural-like score
  mean[z(log1p(SMA)), z(log1p(PDGFRB)), -z(log1p(Caveolin-1))]
computed in the deposited myofibroblast phenotypes.

Discovery/validation split is the Danenberg isValidation field.
Score scaling parameters are estimated in discovery and applied unchanged to validation.
Outcome is breast-cancer-specific status/time from the Danenberg clinical table.
"""
from pathlib import Path
import json,warnings
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from statsmodels.duration.hazard_regression import PHReg

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM17A_VascularSpatialClinicalPhospho_v0_1"
OUT=ST/"results"; AUD=ST/"audit"
P=OUT/"M17A1_patient_spatial_locked_v0_1.tsv"
C=ROOT/"Parallel_Modalities/IMC_METABRIC_Danenberg/audit_temp/IMCClinical_export.csv"

p=pd.read_csv(P,sep="\t")
p=p[p.condition.eq("Tumor")].copy()
p["idnorm"]=p.metabric_id.astype(str).str.replace("-","",regex=False)
c=pd.read_csv(C)
c["idnorm"]=c.metabric_id.astype(str).str.replace("-","",regex=False)
d=p.merge(c,on="idnorm",how="inner",suffixes=("","_clinical"))
assert d.idnorm.nunique()==len(d)
d["cohort"]=np.where(d.isValidation,"validation","discovery")

# Clinical coding fixed before outcome modeling.
size_map={"1-5mm":1,"5-10mm":2,"10-20mm":3,"20-50mm":4,"50+mm":5}
node_map={"0":0,"1":1,"2-3":2,"4-7":3,"7+":4}
d["size_ord"]=d.sizeOrdinal.map(size_map)
d["nodes_ord"]=d.LymphNodesOrdinal.astype(str).map(node_map)
d["ER_pos"]=d.ERStatus.map({"neg":0,"pos":1})
d["HER2_pos"]=d.ERBB2_pos.map({False:0,True:1})
d["Grade_num"]=pd.to_numeric(d.Grade,errors="coerce")
d["time"]=pd.to_numeric(d.yearsToStatus,errors="coerce")
d["event"]=pd.to_numeric(d.DeathBreast,errors="coerce")

# Lock score scaling on discovery, without outcome dependence.
disc=d[d.cohort.eq("discovery")].copy()
params={}
def z_apply(col,newcol,eligible=None):
    q=disc[col].dropna() if eligible is None else disc.loc[eligible,col].dropna()
    mu=float(q.mean()); sd=float(q.std(ddof=1))
    params[newcol]={"source":col,"mean_discovery":mu,"sd_discovery":sd}
    d[newcol]=(d[col]-mu)/sd

z_apply("spatial_decoupling_primary","z_spatial_primary",disc.eligible_primary)
z_apply("spatial_decoupling_graph","z_spatial_graph",disc.eligible_graph)

# Contractile mural-like score using log1p marker means, standardized in discovery.
for marker in ["mean_SMA_myo","mean_PDGFRB_myo","mean_CAV1_myo"]:
    logcol="log1p_"+marker
    d[logcol]=np.log1p(pd.to_numeric(d[marker],errors="coerce"))
    q=d.loc[d.cohort.eq("discovery")&d.eligible_primary,logcol].dropna()
    mu=float(q.mean()); sd=float(q.std(ddof=1))
    zcol="z_"+marker
    d[zcol]=(d[logcol]-mu)/sd
    params[zcol]={"source":logcol,"mean_discovery":mu,"sd_discovery":sd}
d["contractile_score"]=(d["z_mean_SMA_myo"]+d["z_mean_PDGFRB_myo"]-d["z_mean_CAV1_myo"])/3.0

# eligibility for score models
d["score_primary_ok"]=d.eligible_primary & d.z_spatial_primary.notna()
d["score_graph_ok"]=d.eligible_graph & d.z_spatial_graph.notna()
d["score_contractile_ok"]=d.eligible_primary & d.contractile_score.notna()

(OUT/"M17A2_score_scaling_parameters_v0_1.json").write_text(json.dumps(params,indent=2)+"\n")
d.to_csv(OUT/"M17A2_METABRIC_locked_scores_clinical_v0_1.tsv",sep="\t",index=False)

def fit_cox(frame,score,covars,model_name,cohort):
    cols=["time","event",score]+covars
    q=frame[cols].dropna().copy()
    if len(q)<25 or q.event.sum()<10:
        return None,None
    X=q[[score]+covars].astype(float)
    # Drop zero-variance covariates deterministically.
    keep=[x for x in X.columns if X[x].std(ddof=0)>1e-12]
    X=X[keep]
    try:
        mod=PHReg(q.time.to_numpy(float),X.to_numpy(float),status=q.event.to_numpy(int),ties="efron")
        res=mod.fit(disp=0)
    except Exception as e:
        return {"cohort":cohort,"model":model_name,"score":score,"n":len(q),"events":int(q.event.sum()),
                "status":"FIT_FAILED","error":str(e)},None
    names=keep
    rows=[]
    for j,name in enumerate(names):
        beta=float(res.params[j]); se=float(res.bse[j])
        z=beta/se if se>0 else np.nan
        # statsmodels exposes pvalues
        pv=float(res.pvalues[j])
        rows.append(dict(cohort=cohort,model=model_name,score=score,term=name,n=len(q),
                         events=int(q.event.sum()),beta=beta,SE=se,HR=float(np.exp(beta)),
                         CI95_low=float(np.exp(beta-1.96*se)),CI95_high=float(np.exp(beta+1.96*se)),
                         p=pv,status="OK"))
    # lightweight PH diagnostic: Spearman Schoenfeld residual vs log(event time)
    ph=[]
    try:
        sr=np.asarray(res.schoenfeld_residuals)
        ev=q.event.to_numpy(int)==1
        lt=np.log(q.time.to_numpy(float)[ev])
        ser=sr[ev,:] if sr.shape[0]==len(q) else sr
        for j,name in enumerate(names):
            vv=ser[:,j]
            ok=np.isfinite(vv)&np.isfinite(lt)
            rho,pv=spearmanr(vv[ok],lt[ok]) if ok.sum()>=10 else (np.nan,np.nan)
            ph.append(dict(cohort=cohort,model=model_name,term=name,n_events_diag=int(ok.sum()),
                           schoenfeld_spearman_rho=float(rho),schoenfeld_p=float(pv)))
    except Exception:
        pass
    return rows,ph

# Predeclared model hierarchy.
scores=[
 ("z_spatial_primary","score_primary_ok","primary_spatial"),
 ("z_spatial_graph","score_graph_ok","secondary_graph"),
 ("contractile_score","score_contractile_ok","secondary_contractile"),
]
models=[
 ("univariable",[]),
 ("adjusted_core",["Grade_num","size_ord","nodes_ord","ER_pos"]),
 ("adjusted_core_HER2",["Grade_num","size_ord","nodes_ord","ER_pos","HER2_pos"]),
]
rows=[]; phrows=[]
for cohort in ["discovery","validation"]:
  sub=d[d.cohort.eq(cohort)]
  for score,elig,label in scores:
    ss=sub[sub[elig]].copy()
    for mn,covs in models:
      rr,pp=fit_cox(ss,score,covs,mn,cohort)
      if rr is None: continue
      if isinstance(rr,dict): rows.append(rr)
      else: rows.extend(rr)
      if pp: phrows.extend(pp)

R=pd.DataFrame(rows)
R.to_csv(OUT/"M17A2_METABRIC_survival_models_v0_1.tsv",sep="\t",index=False)
PH=pd.DataFrame(phrows)
if len(PH): PH.to_csv(OUT/"M17A2_METABRIC_PH_diagnostics_v0_1.tsv",sep="\t",index=False)

# Primary gate uses score coefficient from core-adjusted model, same locked score and direction.
def score_row(cohort,model,score):
    q=R[(R.cohort==cohort)&(R.model==model)&(R.term==score)&(R.status=="OK")]
    return q.iloc[0] if len(q) else None
a=score_row("discovery","adjusted_core","z_spatial_primary")
b=score_row("validation","adjusted_core","z_spatial_primary")
if a is not None and b is not None:
    same_direction=bool(np.sign(a.beta)==np.sign(b.beta) and a.beta!=0 and b.beta!=0)
    validation_supported=bool(b.p<0.05 and same_direction)
else:
    same_direction=False; validation_supported=False

# Secondary validation multiplicity across three locked scores for adjusted_core.
val=R[(R.cohort=="validation")&(R.model=="adjusted_core")&
      (R.term.isin(["z_spatial_primary","z_spatial_graph","contractile_score"]))].copy()
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=p[o]*len(p)/np.arange(1,len(p)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(p)); z[o]=np.minimum(q,1); return z
if len(val): val["FDR_three_scores"]=bh(val.p.to_numpy(float))
val.to_csv(OUT/"M17A2_validation_score_summary_v0_1.tsv",sep="\t",index=False)

den=d.groupby(["cohort"]).agg(n_clinical=("idnorm","size"),events=("event","sum"),
    n_primary=("score_primary_ok","sum"),n_graph=("score_graph_ok","sum"),
    n_contractile=("score_contractile_ok","sum")).reset_index()
den.to_csv(OUT/"M17A2_survival_denominators_v0_1.tsv",sep="\t",index=False)

gate={
 "status":"M17A2_SURVIVAL_COMPLETE",
 "clinical_overlap_tumour_patients":int(d.idnorm.nunique()),
 "discovery_validation_split_source":"Danenberg isValidation",
 "primary_endpoint":"breast-cancer-specific survival; continuous 1-SD M17A1 spatial decoupling score",
 "primary_adjustment":["Grade","tumour size ordinal","lymph-node burden ordinal","ER status"],
 "HER2_model":"prespecified sensitivity due missing HER2 annotations",
 "primary_discovery_HR":None if a is None else float(a.HR),
 "primary_discovery_p":None if a is None else float(a.p),
 "primary_validation_HR":None if b is None else float(b.HR),
 "primary_validation_p":None if b is None else float(b.p),
 "same_effect_direction":same_direction,
 "validation_gate_pass":validation_supported,
 "interpretation_if_fail":"No prognostic claim; retain survival analysis as negative/exploratory audit only."
}
(AUD/"M17A2_survival_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
print("\nDENOMINATORS\n",den.to_string(index=False))
print("\nVALIDATION SCORE SUMMARY\n",val.to_string(index=False))
