#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd,json

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
SRC=ST/"source_data/supplementary"
RES=ST/"results/s4_mapper"

bp=pd.read_csv(RES/"M11A3_reproduced_patient_confidence_v0_1.tsv",sep="\t")
cal=pd.read_csv(SRC/"FigureS4_degradation_calibration_v0_1.tsv",sep="\t")
null=pd.read_csv(SRC/"FigureS4_full_state_degradation_null_v0_1.tsv",sep="\t")
frac=pd.read_csv(RES/"M15F_reproduced_patient_state_fractions_DCIscorrected_v0_1.tsv",sep="\t")

# Confidence calibration: thinning outputs are frozen; only observed invasive denominator changes.
for lin in cal.lineage.unique():
    q=bp[bp.lineage==lin]
    n=q[q.stage06a_role=="reference_normal"]
    p=q[q.stage06a_role=="core_primary"]
    nmc=float(n.median_max_corr.median()); pmc=float(p.median_max_corr.median())
    nma=float(n.median_margin.median()); pma=float(p.median_margin.median())
    idx=cal.lineage==lin
    # verify baseline normal is exactly reproducible to numerical tolerance
    assert np.allclose(cal.loc[idx,"baseline_normal_median_max_corr"],nmc,rtol=0,atol=1e-7)
    assert np.allclose(cal.loc[idx,"baseline_normal_median_margin"],nma,rtol=0,atol=1e-7)
    cal.loc[idx,"primary_median_max_corr"]=pmc
    cal.loc[idx,"primary_median_margin"]=pma
    cal.loc[idx,"fraction_observed_maxcorr_drop_reproduced"]=(nmc-cal.loc[idx,"degraded_normal_median_max_corr"])/(nmc-pmc)
    cal.loc[idx,"fraction_observed_margin_drop_reproduced"]=(nma-cal.loc[idx,"degraded_normal_median_margin"])/(nma-pma)
cal.to_csv(SRC/"FigureS4_degradation_calibration_v0_1.tsv",sep="\t",index=False)

# Build complete patient fractions from mapper-reproduced assignments, adding zero fractions.
states={"Fibroblast":["FB1","FB2","FB3","FB4"],"Endothelial":["VEV","VEA","VEAT","VEC","LE1","LE2"],"Perivascular":["PV1","PV2","PV3","PV4","PV5"]}
keys=frac[["lineage","stage06a_role","patient_uid"]].drop_duplicates()
wide=frac.pivot_table(index=["lineage","stage06a_role","patient_uid"],columns="repro_state",values="fraction",aggfunc="first",fill_value=0).reset_index()
for lin,ss in states.items():
    for s in ss:
        if s not in wide:wide[s]=0.0

for i,r in null.iterrows():
    lin=r.lineage;st=r.state
    q=wide[wide.lineage==lin]
    bn=float(q[q.stage06a_role=="reference_normal"][st].median())
    pm=float(q[q.stage06a_role=="core_primary"][st].median())
    # same normal mapper baseline as historical degradation experiment
    if not np.isclose(bn,float(r.baseline_normal_median_fraction),atol=1e-7,rtol=0):
        raise RuntimeError((st,"normal baseline mismatch",bn,r.baseline_normal_median_fraction))
    obs=pm-bn
    nul=float(r.degraded_normal_median_fraction)-bn
    ratio=nul/obs if abs(obs)>1e-12 else np.nan
    null.loc[i,"primary_median_fraction"]=pm
    null.loc[i,"observed_primary_minus_normal"]=obs
    null.loc[i,"null_degraded_minus_normal"]=nul
    null.loc[i,"null_to_observed_ratio"]=ratio
    null.loc[i,"direction_concordant"]=bool(obs!=0 and nul!=0 and np.sign(obs)==np.sign(nul))
    null.loc[i,"abs_null_to_observed_ratio"]=abs(ratio) if np.isfinite(ratio) else np.nan
    if not np.isfinite(ratio):
        interp="observed shift ~0; ratio undefined"
    elif abs(ratio)<0.1:
        interp="null small (<10% observed)"
    elif abs(ratio)<1:
        interp="null partial (10-100% observed)"
    else:
        interp="null larger than observed; abundance direction unstable"
    null.loc[i,"interpretation"]=interp
null.to_csv(SRC/"FigureS4_full_state_degradation_null_v0_1.tsv",sep="\t",index=False)

# Source-level patient confidence is already role-corrected; verify no DCIS remains core_primary.
pconf=pd.read_csv(SRC/"FigureS4_patient_confidence_v0_1.tsv",sep="\t")
dcis={"GSE195861:NCCBC2","GSE195861:NCCBC3","GSE195861:NCCBC5","GSE195861:NCCBC6","GSE195861:NCCBC7","GSE195861:NCCBC11","GSE195861:NCCBC14"}
assert not ((pconf.patient_uid.isin(dcis))&(pconf.stage06a_role=="core_primary")).any()

q=cal[cal.condition=="tumor_depth_matched"][["lineage","primary_median_max_corr","fraction_observed_maxcorr_drop_reproduced","primary_median_margin","fraction_observed_margin_drop_reproduced","mean_label_stability"]]
gate={"status":"M15F_S4_DCIscorrected_PASS",
      "invasive_primary_patients_reproduced_mapper":{lin:int(((bp.lineage==lin)&(bp.stage06a_role=="core_primary")).sum()) for lin in states},
      "tumor_depth_matched":q.to_dict("records"),
      "VEAT_null_ratio":float(null.loc[null.state=="VEAT","null_to_observed_ratio"].iloc[0])}
(ST/"audit/M15F_S4_calibration_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
