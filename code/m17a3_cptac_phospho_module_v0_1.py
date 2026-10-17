#!/home/pankaj/venvs/bc_scanpy/bin/python
"""M17A3: predefined CPTAC paired phosphoproteomic module with total-protein adjustment.

Interpretation is deliberately site-agnostic with respect to activation/inhibition.
A positive adjusted delta means phosphosite signal rises more than the corresponding
total-protein signal; it is NOT automatically called activating phosphorylation.
"""
from pathlib import Path
import json,re,warnings
import numpy as np,pandas as pd
from scipy.stats import wilcoxon

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM17A_VascularSpatialClinicalPhospho_v0_1"
OUT=ST/"results"; AUD=ST/"audit"
PDC=ROOT/"Stage11E_Proteomics_v0_1"
RAW=PDC/"raw/pdc/processed"

MODULES={
 "ECM_receptor_context":["COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","ITGA1","ITGB1","PDGFRB"],
 "focal_adhesion":["PTK2","PXN","VCL","TLN1","SRC"],
 "Rho_actomyosin":["RHOA","ROCK1","ROCK2","MYL9","MYH9","ACTA2"],
 "cofilin_axis":["LIMK1","LIMK2","CFL1"],
 "YAP_TAZ_mechanotransduction":["YAP1","WWTR1","TEAD3"],
}
GENES=list(dict.fromkeys(g for gs in MODULES.values() for g in gs))
GENE_MODULE={g:m for m,gs in MODULES.items() for g in gs}
MINPAIRS=5

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    v=p[ok]; o=np.argsort(v); q=v[o]*len(v)/np.arange(1,len(v)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(v)); z[o]=np.minimum(q,1); out[ok]=z
    return out

def load_bio(sid):
    j=json.load(open(PDC/"raw/pdc"/f"{sid}_biospecimen.json"))
    return pd.DataFrame(j["data"]["biospecimenPerStudy"])

def header_id(col,suffix):
    x=col[:-len(suffix)]
    return re.sub(r"\.\d+$","",x)

def map_cols(cols,bio,suffix):
    exact={str(x):r for x,r in bio.set_index("aliquot_submitter_id",drop=False).iterrows()}
    out=[]
    for c in cols:
        x=header_id(c,suffix); r=exact.get(x)
        if r is None: out.append((c,None,None,None,None))
        else: out.append((c,str(r["case_submitter_id"]),str(r["sample_type"]),
                           str(r["aliquot_submitter_id"]),str(r.get("pool"))))
    return pd.DataFrame(out,columns=["column","case_id","sample_type","aliquot_id","pool"])

def collapse_case(X,cmap):
    valid=cmap.case_id.notna() & ~cmap["pool"].astype(str).str.lower().eq("yes") & ~cmap.sample_type.eq("Not Reported")
    cm=cmap.loc[valid].copy(); cols=cm.index.to_numpy()
    keys=cm[["case_id","sample_type"]].astype(str).agg("||".join,axis=1)
    uniq=list(dict.fromkeys(keys))
    arr=np.full((X.shape[0],len(uniq)),np.nan,np.float32); meta=[]
    for j,k in enumerate(uniq):
        ix=np.where(keys.to_numpy()==k)[0]; orig=cols[ix]
        arr[:,j]=np.nanmedian(X[:,orig],axis=1)
        case,typ=k.split("||",1); meta.append((case,typ,len(orig)))
    return arr,pd.DataFrame(meta,columns=["case_id","sample_type","n_technical_columns"])

def load_proteome():
    d=pd.read_csv(RAW/"CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv",sep="\t",low_memory=False)
    d=d[~d["Gene"].isin(["Mean","Median","StdDev"])].copy()
    cols=[c for c in d.columns if c.endswith(" Unshared Log Ratio")]
    cmap=map_cols(cols,load_bio("PDC000120")," Unshared Log Ratio")
    X=d[cols].apply(pd.to_numeric,errors="coerce").to_numpy(np.float32)
    arr,meta=collapse_case(X,cmap)
    return d["Gene"].astype(str).to_numpy(),arr,meta

def load_phospho():
    d=pd.read_csv(RAW/"CPTAC2_Breast_Prospective_Collection_BI_Phosphoproteome.phosphosite.tmt10.tsv",
                  sep="\t",low_memory=False)
    cols=[c for c in d.columns if c.endswith(" Log Ratio")]
    cmap=map_cols(cols,load_bio("PDC000121")," Log Ratio")
    X=d[cols].apply(pd.to_numeric,errors="coerce").to_numpy(np.float32)
    arr,meta=collapse_case(X,cmap)
    return d["Phosphosite"].astype(str).to_numpy(),d["Gene"].astype(str).to_numpy(),arr,meta

def case_indices(meta):
    t={str(r.case_id):int(i) for i,r in meta[meta.sample_type.eq("Primary Tumor")].iterrows()}
    n={str(r.case_id):int(i) for i,r in meta[meta.sample_type.isin(["Solid Tissue Normal","Normal Adjacent Tissue"])].iterrows()}
    return t,n

pg,prot,pmeta=load_proteome()
site,sg,pho,smeta=load_phospho()
pt,pn=case_indices(pmeta); st,sn=case_indices(smeta)
paired=sorted(set(pt)&set(pn)&set(st)&set(sn))
assert len(paired)==18, len(paired)

# Gene-level total protein delta per case, median across duplicate rows if any.
prot_delta={}
for g in GENES:
    ix=np.where(pg==g)[0]
    if not len(ix): continue
    vals={}
    for case in paired:
        tv=prot[ix,pt[case]]; nv=prot[ix,pn[case]]
        dd=tv-nv; dd=dd[np.isfinite(dd)]
        if len(dd): vals[case]=float(np.median(dd))
    prot_delta[g]=vals

global_sites=pd.read_csv(PDC/"results/CPTAC_Prospective_phosphosite_tumor_normal_paired_v0_1.tsv.gz",sep="\t")
global_prot=pd.read_csv(PDC/"results/CPTAC_Prospective_proteome_tumor_normal_paired_v0_1.tsv.gz",sep="\t")
gsite=global_sites.set_index("feature_id",drop=False)

rows=[]; case_rows=[]
for i,(fid,g) in enumerate(zip(site,sg)):
    if g not in GENES: continue
    raw=[]; adj=[]; cases_raw=[]; cases_adj=[]
    for case in paired:
        a=pho[i,st[case]]; b=pho[i,sn[case]]
        if np.isfinite(a) and np.isfinite(b):
            dv=float(a-b); raw.append(dv); cases_raw.append(case)
            pdlt=prot_delta.get(g,{}).get(case,np.nan)
            if np.isfinite(pdlt):
                adj.append(dv-pdlt); cases_adj.append(case)
                case_rows.append([case,GENE_MODULE[g],g,fid,dv,pdlt,dv-pdlt])
    pr=pa=np.nan
    if len(raw)>=MINPAIRS and np.any(np.abs(raw)>0):
        try: pr=float(wilcoxon(raw,alternative="two-sided",zero_method="wilcox").pvalue)
        except Exception: pass
    if len(adj)>=MINPAIRS and np.any(np.abs(adj)>0):
        try: pa=float(wilcoxon(adj,alternative="two-sided",zero_method="wilcox").pvalue)
        except Exception: pass
    old=gsite.loc[fid] if fid in gsite.index else None
    rows.append(dict(module=GENE_MODULE[g],gene=g,phosphosite=fid,
        n_raw=len(raw),raw_median_delta=float(np.median(raw)) if raw else np.nan,raw_p=pr,
        n_protein_adjusted=len(adj),protein_adjusted_median_delta=float(np.median(adj)) if adj else np.nan,
        protein_adjusted_p=pa,
        global_phosphosite_FDR=np.nan if old is None else float(old.paired_FDR) if np.isfinite(old.paired_FDR) else np.nan,
        global_raw_median_delta=np.nan if old is None else float(old.paired_median_delta) if np.isfinite(old.paired_median_delta) else np.nan))
R=pd.DataFrame(rows)
R["raw_module_family_FDR"]=bh(R.raw_p)
R["protein_adjusted_module_family_FDR"]=bh(R.protein_adjusted_p)
R.to_csv(OUT/"M17A3_CPTAC_module_phosphosite_tests_v0_1.tsv",sep="\t",index=False)
pd.DataFrame(case_rows,columns=["case_id","module","gene","phosphosite","site_delta","protein_delta",
                                "site_minus_protein_delta"]).to_csv(
    OUT/"M17A3_CPTAC_case_level_adjusted_deltas_v0_1.tsv.gz",sep="\t",index=False,compression="gzip")

# Gene availability and summary.
avail=[]
for g in GENES:
    q=R[R.gene.eq(g)]
    gp=global_prot[global_prot.gene.eq(g)]
    avail.append(dict(module=GENE_MODULE[g],gene=g,
       phosphosites_detected=len(q),
       phosphosites_raw_testable=int(q.n_raw.ge(MINPAIRS).sum()),
       phosphosites_global_FDR05=int(q.global_phosphosite_FDR.lt(.05).sum()),
       phosphosites_adjusted_testable=int(q.n_protein_adjusted.ge(MINPAIRS).sum()),
       phosphosites_adjusted_module_FDR05=int(q.protein_adjusted_module_family_FDR.lt(.05).sum()),
       protein_n_paired=int(gp.n_paired.max()) if len(gp) and np.isfinite(gp.n_paired.max()) else 0,
       protein_paired_median_delta=float(gp.paired_median_delta.iloc[0]) if len(gp) and np.isfinite(gp.paired_median_delta.iloc[0]) else np.nan,
       protein_global_FDR=float(gp.paired_FDR.iloc[0]) if len(gp) and np.isfinite(gp.paired_FDR.iloc[0]) else np.nan))
A=pd.DataFrame(avail)
A.to_csv(OUT/"M17A3_CPTAC_module_gene_summary_v0_1.tsv",sep="\t",index=False)

# Module summary is descriptive; sites are correlated and have heterogeneous regulatory meaning.
MS=A.groupby("module").agg(n_genes=("gene","size"),
    genes_with_testable_adjusted_site=("phosphosites_adjusted_testable",lambda x:int((x>0).sum())),
    genes_with_global_FDR05_site=("phosphosites_global_FDR05",lambda x:int((x>0).sum())),
    genes_with_adjusted_module_FDR05_site=("phosphosites_adjusted_module_FDR05",lambda x:int((x>0).sum()))
    ).reset_index()
MS.to_csv(OUT/"M17A3_CPTAC_module_summary_v0_1.tsv",sep="\t",index=False)

# Highlight sites only by prespecified module + statistics, never by functional "activation".
H=R[(R.n_protein_adjusted>=MINPAIRS)].sort_values(
    ["protein_adjusted_module_family_FDR","global_phosphosite_FDR","gene","phosphosite"])
H.to_csv(OUT/"M17A3_CPTAC_testable_adjusted_sites_ranked_v0_1.tsv",sep="\t",index=False)

n_test=int(R.n_protein_adjusted.ge(MINPAIRS).sum())
n_sig=int(R.protein_adjusted_module_family_FDR.lt(.05).sum())
n_genes_sig=int(R.loc[R.protein_adjusted_module_family_FDR.lt(.05),"gene"].nunique())
gate={
 "status":"M17A3_CPTAC_PHOSPHO_COMPLETE",
 "paired_cases":len(paired),
 "predefined_modules":MODULES,
 "protein_adjustment":"case-matched phosphosite tumour-normal delta minus case-matched total-protein tumour-normal delta for same gene",
 "testable_adjusted_sites":n_test,
 "adjusted_sites_FDR05_within_predefined_module_family":n_sig,
 "genes_with_adjusted_site_FDR05":n_genes_sig,
 "global_FDR_reference":"Stage11E BH across all 27,755 prospective phosphosite rows retained alongside module-family adjusted FDR",
 "interpretation_lock":"Significant site remodeling is not labelled activating/inhibitory without independent site-specific functional evidence.",
}
(AUD/"M17A3_phospho_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
print("\nMODULE SUMMARY\n",MS.to_string(index=False))
print("\nTOP ADJUSTED SITES\n",H.head(30).to_string(index=False))
