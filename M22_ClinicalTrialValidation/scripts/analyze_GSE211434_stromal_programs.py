import gzip,csv,json,math,hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1/StageM22_ClinicalTrialStromalValidation_v0_1")
ANN=ROOT/"raw/GSE211434_cell_types.txt.gz"
CNT=ROOT/"raw/GSE211434_counts.txt.gz"
OUT=ROOT/"derived"

PROGRAMS={
 "Fibroblasts":("fibroblast_matrix_contractile",["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","FN1","ACTA2","ITGA11","POSTN","JAG1"]),
 "Pericytes":("mural_matrix_contractile",["COL1A1","COL1A2","COL3A1","ACTA2","JAG1"]),
 "Endothelial cells":("endothelial_NOTCH_ECM",["JAG2","COL4A1","COL4A2","NOTCH4","ITGA1"])
}
TARGETS=sorted(set(g for _,genes in PROGRAMS.values() for g in genes))

def bh(pvals):
    p=np.asarray(pvals,float); n=len(p); q=np.full(n,np.nan)
    good=np.where(np.isfinite(p))[0]
    if not len(good): return q
    pg=p[good]; o=np.argsort(pg); adj=np.empty(len(good)); prev=1.0
    for k in range(len(good)-1,-1,-1):
        idx=o[k]; rank=k+1
        val=min(prev,pg[idx]*len(good)/rank); adj[idx]=val; prev=val
    q[good]=adj
    return q

# metadata
meta=pd.read_csv(ANN,sep="\t",dtype={"Day":int})
assert meta.CellID.is_unique
assert meta.Patient.nunique()==27
meta["group"]=meta.Patient.astype(str)+"|"+meta.Day.astype(str)+"|"+meta.CellType.astype(str)
group_info=(meta.groupby("group",sort=False)
    .agg(Patient=("Patient","first"),Day=("Day","first"),CellType=("CellType","first"),
         Treatment=("Treatment","first"),TreatmentArm=("TreatmentArm","first"),
         Response=("Response","first"),ResponseYN=("ResponseYN","first"),
         n_cells=("CellID","size"),library_size=("nCount_RNA","sum"))
    .reset_index())
groups=group_info["group"].tolist()
gindex={g:i for i,g in enumerate(groups)}
cell_group=np.array([gindex[g] for g in meta["group"]],dtype=np.int32)
cell_to_group=dict(zip(meta.CellID,cell_group))
ng=len(groups)

# count header and exact cell join
found={}
with gzip.open(CNT,"rt",errors="strict") as fh:
    header=fh.readline().rstrip("\n").split("\t")
    assert header[:3]==["Gene","Gene.ID","Gene.Symbol"],header[:3]
    cells=header[3:]
    assert len(cells)==len(meta),(len(cells),len(meta))
    assert len(set(cells))==len(cells)
    ann_ids=set(meta.CellID)
    count_ids=set(cells)
    assert ann_ids==count_ids,(len(ann_ids-count_ids),len(count_ids-ann_ids))
    codes=np.fromiter((cell_to_group[c] for c in cells),dtype=np.int32,count=len(cells))
    for line in fh:
        p1=line.find("\t"); p2=line.find("\t",p1+1); p3=line.find("\t",p2+1)
        sym=line[p2+1:p3]
        if sym not in TARGETS: continue
        vals=np.fromstring(line[p3+1:],sep="\t",dtype=np.float64)
        assert len(vals)==len(cells),(sym,len(vals),len(cells))
        agg=np.bincount(codes,weights=vals,minlength=ng)
        if sym in found:
            found[sym]=found[sym]+agg
        else:
            found[sym]=agg
        print("FOUND",sym,flush=True)

availability={g:(g in found) for g in TARGETS}
for lineage,(pname,genes) in PROGRAMS.items():
    present=[g for g in genes if g in found]
    assert len(present)/len(genes)>=0.70,(lineage,present)

# Pseudobulk gene table
pb=group_info.copy()
for gene in TARGETS:
    pb[gene+"_count"]=found.get(gene,np.full(ng,np.nan))
    if gene in found:
        pb[gene+"_log2CPM"]=np.log2(pb[gene+"_count"]/pb.library_size*1e6+1)

# program scores by threshold, z-scored within lineage across all eligible patient-timepoint PBs
all_scores=[]
all_tests=[]
loo=[]
for thr in [20,10]:
    for lineage,(pname,genes) in PROGRAMS.items():
        genes=[g for g in genes if g in found]
        d=pb[(pb.CellType==lineage)&(pb.n_cells>=thr)].copy()
        for g in genes:
            col=g+"_log2CPM"
            mu=d[col].mean(); sd=d[col].std(ddof=1)
            d[g+"_z"]=(d[col]-mu)/sd if sd>0 else 0.0
        d["program"]=pname
        d["threshold"]=thr
        d["program_score"]=d[[g+"_z" for g in genes]].mean(axis=1)
        all_scores.append(d)
        for day in [14,180]:
            wide=d[d.Day.isin([0,day])].pivot(index="Patient",columns="Day",values="program_score").dropna()
            patients=wide.index.tolist()
            a=wide[0].to_numpy(); b=wide[day].to_numpy(); delta=b-a
            p=wilcoxon(b,a,alternative="two-sided",zero_method="wilcox").pvalue if len(delta) else np.nan
            fullmed=float(np.median(delta)) if len(delta) else np.nan
            all_tests.append(dict(threshold=thr,lineage=lineage,program=pname,contrast=f"D0_D{day}",
                                  test_family="paired_change",n=len(delta),median_change=fullmed,p=float(p)))
            # leave-one-patient-out direction
            if len(delta)>=3:
                signs=[]
                for i,pt in enumerate(patients):
                    med=float(np.median(np.delete(delta,i)))
                    signs.append(np.sign(med))
                    loo.append(dict(threshold=thr,lineage=lineage,program=pname,contrast=f"D0_D{day}",
                                    omitted_patient=pt,loo_median_change=med,full_median_change=fullmed,
                                    same_direction=bool(np.sign(med)==np.sign(fullmed) or med==0 or fullmed==0)))
            # patient change metadata
            change=pd.DataFrame({"Patient":patients,"delta":delta})
            patient_meta=d.drop_duplicates("Patient").set_index("Patient")
            for family,field,A,B in [
                ("response_change","ResponseYN","Response","Non-response"),
                ("treatment_change","Treatment","Combination ribociclib","Letrozole alone")]:
                x=[row.delta for row in change.itertuples() if patient_meta.loc[row.Patient,field]==A]
                y=[row.delta for row in change.itertuples() if patient_meta.loc[row.Patient,field]==B]
                pu=mannwhitneyu(x,y,alternative="two-sided").pvalue if x and y else np.nan
                all_tests.append(dict(threshold=thr,lineage=lineage,program=pname,contrast=f"D0_D{day}",
                                      test_family=family,n=len(delta),n_A=len(x),n_B=len(y),group_A=A,group_B=B,
                                      median_change_A=float(np.median(x)) if x else np.nan,
                                      median_change_B=float(np.median(y)) if y else np.nan,
                                      difference_medians=float(np.median(x)-np.median(y)) if x and y else np.nan,
                                      p=float(pu) if np.isfinite(pu) else np.nan))

scores=pd.concat(all_scores,ignore_index=True)
tests=pd.DataFrame(all_tests)
tests["FDR"]=np.nan
for (thr,fam),idx in tests.groupby(["threshold","test_family"]).groups.items():
    ix=list(idx); tests.loc[ix,"FDR"]=bh(tests.loc[ix,"p"].to_numpy())
loo_df=pd.DataFrame(loo)
if len(loo_df):
    robust=(loo_df.groupby(["threshold","lineage","program","contrast"])["same_direction"].all()
            .rename("LOO_direction_robust").reset_index())
    tests=tests.merge(robust,on=["threshold","lineage","program","contrast"],how="left")

scores.to_csv(OUT/"GSE211434_stromal_program_patient_values.tsv.gz",sep="\t",index=False,compression="gzip")
tests.to_csv(OUT/"GSE211434_stromal_program_stats.tsv",sep="\t",index=False)
loo_df.to_csv(OUT/"GSE211434_stromal_program_LOO.tsv",sep="\t",index=False)
pd.DataFrame({"gene":TARGETS,"available":[availability[g] for g in TARGETS]}).to_csv(OUT/"GSE211434_program_gene_availability.tsv",sep="\t",index=False)

qc={
 "count_cells":len(cells),"annotation_cells":len(meta),"cell_id_set_exact":True,
 "target_genes_requested":len(TARGETS),"target_genes_found":sum(availability.values()),
 "target_gene_availability":availability,
 "program_scoreability":{lin:sum(g in found for g in genes)/len(genes) for lin,(_,genes) in PROGRAMS.items()}
}
(OUT/"M22_EXPRESSION_QC.json").write_text(json.dumps(qc,indent=2))
print("QC",json.dumps(qc),flush=True)
print("\nPRIMARY TESTS")
print(tests[tests.threshold.eq(20)].sort_values(["test_family","lineage","contrast"]).to_string(index=False),flush=True)
print("\nPRIMARY FDR<0.05")
print(tests[(tests.threshold.eq(20))&(tests.FDR<0.05)].to_string(index=False),flush=True)
