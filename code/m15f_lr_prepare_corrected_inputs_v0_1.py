#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
BASE=ROOT/"Stage07E_CCC_v0_1"
OUT=ST/"results/lr_corrected/stage07e"; OUT.mkdir(parents=True,exist_ok=True)

src=BASE/"WAVE_A_CCC_4COMP_LR_COUNTS_v0_1.h5ad"
a=ad.read_h5ad(src)
manifest=pd.read_csv(ROOT/"Stage06B/manifests/stage06a_sample_manifest.tsv",sep="\t",dtype=str).fillna("")
dc=manifest[(manifest.gse=="GSE195861")&(manifest.cell_population.str.upper()=="DCIS")][["gsm","title"]]
sample2uid=dict(zip(dc.gsm,"GSE195861:"+dc.title))
mask=a.obs.sample_id.astype(str).isin(sample2uid)
assert int(mask.sum())==14270,int(mask.sum())
assert a.obs.loc[mask,"role"].astype(str).eq("core_primary").all()
a.obs["role"]=a.obs["role"].astype(str)
a.obs["patient_uid"]=a.obs["patient_uid"].astype(str)
a.obs.loc[mask,"role"]="reference_dcis"
a.obs.loc[mask,"patient_uid"]=a.obs.loc[mask,"sample_id"].astype(str).map(sample2uid)
assert a.obs.loc[mask,"patient_uid"].nunique()==7
# Save corrected metadata/object for provenance and downstream deterministic rebuild.
corr=OUT/"WAVE_A_CCC_4COMP_LR_COUNTS_DCIscorrected_v0_1.h5ad"
a.write_h5ad(corr,compression="gzip")

totals=np.load(BASE/"full_library_total_counts_4comp_v0_1.npy")
assert len(totals)==a.n_obs
X=a.X.tocsr();obs=a.obs.copy();genes=a.var_names.astype(str).to_numpy()

# Patient x broad compartment pseudobulk.
meta=[];mats=[]
for (role,pid,ct),idx in obs.groupby(["role","patient_uid","cell_type"],observed=True).indices.items():
 idx=np.asarray(idx,dtype=int);lib=float(totals[idx].sum());counts=np.asarray(X[idx,:].sum(axis=0)).ravel().astype(float)
 cpm=(counts/lib*1e6 if lib>0 else np.zeros_like(counts)).astype(np.float32)
 gse=str(obs.iloc[idx[0]].gse);meta.append([role,pid,gse,ct,len(idx),lib]);mats.append(cpm)
M=np.vstack(mats).astype(np.float32)
md=pd.DataFrame(meta,columns=["role","patient_uid","gse","cell_type","n_cells","library_total_counts"])
md.to_csv(OUT/"patient_compartment_pseudobulk_metadata_v0_1.tsv",sep="\t",index=False)
np.savez_compressed(OUT/"patient_compartment_pseudobulk_CPM_v0_1.npz",X=M,genes=genes)

# Patient x reference-state pseudobulk over same LR gene matrix.
gcols=["role","patient_uid","gse","cell_type","state_affinity"]
groups=obs[gcols].drop_duplicates().reset_index(drop=True);groups["group_id"]=np.arange(len(groups),dtype=int)
mm=obs.reset_index().merge(groups,on=gcols,how="left",validate="many_to_one")
gid=mm.group_id.to_numpy();H=sparse.csr_matrix((np.ones(len(gid)),(gid,np.arange(len(gid)))),shape=(len(groups),len(gid)))
agg=H@X;lib=np.asarray(H@totals.reshape(-1,1)).ravel();cells=np.bincount(gid,minlength=len(groups))
groups["n_cells"]=cells;groups["library_total_counts"]=lib
groups.to_csv(OUT/"state_pseudobulk_metadata_v0_1.tsv",sep="\t",index=False)
np.savez_compressed(OUT/"state_pseudobulk_counts_v0_1.npz",X=agg.toarray().astype(np.float32),genes=genes)

inv=(obs.groupby(["role","gse","cell_type"],observed=True).size().rename("n_cells").reset_index())
inv.to_csv(OUT/"corrected_role_cell_inventory_v0_1.tsv",sep="\t",index=False)
audit={
 "status":"M15F_LR_INPUT_CORRECTION_PASS",
 "total_cells":int(a.n_obs),
 "relabelled_DCIScells":int(mask.sum()),
 "DCIS_patients":int(a.obs.loc[mask,"patient_uid"].nunique()),
 "DCIS_by_compartment":a.obs.loc[mask,"cell_type"].value_counts().astype(int).to_dict(),
 "GSE195861_core_primary_patients_after":int(a.obs[(a.obs.gse.astype(str)=="GSE195861")&(a.obs.role.astype(str)=="core_primary")].patient_uid.nunique()),
 "GSE195861_DCISpatients_after":int(a.obs[(a.obs.gse.astype(str)=="GSE195861")&(a.obs.role.astype(str)=="reference_dcis")].patient_uid.nunique())
}
(ST/"audit/M15F_LR_input_correction_gate_v0_1.json").write_text(json.dumps(audit,indent=2)+"\n")
print(json.dumps(audit,indent=2))
