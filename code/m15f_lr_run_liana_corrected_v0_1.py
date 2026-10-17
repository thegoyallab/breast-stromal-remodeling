#!/home/pankaj/venvs/bc_scanpy/bin/python
import os,sys,json,warnings
os.environ["OMP_NUM_THREADS"]="20";os.environ["OPENBLAS_NUM_THREADS"]="20";os.environ["MKL_NUM_THREADS"]="20";os.environ["NUMBA_NUM_THREADS"]="20"
from pathlib import Path
ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
BASE=ROOT/"Stage07E_CCC_v0_1"
OUT=ST/"results/lr_corrected/stage07e"
sys.path.insert(0,str(BASE/"vendor"))
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import liana as li
warnings.filterwarnings("ignore")
rng=np.random.default_rng(20261009)
a=ad.read_h5ad(OUT/"WAVE_A_CCC_4COMP_LR_COUNTS_DCIscorrected_v0_1.h5ad")
totals=np.load(BASE/"full_library_total_counts_4comp_v0_1.npy")
res=pd.read_csv(BASE/"LIANA_consensus_resource_v1_10_0.tsv",sep="\t").drop_duplicates()
sel=[];rows=[]
# reference_dcis deliberately not used to construct invasive-vs-normal discovery scores
for (role,pid,ct),idx in a.obs.groupby(["role","patient_uid","cell_type"],observed=True).indices.items():
    if str(role)=="reference_dcis": continue
    idx=np.asarray(idx,dtype=int);n=len(idx)
    if n<10:continue
    take=idx if n<=200 else np.sort(rng.choice(idx,200,replace=False))
    sel.extend(take.tolist());rows.append([role,pid,ct,n,len(take)])
sel=np.array(sorted(sel),dtype=int)
inv=pd.DataFrame(rows,columns=["role","patient_uid","cell_type","available_cells","sampled_cells"])
inv.to_csv(OUT/"liana_balanced_sampling_inventory_v0_1.tsv",sep="\t",index=False)
X=a.X[sel,:].tocsr().astype(float);scale=1e4/np.maximum(totals[sel],1.0);X=sparse.diags(scale)@X;X.data=np.log1p(X.data);X=X.astype(np.float32)
b=ad.AnnData(X=X,obs=a.obs.iloc[sel].copy(),var=a.var.copy())
outs=[]
for role in ["reference_normal","reference_preneoplastic","core_primary","linked_progression_LN"]:
 sub=b[b.obs["role"].astype(str).to_numpy()==role].copy()
 print("ROLE",role,"cells",sub.n_obs,"patients",sub.obs.patient_uid.nunique(),"counts",sub.obs.cell_type.value_counts().to_dict(),flush=True)
 rr=li.mt.rank_aggregate(sub,groupby="cell_type",resource=res,use_raw=False,expr_prop=.10,min_cells=10,n_perms=100,n_jobs=20,return_all_lrs=True,inplace=False,verbose=False)
 rr["role"]=role;rr.to_csv(OUT/f"liana_consensus_{role}_v0_1.tsv.gz",sep="\t",index=False,compression="gzip");outs.append(rr)
allr=pd.concat(outs,ignore_index=True);allr.to_csv(OUT/"liana_consensus_all_roles_v0_1.tsv.gz",sep="\t",index=False,compression="gzip")
summ=inv.groupby(["role","cell_type"],observed=True).agg(patients=("patient_uid","nunique"),cells=("sampled_cells","sum")).reset_index()
summ.to_csv(OUT/"liana_balanced_role_summary_v0_1.tsv",sep="\t",index=False)
manifest={"liana_version":li.__version__,"resource_pairs_input":int(len(res)),"sampled_cells":int(len(sel)),"min_cells_patient_compartment":10,"cap_cells_patient_compartment":200,"expr_prop":.10,"n_perms":100,"n_jobs":20,"result_rows":int(len(allr)),"DCIS_excluded_from_invasive_role":True}
(OUT/"liana_manifest_v0_1.json").write_text(json.dumps(manifest,indent=2)+"\n")
print(json.dumps(manifest,indent=2),flush=True)
