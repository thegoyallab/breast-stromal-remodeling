#!/usr/bin/env python3
import os
from pathlib import Path
import json, time
import h5py
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT=Path(os.environ.get('BC_PROJECT_ROOT', Path.cwd())).resolve()
BR=ROOT/'StageM15F_FinalScientificFreeze_v0_1'
OUT=BR/'results/s4_mapper'; LOG=BR/'audit'
OUT.mkdir(parents=True,exist_ok=True); LOG.mkdir(parents=True,exist_ok=True)
HB=ROOT/'Reference_Arms/Normal_Reference/HBCA_Reed_Lineages/HBCA_stroma_be35c935.h5ad'
WA=ROOT/'Stage06B/objects/WAVE_A_COMMON_32229_COUNTS.h5ad'
MAP=ROOT/'Reference_Arms/Normal_Reference/Mapping_Benchmark/stroma_fine_v0_1/WaveA_stroma_fine_mapping.tsv.gz'
MRK=ROOT/'Reference_Arms/Normal_Reference/Mapping_Benchmark/stroma_fine_v0_1/fine_state_markers.tsv'
SEED=20261011
RNG=np.random.default_rng(SEED)
CHUNK=5000
RETENTIONS=[1.0,0.75,0.50,0.25,0.10]
NREP=5

def pearson_rows(X, C):
    X=np.asarray(X,dtype=np.float32); C=np.asarray(C,dtype=np.float32)
    Xc=X-X.mean(1,keepdims=True)
    Cc=C-C.mean(1,keepdims=True)
    xn=np.sqrt((Xc*Xc).sum(1)); cn=np.sqrt((Cc*Cc).sum(1))
    den=xn[:,None]*cn[None,:]
    out=(Xc@Cc.T)/np.where(den>0,den,np.nan)
    return out

def scores_from_counts(Csp,total,marker_idx,centroids,states):
    X=Csp[:,marker_idx].toarray().astype(np.float32)
    t=np.asarray(total,dtype=np.float64)
    ok=t>0
    X[ok]=np.log1p(X[ok]*(10000.0/t[ok,None])).astype(np.float32)
    X[~ok]=0
    R=pearson_rows(X,centroids)
    order=np.argsort(R,axis=1)
    best=order[:,-1]; second=order[:,-2]
    maxc=R[np.arange(len(R)),best]
    margin=maxc-R[np.arange(len(R)),second]
    pred=np.asarray(states,dtype=object)[best]
    return pred,maxc,margin

def thin_counts(Csp,total,p,rng):
    C=Csp.tocsr(copy=False)
    if np.ndim(p)==0:
        td=rng.binomial(C.data.astype(np.int64),float(p)).astype(np.int32)
        nonmarker=np.maximum(np.asarray(total,dtype=np.int64)-np.asarray(C.sum(1)).ravel().astype(np.int64),0)
        other=rng.binomial(nonmarker,float(p)).astype(np.int64)
    else:
        p=np.asarray(p,dtype=float)
        row=np.repeat(np.arange(C.shape[0]),np.diff(C.indptr))
        td=rng.binomial(C.data.astype(np.int64),p[row]).astype(np.int32)
        nonmarker=np.maximum(np.asarray(total,dtype=np.int64)-np.asarray(C.sum(1)).ravel().astype(np.int64),0)
        other=rng.binomial(nonmarker,p).astype(np.int64)
    T=sparse.csr_matrix((td,C.indices.copy(),C.indptr.copy()),shape=C.shape)
    T.eliminate_zeros()
    tt=np.asarray(T.sum(1)).ravel().astype(np.int64)+other
    return T,tt

print('M11A3 start',time.ctime(),flush=True)
markers=pd.read_csv(MRK,sep='\t')
mapping=pd.read_csv(MAP,sep='\t')
dcis={'GSE195861:NCCBC2','GSE195861:NCCBC3','GSE195861:NCCBC5','GSE195861:NCCBC6','GSE195861:NCCBC7','GSE195861:NCCBC11','GSE195861:NCCBC14'}
mapping.loc[mapping.patient_uid.astype(str).isin(dcis),'stage06a_role']='reference_dcis'
mapping=mapping[mapping.stage06a_role.isin(['reference_normal','core_primary'])].copy().reset_index(drop=True)

# Common marker universe in deterministic HBCA-column order.
hbmeta=ad.read_h5ad(HB,backed='r')
wameta=ad.read_h5ad(WA,backed='r')
hbpos={g:i for i,g in enumerate(hbmeta.var_names.astype(str))}
wapos={g:i for i,g in enumerate(wameta.var_names.astype(str))}
requested=set(markers.ensembl_id.astype(str))
common=[g for g in requested if g in hbpos and g in wapos]
common=sorted(common,key=lambda g:hbpos[g])
g2m={g:i for i,g in enumerate(common)}
markers=markers[markers.ensembl_id.astype(str).isin(g2m)].copy()
markers['midx']=markers.ensembl_id.astype(str).map(g2m).astype(int)
print('common markers',len(common),'of',len(requested),flush=True)

# Build donor-state pseudobulk HBCA centroids, matching the original LODO geometry.
obs=hbmeta.obs[['donor_id','level2']].copy()
valid_states=set(markers.state.astype(str))
valid=obs.level2.astype(str).isin(valid_states).to_numpy()
group_key=np.where(valid,obs.donor_id.astype(str).to_numpy()+'|'+obs.level2.astype(str).to_numpy(),'')
ug=pd.Index(pd.unique(group_key[valid]))
gmap={g:i for i,g in enumerate(ug)}
gcode=np.full(len(obs),-1,dtype=np.int32)
for i,g in enumerate(group_key):
    if g: gcode[i]=gmap[g]
ng=len(ug); nm=len(common)
gsum=np.zeros((ng,nm),dtype=np.float64); gcount=np.zeros(ng,dtype=np.int64)
hbcols=np.array([hbpos[g] for g in common],dtype=np.int64)
with h5py.File(HB,'r') as f:
    indptr=f['X/indptr'][:]
    data=f['X/data']; indices=f['X/indices']; nvar=int(f['X'].attrs['shape'][1])
    for r0 in range(0,len(obs),CHUNK):
        r1=min(r0+CHUNK,len(obs)); s=int(indptr[r0]); e=int(indptr[r1])
        dat=data[s:e]; idx=indices[s:e]; ptr=indptr[r0:r1+1]-s
        X=sparse.csr_matrix((dat,idx,ptr),shape=(r1-r0,nvar))
        Xm=X[:,hbcols]
        gc=gcode[r0:r1]; vr=np.flatnonzero(gc>=0)
        if len(vr):
            A=sparse.csr_matrix((np.ones(len(vr),dtype=np.float32),(gc[vr],vr)),shape=(ng,r1-r0))
            gsum += (A@Xm).toarray()
            gcount += np.bincount(gc[vr],minlength=ng)
        if r0 % 50000==0: print('HBCA rows',r0,'/',len(obs),flush=True)
gmean=gsum/np.maximum(gcount[:,None],1)
group_state=np.array([g.split('|',1)[1] for g in ug],dtype=object)
states_by={}
centroids_by={}
for lin,g in markers.groupby('lineage'):
    states=list(dict.fromkeys(g.sort_values(['state','rank']).state.astype(str).tolist()))
    midx=np.array(sorted(g.midx.unique()),dtype=int)
    C=[]
    for st in states:
        C.append(gmean[group_state==st][:,midx].mean(0))
    states_by[lin]=states
    centroids_by[lin]=(midx,np.vstack(C).astype(np.float32))

# Verify saved marker class means are reproduced by donor-state centroid construction.
errs=[]
for _,r in markers.iterrows():
    lin=str(r.lineage); st=str(r.state); gm=int(r.midx)
    midx,C=centroids_by[lin]
    loc=np.where(midx==gm)[0]
    if not len(loc): continue
    val=C[states_by[lin].index(st),loc[0]]
    errs.append(float(val)-float(r.class_mean))
centroid_audit={'n_marker_rows':len(errs),'rmse':float(np.sqrt(np.mean(np.square(errs)))),'max_abs_error':float(np.max(np.abs(errs)))}
print('centroid audit',centroid_audit,flush=True)

# Extract Wave-A cells in mapping table, their total UMIs and counts for common markers.
wa_index=pd.Index(wameta.obs_names.astype(str))
rowidx=wa_index.get_indexer(mapping.cell_id.astype(str))
if (rowidx<0).any(): raise RuntimeError(f'missing mapped cells: {(rowidx<0).sum()}')
row_to_map=np.full(wameta.n_obs,-1,dtype=np.int32); row_to_map[rowidx]=np.arange(len(mapping),dtype=np.int32)
wacols=np.array([wapos[g] for g in common],dtype=np.int64)
rows_acc=[]; cols_acc=[]; dat_acc=[]
tot=np.zeros(len(mapping),dtype=np.int64)
with h5py.File(WA,'r') as f:
    indptr=f['X/indptr'][:]
    data=f['X/data']; indices=f['X/indices']; nvar=int(f['X'].attrs['shape'][1])
    for r0 in range(0,wameta.n_obs,CHUNK):
        r1=min(r0+CHUNK,wameta.n_obs); mi=row_to_map[r0:r1]; sel=np.flatnonzero(mi>=0)
        if not len(sel): continue
        s=int(indptr[r0]); e=int(indptr[r1])
        dat=data[s:e]; idx=indices[s:e]; ptr=indptr[r0:r1+1]-s
        X=sparse.csr_matrix((dat,idx,ptr),shape=(r1-r0,nvar))
        ms=mi[sel]
        tot[ms]=np.asarray(X[sel].sum(1)).ravel().astype(np.int64)
        Xm=X[sel][:,wacols].tocoo()
        rows_acc.append(ms[Xm.row]); cols_acc.append(Xm.col.astype(np.int32)); dat_acc.append(Xm.data.astype(np.int32))
        if r0 % 50000==0: print('WaveA rows',r0,'/',wameta.n_obs,flush=True)
Crow=np.concatenate(rows_acc); Ccol=np.concatenate(cols_acc); Cdat=np.concatenate(dat_acc)
WC=sparse.csr_matrix((Cdat,(Crow,Ccol)),shape=(len(mapping),nm),dtype=np.int32)
print('WaveA mapped matrix',WC.shape,'nnz',WC.nnz,flush=True)

# Reproduce original mapping numerically.
rep_rows=[]; base_pred=np.empty(len(mapping),dtype=object); base_corr=np.full(len(mapping),np.nan); base_margin=np.full(len(mapping),np.nan)
for lin in states_by:
    ii=np.flatnonzero(mapping.lineage.astype(str).to_numpy()==lin)
    midx,C=centroids_by[lin]; st=states_by[lin]
    pr,mc,ma=scores_from_counts(WC[ii],tot[ii],midx,C,st)
    base_pred[ii]=pr; base_corr[ii]=mc; base_margin[ii]=ma
    orig=mapping.iloc[ii]
    concord=(pr==orig.nearest_normal_state.astype(str).to_numpy()).mean()
    rho_mc=spearmanr(mc,orig.max_corr.to_numpy(float),nan_policy='omit').statistic
    rho_ma=spearmanr(ma,orig.margin.to_numpy(float),nan_policy='omit').statistic
    rep_rows.append([lin,len(ii),concord,rho_mc,float(np.nanmedian(np.abs(mc-orig.max_corr.to_numpy(float)))),rho_ma,float(np.nanmedian(np.abs(ma-orig.margin.to_numpy(float))))])
repro=pd.DataFrame(rep_rows,columns=['lineage','n_cells','state_concordance','spearman_max_corr','median_abs_error_max_corr','spearman_margin','median_abs_error_margin'])
repro.to_csv(OUT/'M11A3_mapping_reproduction_audit_v0_1.tsv',sep='\t',index=False)
print(repro.to_string(index=False),flush=True)

# Baseline patient metrics and state composition under the reproduced mapper.
mapping['repro_state']=base_pred; mapping['repro_max_corr']=base_corr; mapping['repro_margin']=base_margin
base_patient=(mapping.groupby(['lineage','stage06a_role','patient_uid'])
              .agg(n_cells=('cell_id','size'),median_max_corr=('repro_max_corr','median'),median_margin=('repro_margin','median'))
              .reset_index())
base_patient.to_csv(OUT/'M11A3_reproduced_patient_confidence_v0_1.tsv',sep='\t',index=False)


# Save invasive-only reproduced state fractions for M15F correction.
frac=(mapping.groupby(['lineage','stage06a_role','patient_uid','repro_state']).size().rename('n').reset_index())
tot=frac.groupby(['lineage','stage06a_role','patient_uid']).n.transform('sum')
frac['fraction']=frac.n/tot
frac.to_csv(OUT/'M15F_reproduced_patient_state_fractions_DCIscorrected_v0_1.tsv',sep='\t',index=False)
print('M15F mapper-only PASS',len(mapping),flush=True)
