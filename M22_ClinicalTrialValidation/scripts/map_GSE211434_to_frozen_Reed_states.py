#!/usr/bin/env python3
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[k]='4'
from pathlib import Path
import gzip,json,hashlib,time
import h5py
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
OUT=ROOT/'StageM22_ClinicalTrialStromalValidation_v0_1'
RAW=OUT/'raw'; DER=OUT/'derived'
HB=ROOT/'Reference_Arms/Normal_Reference/HBCA_Reed_Lineages/HBCA_stroma_be35c935.h5ad'
WA=ROOT/'Stage06B/objects/WAVE_A_COMMON_32229_COUNTS.h5ad'
MRK=ROOT/'Reference_Arms/Normal_Reference/Mapping_Benchmark/stroma_fine_v0_1/fine_state_markers.tsv'
ANN=RAW/'GSE211434_cell_types.txt.gz'
CNT=RAW/'GSE211434_counts.txt.gz'
CHUNK=5000

def sha(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def pearson_rows(X,C):
    X=np.asarray(X,dtype=np.float32); C=np.asarray(C,dtype=np.float32)
    Xc=X-X.mean(1,keepdims=True); Cc=C-C.mean(1,keepdims=True)
    xn=np.sqrt((Xc*Xc).sum(1)); cn=np.sqrt((Cc*Cc).sum(1))
    den=xn[:,None]*cn[None,:]
    return (Xc@Cc.T)/np.where(den>0,den,np.nan)

print('M22B_MAP_START',time.ctime(),flush=True)
markers=pd.read_csv(MRK,sep='\t')
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
print('FROZEN_COMMON_MARKERS',len(common),'of',len(requested),flush=True)

# Rebuild the frozen Reed donor-state centroid geometry exactly as in the verified historical mapper.
obs=hbmeta.obs[['donor_id','level2']].copy()
valid_states=set(markers.state.astype(str))
valid=obs.level2.astype(str).isin(valid_states).to_numpy()
group_key=np.where(valid,obs.donor_id.astype(str).to_numpy()+'|'+obs.level2.astype(str).to_numpy(),'')
ug=pd.Index(pd.unique(group_key[valid]))
gmap={g:i for i,g in enumerate(ug)}
gcode=np.full(len(obs),-1,dtype=np.int32)
for i,g in enumerate(group_key):
    if g:gcode[i]=gmap[g]
ng=len(ug); nm=len(common)
gsum=np.zeros((ng,nm),dtype=np.float64); gcount=np.zeros(ng,dtype=np.int64)
hbcols=np.array([hbpos[g] for g in common],dtype=np.int64)
with h5py.File(HB,'r') as f:
    indptr=f['X/indptr'][:]; data=f['X/data']; indices=f['X/indices']; nvar=int(f['X'].attrs['shape'][1])
    for r0 in range(0,len(obs),CHUNK):
        r1=min(r0+CHUNK,len(obs)); s=int(indptr[r0]); e=int(indptr[r1])
        dat=data[s:e]; idx=indices[s:e]; ptr=indptr[r0:r1+1]-s
        X=sparse.csr_matrix((dat,idx,ptr),shape=(r1-r0,nvar))
        Xm=X[:,hbcols]
        gc=gcode[r0:r1]; vr=np.flatnonzero(gc>=0)
        if len(vr):
            A=sparse.csr_matrix((np.ones(len(vr),dtype=np.float32),(gc[vr],vr)),shape=(ng,r1-r0))
            gsum+=(A@Xm).toarray(); gcount+=np.bincount(gc[vr],minlength=ng)
        if r0%100000==0:print('REED_ROWS',r0,'/',len(obs),flush=True)
gmean=gsum/np.maximum(gcount[:,None],1)
group_state=np.array([g.split('|',1)[1] for g in ug],dtype=object)
states_by={}; centroids_by={}
for lin,g in markers.groupby('lineage',sort=False):
    states=list(dict.fromkeys(g.sort_values(['state','rank']).state.astype(str).tolist()))
    midx=np.array(sorted(g.midx.unique()),dtype=int)
    C=np.vstack([gmean[group_state==st][:,midx].mean(0) for st in states]).astype(np.float32)
    states_by[lin]=states; centroids_by[lin]=(midx,C)
print('STATES',json.dumps(states_by),flush=True)

# Same centroid audit used in historical reconstruction.
errs=[]
for _,r in markers.iterrows():
    lin=str(r.lineage); st=str(r.state); gm=int(r.midx)
    midx,C=centroids_by[lin]; loc=np.where(midx==gm)[0]
    if len(loc):
        errs.append(float(C[states_by[lin].index(st),loc[0]])-float(r.class_mean))
audit={'n_marker_rows':len(errs),'rmse':float(np.sqrt(np.mean(np.square(errs)))),'max_abs_error':float(np.max(np.abs(errs)))}
print('CENTROID_AUDIT',json.dumps(audit),flush=True)

# Source metadata, with no response/timepoint information used as mapping features.
meta=pd.read_csv(ANN,sep='\t',dtype={'Day':int})
assert meta.CellID.is_unique and meta.Patient.nunique()==27
cross={'Fibroblasts':'Fibroblast','Pericytes':'Perivascular','Endothelial cells':'Endothelial'}
q=meta[meta.CellType.isin(cross)].copy().reset_index(drop=True)
q['mapper_lineage']=q.CellType.map(cross)
qids=q.CellID.astype(str).tolist()
qid_to_row={x:i for i,x in enumerate(qids)}
nq=len(q)
print('QUERY_STROMAL_CELLS',nq,q.CellType.value_counts().to_dict(),flush=True)

# Stream the deposited dense count table once; retain only frozen marker genes and only stromal query cells.
rows_acc=[]; cols_acc=[]; data_acc=[]
seen=set()
with gzip.open(CNT,'rt',errors='strict') as fh:
    header=fh.readline().rstrip('\n').split('\t')
    assert header[:3]==['Gene','Gene.ID','Gene.Symbol'],header[:3]
    cells=header[3:]
    assert len(cells)==len(meta)
    assert set(cells)==set(meta.CellID.astype(str))
    colpos={c:i for i,c in enumerate(cells)}
    selected=np.array([colpos[c] for c in qids],dtype=np.int64)
    for li,line in enumerate(fh,1):
        p1=line.find('\t'); p2=line.find('\t',p1+1); p3=line.find('\t',p2+1)
        ens=line[:p1]
        if ens not in g2m:continue
        vals=np.fromstring(line[p3+1:],sep='\t',dtype=np.int32)
        if len(vals)!=len(cells):raise RuntimeError((ens,len(vals),len(cells)))
        v=vals[selected]
        nz=np.flatnonzero(v)
        if len(nz):
            rows_acc.append(nz.astype(np.int32)); cols_acc.append(np.full(len(nz),g2m[ens],dtype=np.int32)); data_acc.append(v[nz].astype(np.int32))
        seen.add(ens)
        if len(seen)%200==0:print('QUERY_MARKERS_READ',len(seen),'/',len(common),flush=True)
assert seen==set(common),(len(set(common)-seen),list(sorted(set(common)-seen))[:10])
Crow=np.concatenate(rows_acc); Ccol=np.concatenate(cols_acc); Cdat=np.concatenate(data_acc)
QC=sparse.csr_matrix((Cdat,(Crow,Ccol)),shape=(nq,nm),dtype=np.int32)
print('QUERY_MARKER_MATRIX',QC.shape,'nnz',QC.nnz,flush=True)

# Source nCount_RNA is the deposited per-cell total-count field used as the denominator.
total=q.nCount_RNA.to_numpy(dtype=np.float64)
assert np.all(total>0)

pred=np.empty(nq,dtype=object); maxcorr=np.full(nq,np.nan); margin=np.full(nq,np.nan)
for lin in ['Fibroblast','Perivascular','Endothelial']:
    ii=np.flatnonzero(q.mapper_lineage.to_numpy()==lin)
    midx,C=centroids_by[lin]
    for s0 in range(0,len(ii),5000):
        jj=ii[s0:s0+5000]
        X=QC[jj][:,midx].toarray().astype(np.float32)
        t=total[jj]
        X=np.log1p(X*(10000.0/t[:,None])).astype(np.float32)
        R=pearson_rows(X,C)
        order=np.argsort(R,axis=1)
        b=order[:,-1]; sec=order[:,-2]
        pred[jj]=np.asarray(states_by[lin],dtype=object)[b]
        maxcorr[jj]=R[np.arange(len(jj)),b]
        margin[jj]=maxcorr[jj]-R[np.arange(len(jj)),sec]
    print('MAPPED',lin,len(ii),flush=True)

q['nearest_normal_state']=pred
q['max_corr']=maxcorr
q['margin']=margin

# Mapping outputs. Clinical columns are retained only now, after prediction.
cols=['CellID','Patient','Day','Treatment','TreatmentArm','Response','ResponseYN','CellType','CellTypeSubtype','Quality','nCount_RNA','nFeature_RNA','mapper_lineage','nearest_normal_state','max_corr','margin']
q[cols].to_csv(DER/'GSE211434_Reed_state_mapping.tsv.gz',sep='\t',index=False,compression='gzip')

summ=(q.groupby(['mapper_lineage','Day'])
      .agg(n_cells=('CellID','size'),n_patients=('Patient','nunique'),median_max_corr=('max_corr','median'),median_margin=('margin','median'))
      .reset_index())
summ.to_csv(DER/'GSE211434_Reed_mapping_confidence_summary.tsv',sep='\t',index=False)

state=(q.groupby(['Patient','Day','mapper_lineage','nearest_normal_state']).size().rename('n').reset_index())
tot=state.groupby(['Patient','Day','mapper_lineage']).n.transform('sum')
state['fraction']=state.n/tot
state.to_csv(DER/'GSE211434_Reed_patient_state_counts.tsv',sep='\t',index=False)

qc={
 'reference_h5ad':str(HB),'reference_sha256':sha(HB),
 'fine_state_markers':str(MRK),'fine_state_markers_sha256':sha(MRK),
 'annotation_sha256':sha(ANN),'counts_sha256':sha(CNT),
 'frozen_common_markers':len(common),'all_frozen_markers_present_in_GSE211434':True,
 'query_stromal_cells':nq,'lineage_cells':q.CellType.value_counts().to_dict(),
 'centroid_audit':audit,'response_or_treatment_used_for_mapping':False,
 'source_total_count_field':'nCount_RNA'
}
(DER/'M22B_MAPPING_QC.json').write_text(json.dumps(qc,indent=2))
print('M22B_MAPPING_PASS',json.dumps(qc),flush=True)
