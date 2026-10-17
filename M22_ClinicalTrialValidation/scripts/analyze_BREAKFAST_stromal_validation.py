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
from scipy.stats import wilcoxon

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
OUT=ROOT/'StageM22_ClinicalTrialStromalValidation_v0_1'
DER=OUT/'derived'
QH5=ROOT/'AtlasExpansion_Readiness_20261014_v0_9/GSE279219_PRIMARY_QC_COMMON35371_UNINTEGRATED.h5ad'
LIN=ROOT/'AtlasExpansion_Lineage_20261014_v0_10/GSE279219_broad_lineage_annotations.tsv.gz'
CW=ROOT/'AtlasExpansion_BREAKFAST_20261015_v0_31/BREAKFAST_GSE279219_patient_timepoint_crosswalk_v0_1.tsv'
FEATURE=ROOT/'AtlasExpansion_Acquisition_20261014_v0_6/GSE279219/GSM8564809_BRF14_features.tsv.gz'
HB=ROOT/'Reference_Arms/Normal_Reference/HBCA_Reed_Lineages/HBCA_stroma_be35c935.h5ad'
WA=ROOT/'Stage06B/objects/WAVE_A_COMMON_32229_COUNTS.h5ad'
MRK=ROOT/'Reference_Arms/Normal_Reference/Mapping_Benchmark/stroma_fine_v0_1/fine_state_markers.tsv'
CHUNK=5000

PROGRAMS={
 'Fibroblast':('fibroblast_matrix_contractile',['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','FN1','ACTA2','ITGA11','POSTN','JAG1']),
 'Perivascular':('mural_matrix_contractile',['COL1A1','COL1A2','COL3A1','ACTA2','JAG1']),
 'Endothelial':('endothelial_NOTCH_ECM',['JAG2','COL4A1','COL4A2','NOTCH4','ITGA1'])
}
GATES={
 'Fibroblast':(0.3168938385,0.00269698504),
 'Endothelial':(0.3484501570,0.0020778835),
 'Perivascular':(0.3765755535,0.003566771765)
}

def sha(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); prev=1.
    for k in range(n-1,-1,-1):
        i=o[k]; rank=k+1; v=min(prev,p[i]*n/rank); q[i]=v; prev=v
    return q
def pearson_rows(X,C):
    X=np.asarray(X,dtype=np.float32); C=np.asarray(C,dtype=np.float32)
    Xc=X-X.mean(1,keepdims=True); Cc=C-C.mean(1,keepdims=True)
    xn=np.sqrt((Xc*Xc).sum(1)); cn=np.sqrt((Cc*Cc).sum(1))
    den=xn[:,None]*cn[None,:]
    return (Xc@Cc.T)/np.where(den>0,den,np.nan)

print('M22C_START',time.ctime(),flush=True)

# Frozen pair crosswalk.
cw=pd.read_csv(CW,sep='\t')
cw=cw[cw.pair_status.eq('COMPLETE_T0_T1')].copy()
assert cw.patient_label.nunique()==10 and len(cw)==20
sample_to_tp=dict(zip(cw.sample_label,cw.timepoint_code))
sample_to_pt=dict(zip(cw.sample_label,cw.patient_label))

# Frozen broad lineage calls.
lin=pd.read_csv(LIN,sep='\t')
lin=lin[lin['sample'].isin(sample_to_tp)].copy()
lin['Patient']=lin['sample'].map(sample_to_pt); lin['Timepoint']=lin['sample'].map(sample_to_tp)
nonamb=lin[~lin.provisional_joint_lineage.eq('Ambiguous')].copy()
targets={'Fibroblast':'Fibroblast','Endothelial':'Endothelial','Perivascular':'Perivascular'}

# Family A: broad abundance.
abvals=[]; abtests=[]; abloo=[]
for L in targets:
    for (pt,tp),g in nonamb.groupby(['Patient','Timepoint']):
        total=len(g); n=int((g.provisional_joint_lineage==L).sum())
        abvals.append(dict(Patient=pt,Timepoint=tp,lineage=L,n_lineage=n,total_nonambiguous=total,
                           fraction=n/total,balance=np.log((n+0.5)/(total-n+0.5))))
    z=pd.DataFrame([x for x in abvals if x['lineage']==L])
    w=z.pivot(index='Patient',columns='Timepoint',values='balance').dropna()
    f=z.pivot(index='Patient',columns='Timepoint',values='fraction').reindex(w.index)
    delta=(w['T1']-w['T0']).to_numpy()
    p=float(wilcoxon(w['T1'],w['T0'],alternative='two-sided').pvalue)
    med=float(np.median(delta))
    abtests.append(dict(lineage=L,n=len(w),median_balance_change=med,
                        median_fraction_T0=float(np.median(f['T0'])),median_fraction_T1=float(np.median(f['T1'])),p=p))
    for i,pt in enumerate(w.index):
        lm=float(np.median(np.delete(delta,i)))
        abloo.append(dict(lineage=L,omitted_patient=pt,full_median_change=med,loo_median_change=lm,
                          same_direction=bool(np.sign(lm)==np.sign(med) or lm==0 or med==0)))
abtests=pd.DataFrame(abtests);abtests['FDR']=bh(abtests.p.to_numpy())
pd.DataFrame(abvals).to_csv(DER/'BREAKFAST_broad_abundance_patient_values.tsv',sep='\t',index=False)
abtests.to_csv(DER/'BREAKFAST_broad_abundance_stats.tsv',sep='\t',index=False)
pd.DataFrame(abloo).to_csv(DER/'BREAKFAST_broad_abundance_LOO.tsv',sep='\t',index=False)

# Query H5AD + exact lineages.
qh=ad.read_h5ad(QH5,backed='r')
assert qh.obs_names.astype(str).equals(pd.Index(qh.obs.cell_id.astype(str)))
qvar=list(qh.var_names.astype(str)); qpos={g:i for i,g in enumerate(qvar)}
qobs=qh.obs[['cell_id','sample','total_counts']].copy().reset_index(drop=True)
qobs=qobs.merge(lin[['cell_id','provisional_joint_lineage','Patient','Timepoint']],
                 on='cell_id',how='left',validate='one_to_one')
q=qobs[qobs.provisional_joint_lineage.isin(targets)].copy().reset_index(drop=True)
q['mapper_lineage']=q.provisional_joint_lineage
qrows=pd.Index(qh.obs_names.astype(str)).get_indexer(q.cell_id.astype(str))
assert (qrows>=0).all()
q['h5row']=qrows
print('BREAKFAST_TARGET_CELLS',len(q),q.mapper_lineage.value_counts().to_dict(),flush=True)

# Frozen Reed marker universe and centroids.
markers=pd.read_csv(MRK,sep='\t')
hbmeta=ad.read_h5ad(HB,backed='r'); wameta=ad.read_h5ad(WA,backed='r')
hbpos={g:i for i,g in enumerate(hbmeta.var_names.astype(str))}
wapos={g:i for i,g in enumerate(wameta.var_names.astype(str))}
requested=set(markers.ensembl_id.astype(str))
common=[g for g in requested if g in hbpos and g in wapos]
common=sorted(common,key=lambda g:hbpos[g])
assert set(common)<=set(qvar),(len(set(common)-set(qvar)),list(set(common)-set(qvar))[:10])
g2m={g:i for i,g in enumerate(common)}
markers=markers[markers.ensembl_id.astype(str).isin(g2m)].copy()
markers['midx']=markers.ensembl_id.astype(str).map(g2m).astype(int)

obs=hbmeta.obs[['donor_id','level2']].copy()
valid_states=set(markers.state.astype(str)); valid=obs.level2.astype(str).isin(valid_states).to_numpy()
group_key=np.where(valid,obs.donor_id.astype(str).to_numpy()+'|'+obs.level2.astype(str).to_numpy(),'')
ug=pd.Index(pd.unique(group_key[valid])); gmap={g:i for i,g in enumerate(ug)}
gcode=np.full(len(obs),-1,dtype=np.int32)
for i,g in enumerate(group_key):
    if g:gcode[i]=gmap[g]
ng=len(ug); nm=len(common); gsum=np.zeros((ng,nm),dtype=np.float64); gcount=np.zeros(ng,dtype=np.int64)
hbcols=np.array([hbpos[g] for g in common],dtype=np.int64)
with h5py.File(HB,'r') as f:
    indptr=f['X/indptr'][:]; data=f['X/data']; indices=f['X/indices']; nvar=int(f['X'].attrs['shape'][1])
    for r0 in range(0,len(obs),CHUNK):
        r1=min(r0+CHUNK,len(obs)); s=int(indptr[r0]); e=int(indptr[r1])
        X=sparse.csr_matrix((data[s:e],indices[s:e],indptr[r0:r1+1]-s),shape=(r1-r0,nvar))
        gc=gcode[r0:r1]; vr=np.flatnonzero(gc>=0)
        if len(vr):
            A=sparse.csr_matrix((np.ones(len(vr),dtype=np.float32),(gc[vr],vr)),shape=(ng,r1-r0))
            gsum+=(A@X[:,hbcols]).toarray();gcount+=np.bincount(gc[vr],minlength=ng)
gmean=gsum/np.maximum(gcount[:,None],1)
group_state=np.array([g.split('|',1)[1] for g in ug],dtype=object)
states_by={};centroids_by={}
for L,g in markers.groupby('lineage',sort=False):
    states=list(dict.fromkeys(g.sort_values(['state','rank']).state.astype(str).tolist()))
    midx=np.array(sorted(g.midx.unique()),dtype=int)
    C=np.vstack([gmean[group_state==st][:,midx].mean(0) for st in states]).astype(np.float32)
    states_by[L]=states;centroids_by[L]=(midx,C)

# Gene-symbol -> Ensembl map from original GSE279219 feature file.
sym2ens={}
with gzip.open(FEATURE,'rt') as f:
    for line in f:
        parts=line.rstrip('\n').split('\t')
        if len(parts)>=2:sym2ens.setdefault(parts[1],[]).append(parts[0])
all_program_symbols=sorted(set(x for _,genes in PROGRAMS.values() for x in genes))
program_ens={s:[e for e in sym2ens.get(s,[]) if e in qpos] for s in all_program_symbols}
assert all(program_ens[s] for s in all_program_symbols),program_ens
prog_genes=sorted(set(e for es in program_ens.values() for e in es),key=lambda e:qpos[e])
prog_pos={e:i for i,e in enumerate(prog_genes)}

# Extract query marker + program matrices directly from sparse H5AD counts.
full_to_q=np.full(qh.n_obs,-1,dtype=np.int32);full_to_q[qrows]=np.arange(len(q),dtype=np.int32)
mark_cols=np.array([qpos[g] for g in common],dtype=np.int64)
prog_cols=np.array([qpos[g] for g in prog_genes],dtype=np.int64)
mr=[];mc=[];md=[];pr=[];pc=[];pdv=[]
with h5py.File(QH5,'r') as f:
    indptr=f['X/indptr'][:]; data=f['X/data']; indices=f['X/indices']; nvar=int(f['X'].attrs['shape'][1])
    for r0 in range(0,qh.n_obs,CHUNK):
        r1=min(r0+CHUNK,qh.n_obs); qi=full_to_q[r0:r1]; sel=np.flatnonzero(qi>=0)
        if not len(sel):continue
        s=int(indptr[r0]);e=int(indptr[r1])
        X=sparse.csr_matrix((data[s:e],indices[s:e],indptr[r0:r1+1]-s),shape=(r1-r0,nvar))
        qr=qi[sel]
        A=X[sel][:,mark_cols].tocoo()
        if A.nnz:
            mr.append(qr[A.row]);mc.append(A.col.astype(np.int32));md.append(A.data.astype(np.float32))
        P=X[sel][:,prog_cols].tocoo()
        if P.nnz:
            pr.append(qr[P.row]);pc.append(P.col.astype(np.int32));pdv.append(P.data.astype(np.float32))
Qm=sparse.csr_matrix((np.concatenate(md),(np.concatenate(mr),np.concatenate(mc))),shape=(len(q),len(common)),dtype=np.float32)
Qp=sparse.csr_matrix((np.concatenate(pdv),(np.concatenate(pr),np.concatenate(pc))),shape=(len(q),len(prog_genes)),dtype=np.float32)
print('QUERY_MATRICES',Qm.shape,Qm.nnz,Qp.shape,Qp.nnz,flush=True)

# Family C: map to frozen Reed states.
pred=np.empty(len(q),dtype=object);maxcorr=np.full(len(q),np.nan);margin=np.full(len(q),np.nan)
tot=q.total_counts.to_numpy(float)
for L in ['Fibroblast','Perivascular','Endothelial']:
    ii=np.flatnonzero(q.mapper_lineage.to_numpy()==L);midx,C=centroids_by[L]
    for s0 in range(0,len(ii),5000):
        jj=ii[s0:s0+5000]
        X=Qm[jj][:,midx].toarray().astype(np.float32)
        X=np.log1p(X*(10000.0/tot[jj,None])).astype(np.float32)
        R=pearson_rows(X,C);o=np.argsort(R,axis=1);b=o[:,-1];sec=o[:,-2]
        pred[jj]=np.asarray(states_by[L],object)[b]
        maxcorr[jj]=R[np.arange(len(jj)),b];margin[jj]=maxcorr[jj]-R[np.arange(len(jj)),sec]
q['nearest_normal_state']=pred;q['max_corr']=maxcorr;q['margin']=margin

app=[]
for L,g in q.groupby('mapper_lineage'):
    cm=float(g.max_corr.median());mm=float(g.margin.median());cg,mg=GATES[L]
    app.append(dict(lineage=L,n_cells=len(g),median_max_corr=cm,max_corr_gate=cg,
                    median_margin=mm,margin_gate=mg,applicability_pass=bool(cm>=cg and mm>=mg)))
app=pd.DataFrame(app);app.to_csv(DER/'BREAKFAST_Reed_applicability_gate.tsv',sep='\t',index=False)

q[['cell_id','sample','Patient','Timepoint','mapper_lineage','nearest_normal_state','max_corr','margin','total_counts']].to_csv(
    DER/'BREAKFAST_Reed_state_mapping.tsv.gz',sep='\t',index=False,compression='gzip')

statevals=[];statetests=[];stateloo=[]
for L,g in q.groupby('mapper_lineage'):
    st=states_by[L]
    for (pt,tp),z in g.groupby(['Patient','Timepoint']):
        c=z.nearest_normal_state.value_counts();d={'Patient':pt,'Timepoint':tp,'lineage':L,'n_cells':len(z)}
        for s in st:d[s]=int(c.get(s,0));d[s+'_fraction']=d[s]/len(z)
        if L=='Fibroblast':
            d['endpoint']='FB13_vs_FB24';d['balance']=np.log((d['FB1']+d['FB3']+0.5)/(d['FB2']+d['FB4']+0.5))
        elif L=='Perivascular':
            d['endpoint']='PV3_vs_rest';d['balance']=np.log((d['PV3']+0.5)/(d['PV1']+d['PV2']+d['PV4']+d['PV5']+0.5))
        else:
            d['endpoint']='LE2_vs_LE1';d['balance']=np.log((d['LE2']+0.5)/(d['LE1']+0.5))
        statevals.append(d)
sv=pd.DataFrame(statevals)
for L,z in sv.groupby('lineage'):
    w=z.pivot(index='Patient',columns='Timepoint',values='balance').dropna()
    delta=(w['T1']-w['T0']).to_numpy();p=float(wilcoxon(w['T1'],w['T0'],alternative='two-sided').pvalue);med=float(np.median(delta))
    statetests.append(dict(lineage=L,endpoint=z.endpoint.iloc[0],n=len(w),median_balance_change=med,p=p))
    for i,pt in enumerate(w.index):
        lm=float(np.median(np.delete(delta,i)));stateloo.append(dict(lineage=L,omitted_patient=pt,full_median_change=med,loo_median_change=lm,
                same_direction=bool(np.sign(lm)==np.sign(med) or lm==0 or med==0)))
statetests=pd.DataFrame(statetests);statetests['FDR']=bh(statetests.p.to_numpy())
statetests=statetests.merge(app[['lineage','applicability_pass']],on='lineage',how='left')
sv.to_csv(DER/'BREAKFAST_Reed_state_patient_values.tsv',sep='\t',index=False)
statetests.to_csv(DER/'BREAKFAST_Reed_state_stats.tsv',sep='\t',index=False)
pd.DataFrame(stateloo).to_csv(DER/'BREAKFAST_Reed_state_LOO.tsv',sep='\t',index=False)

# Family B: fixed program pseudobulks.
group_key=q.Patient.astype(str)+'|'+q.Timepoint.astype(str)+'|'+q.mapper_lineage.astype(str)
ug=pd.Index(pd.unique(group_key));gmap={g:i for i,g in enumerate(ug)}
gc=np.array([gmap[g] for g in group_key],dtype=np.int32)
A=sparse.csr_matrix((np.ones(len(q),dtype=np.float32),(gc,np.arange(len(q)))),shape=(len(ug),len(q)))
gsum=(A@Qp).toarray()
glib=np.bincount(gc,weights=tot,minlength=len(ug));gn=np.bincount(gc,minlength=len(ug))
gi=pd.DataFrame([x.split('|') for x in ug],columns=['Patient','Timepoint','lineage'])
gi['n_cells']=gn;gi['library_size']=glib
for sym,enss in program_ens.items():
    cols=[prog_pos[e] for e in enss]
    count=gsum[:,cols].sum(1)
    gi[sym+'_count']=count
    gi[sym+'_log2CPM']=np.log2(count/glib*1e6+1)

progvals=[];progtests=[];progloo=[]
for L,(pname,genes) in PROGRAMS.items():
    z=gi[gi.lineage.eq(L)].copy()
    for gene in genes:
        c=gene+'_log2CPM';sd=z[c].std(ddof=1);z[gene+'_z']=(z[c]-z[c].mean())/sd if sd>0 else 0.
    z['program']=pname;z['program_score']=z[[x+'_z' for x in genes]].mean(1)
    progvals.append(z)
    w=z.pivot(index='Patient',columns='Timepoint',values='program_score').dropna()
    delta=(w['T1']-w['T0']).to_numpy();p=float(wilcoxon(w['T1'],w['T0'],alternative='two-sided').pvalue);med=float(np.median(delta))
    progtests.append(dict(lineage=L,program=pname,n=len(w),median_score_change=med,p=p))
    for i,pt in enumerate(w.index):
        lm=float(np.median(np.delete(delta,i)));progloo.append(dict(lineage=L,omitted_patient=pt,full_median_change=med,loo_median_change=lm,
                same_direction=bool(np.sign(lm)==np.sign(med) or lm==0 or med==0)))
progvals=pd.concat(progvals,ignore_index=True);progtests=pd.DataFrame(progtests);progtests['FDR']=bh(progtests.p.to_numpy())
progvals.to_csv(DER/'BREAKFAST_program_patient_values.tsv.gz',sep='\t',index=False,compression='gzip')
progtests.to_csv(DER/'BREAKFAST_program_stats.tsv',sep='\t',index=False)
pd.DataFrame(progloo).to_csv(DER/'BREAKFAST_program_LOO.tsv',sep='\t',index=False)

gene_map=[]
for s in all_program_symbols:
    gene_map.append(dict(symbol=s,ensembl_ids=';'.join(program_ens[s]),n_ids=len(program_ens[s])))
pd.DataFrame(gene_map).to_csv(DER/'BREAKFAST_program_gene_map.tsv',sep='\t',index=False)

qc={
 'paired_patients':10,'paired_samples':20,'target_cells':len(q),
 'target_lineage_cells':q.mapper_lineage.value_counts().to_dict(),
 'frozen_common_markers':len(common),'all_frozen_markers_present':True,
 'applicability':app.to_dict(orient='records'),
 'response_analysis_performed':False,
 'inputs':{'query_h5ad_sha256':sha(QH5),'lineage_sha256':sha(LIN),'crosswalk_sha256':sha(CW),
           'feature_sha256':sha(FEATURE),'Reed_sha256':sha(HB),'marker_sha256':sha(MRK)}
}
(DER/'M22C_BREAKFAST_QC.json').write_text(json.dumps(qc,indent=2))
print('\nAPPLICABILITY');print(app.to_string(index=False))
print('\nABUNDANCE');print(abtests.to_string(index=False))
print('\nPROGRAMS');print(progtests.to_string(index=False))
print('\nREED_STATES');print(statetests.to_string(index=False))
print('\nM22C_PASS',json.dumps(qc),flush=True)
