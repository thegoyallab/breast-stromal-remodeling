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
OUT=BR/'results/s4_full'; LOG=BR/'audit'
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

# Library-depth audit.
depth=(mapping.assign(total_umi=tot).groupby(['lineage','stage06a_role'])
       .agg(n_cells=('cell_id','size'),n_patients=('patient_uid','nunique'),median_total_umi=('total_umi','median'),
            q25_total_umi=('total_umi',lambda x:x.quantile(.25)),q75_total_umi=('total_umi',lambda x:x.quantile(.75)))
       .reset_index())
depth.to_csv(OUT/'M11A3_library_depth_by_lineage_role_v0_1.tsv',sep='\t',index=False)
print('depth\n',depth.to_string(index=False),flush=True)

cal_pat=[]; cal_state=[]; stability=[]
for lin in states_by:
    il=np.flatnonzero(mapping.lineage.astype(str).to_numpy()==lin)
    locmap=mapping.iloc[il].reset_index(drop=True)
    Cline=WC[il].tocsr(); tline=tot[il]
    midx,cent=centroids_by[lin]; sts=states_by[lin]
    norm=np.flatnonzero(locmap.stage06a_role.astype(str).to_numpy()=='reference_normal')
    tum=np.flatnonzero(locmap.stage06a_role.astype(str).to_numpy()=='core_primary')
    Cn=Cline[norm]; tn=tline[norm]; mn=locmap.iloc[norm].reset_index(drop=True)
    baseline_state=base_pred[il][norm]
    tumor_depth=tline[tum]
    conditions=[]
    for p in RETENTIONS:
        conditions.append((f'retain_{int(round(p*100)):03d}',p))
    conditions.append(('tumor_depth_matched',None))
    for cond,p in conditions:
        for rep in range(NREP if p!=1.0 else 1):
            rng=np.random.default_rng(SEED+100000*list(states_by).index(lin)+1000*conditions.index((cond,p))+rep)
            if p is None:
                targets=rng.choice(tumor_depth,size=len(tn),replace=True)
                pv=np.minimum(1.0,targets/np.maximum(tn,1))
                Ct,tt=thin_counts(Cn,tn,pv,rng)
                pmedian=float(np.median(pv))
            elif p==1.0:
                Ct,tt=Cn,tn.copy(); pmedian=1.0
            else:
                Ct,tt=thin_counts(Cn,tn,p,rng); pmedian=float(p)
            pr,mc,ma=scores_from_counts(Ct,tt,midx,cent,sts)
            stability.append([lin,cond,rep,pmedian,len(pr),float(np.mean(pr==baseline_state)),float(np.median(mc)),float(np.median(ma)),float(np.median(tt))])
            tmp=pd.DataFrame({'patient_uid':mn.patient_uid.astype(str),'state':pr,'max_corr':mc,'margin':ma})
            for pid,gp in tmp.groupby('patient_uid'):
                cal_pat.append([lin,cond,rep,pmedian,pid,len(gp),float(np.median(gp.max_corr)),float(np.median(gp.margin))])
                for st in sts:
                    cal_state.append([lin,cond,rep,pmedian,pid,st,float(np.mean(gp.state==st))])

CP=pd.DataFrame(cal_pat,columns=['lineage','condition','replicate','median_retention','patient_uid','n_cells','median_max_corr','median_margin'])
CS=pd.DataFrame(cal_state,columns=['lineage','condition','replicate','median_retention','patient_uid','state','state_fraction'])
ST=pd.DataFrame(stability,columns=['lineage','condition','replicate','median_retention','n_cells','label_stability_vs_baseline','cell_median_max_corr','cell_median_margin','cell_median_total_umi'])
CP.to_csv(OUT/'M11A3_degradation_patient_confidence_v0_1.tsv',sep='\t',index=False)
CS.to_csv(OUT/'M11A3_degradation_patient_state_fractions_v0_1.tsv',sep='\t',index=False)
ST.to_csv(OUT/'M11A3_degradation_cell_summary_v0_1.tsv',sep='\t',index=False)

# Collapse replicate simulations to patient means, compare with reproduced primary-vs-normal shift.
summary=[]; state_summary=[]; paired_tests=[]
for lin in states_by:
    bp=base_patient[base_patient.lineage==lin]
    n0=bp[bp.stage06a_role=='reference_normal']; tp=bp[bp.stage06a_role=='core_primary']
    nmc=float(n0.median_max_corr.median()); tmc=float(tp.median_max_corr.median())
    nma=float(n0.median_margin.median()); tma=float(tp.median_margin.median())
    for cond,g in CP[CP.lineage==lin].groupby('condition'):
        gp=g.groupby('patient_uid')[['median_max_corr','median_margin']].mean()
        cmc=float(gp.median_max_corr.median()); cma=float(gp.median_margin.median())
        frac_mc=(nmc-cmc)/(nmc-tmc) if nmc!=tmc else np.nan
        frac_ma=(nma-cma)/(nma-tma) if nma!=tma else np.nan
        summary.append([lin,cond,len(gp),nmc,tmc,cmc,frac_mc,nma,tma,cma,frac_ma,
                        float(ST[(ST.lineage==lin)&(ST.condition==cond)].label_stability_vs_baseline.mean()),
                        float(ST[(ST.lineage==lin)&(ST.condition==cond)].median_retention.mean())])
        # paired normal baseline vs degraded patient metrics
        base=n0.set_index('patient_uid')
        commonp=base.index.intersection(gp.index)
        for metric in ['median_max_corr','median_margin']:
            x=base.loc[commonp,metric].to_numpy(float); y=gp.loc[commonp,metric].to_numpy(float)
            try: stat,pv=wilcoxon(x,y,zero_method='wilcox',alternative='two-sided')
            except Exception: stat,pv=np.nan,np.nan
            paired_tests.append([lin,cond,metric,len(commonp),float(np.median(x)),float(np.median(y)),float(np.median(y-x)),stat,pv])
    # state fractions: reproduced baseline normal and primary
    sub=mapping[mapping.lineage==lin]
    for st in states_by[lin]:
        nfrac=sub[sub.stage06a_role=='reference_normal'].groupby('patient_uid').repro_state.apply(lambda x:float(np.mean(np.asarray(x)==st)))
        tfrac=sub[sub.stage06a_role=='core_primary'].groupby('patient_uid').repro_state.apply(lambda x:float(np.mean(np.asarray(x)==st)))
        bn=float(nfrac.median()); bt=float(tfrac.median()); obs=bt-bn
        for cond,g in CS[(CS.lineage==lin)&(CS.state==st)].groupby('condition'):
            gp=g.groupby('patient_uid').state_fraction.mean()
            bd=float(gp.median()); null=bd-bn
            ratio=null/obs if abs(obs)>1e-12 else np.nan
            state_summary.append([lin,st,cond,len(gp),bn,bt,obs,bd,null,ratio,np.sign(null)==np.sign(obs) if obs!=0 and null!=0 else False])
SUM=pd.DataFrame(summary,columns=['lineage','condition','n_normal_patients','baseline_normal_median_max_corr','primary_median_max_corr','degraded_normal_median_max_corr','fraction_observed_maxcorr_drop_reproduced','baseline_normal_median_margin','primary_median_margin','degraded_normal_median_margin','fraction_observed_margin_drop_reproduced','mean_label_stability','mean_median_retention'])
SS=pd.DataFrame(state_summary,columns=['lineage','state','condition','n_normal_patients','baseline_normal_median_fraction','primary_median_fraction','observed_primary_minus_normal','degraded_normal_median_fraction','null_degraded_minus_normal','null_to_observed_ratio','direction_concordant'])
PT=pd.DataFrame(paired_tests,columns=['lineage','condition','metric','n_pairs','median_baseline','median_degraded','median_paired_delta','W','p'])
PT['FDR']=np.nan
for (lin,metric),ix in PT.groupby(['lineage','metric']).groups.items():
    good=PT.loc[ix,'p'].notna()
    if good.any(): PT.loc[np.array(ix)[good.to_numpy()],'FDR']=multipletests(PT.loc[np.array(ix)[good.to_numpy()],'p'],method='fdr_bh')[1]
SUM.to_csv(OUT/'M11A3_degradation_calibration_summary_v0_1.tsv',sep='\t',index=False)
SS.to_csv(OUT/'M11A3_degradation_state_artifact_summary_v0_1.tsv',sep='\t',index=False)
PT.to_csv(OUT/'M11A3_degradation_paired_tests_v0_1.tsv',sep='\t',index=False)

manifest={
 'seed':SEED,'n_replicates':NREP,'fixed_retention_fractions':RETENTIONS,
 'centroid_definition':'mean across HBCA donor-by-state pseudobulk profiles on the frozen fine-state marker union',
 'wave_normalization':'log1p(count/total_UMI*10000) before Pearson correlation to state centroids',
 'degradation_null':'independent binomial UMI thinning; marker counts thinned gene-wise and all remaining UMIs thinned in aggregate',
 'depth_matched_null':'each normal cell assigned a randomly sampled primary-tumour total-UMI target within lineage; retention=min(1,target/normal_total)',
 'centroid_audit':centroid_audit
}
(OUT/'M11A3_degradation_calibration_manifest_v0_1.json').write_text(json.dumps(manifest,indent=2))
print('\nCALIBRATION SUMMARY\n',SUM.to_string(index=False),flush=True)
print('\nKEY STATE ARTIFACTS\n',SS[SS.state.isin(['FB1','FB2','FB3','FB4','VEV','VEA','VEAT','VEC','PV1','PV2','PV3','PV5'])].to_string(index=False),flush=True)
print('M11A3 PASS',time.ctime(),flush=True)
