#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import shutil, json
import numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM18_IntegratedRevision_v0_1"
M17=ROOT/"StageM17A_VascularSpatialClinicalPhospho_v0_1"
SRC5=ST/"source_data/main/Figure5"
SRCS=ST/"source_data/supplementary"
TAB=ST/"tables"; OUTM=ST/"figures/main"; OUTS=ST/"figures/supplementary"; AUD=ST/"audit"
for p in [SRC5,SRCS,TAB,OUTM,OUTS,AUD]: p.mkdir(parents=True,exist_ok=True)

# ---------- propagate source tables ----------
shutil.copy2(M17/"results/M17A4_IMC_murallike_patient_marker_medians_v0_1.tsv",
             SRC5/"Figure5C_IMC_murallike_patient_marker_medians_v0_1.tsv")
shutil.copy2(M17/"results/M17A4_IMC_murallike_marker_tests_v0_1.tsv",
             SRC5/"Figure5C_IMC_murallike_marker_tests_v0_1.tsv")
shutil.copy2(M17/"results/M17A3_CPTAC_module_phosphosite_tests_v0_1.tsv",
             SRCS/"FigureS13_CPTAC_phosphoproteomic_module_sites_v0_1.tsv")
shutil.copy2(M17/"results/M17A3_CPTAC_module_summary_v0_1.tsv",
             SRCS/"FigureS13_CPTAC_phosphoproteomic_module_summary_v0_1.tsv")
shutil.copy2(M17/"results/M17A3_CPTAC_module_phosphosite_tests_v0_1.tsv",
             TAB/"TableS18_CPTAC_phosphoproteomic_module_v0_1.tsv")

# Append myofibroblast-defined IMC tests to Supplementary Table 9.
old=pd.read_csv(ST/"tables/TableS9_IMC_gate_sensitivity.tsv",sep="\t")
myo=pd.read_csv(M17/"results/M17A4_IMC_murallike_marker_tests_v0_1.tsv",sep="\t")
add=pd.DataFrame({
 "gate":"mural_like_cellPhenotype",
 "marker":myo.marker,
 "contrast":myo.contrast,
 "n_a":myo.n_a,"n_b":myo.n_b,
 "median_a":myo.median_a,"median_b":myo.median_b,
 "delta_a_minus_b":myo.delta_a_minus_b,"U":myo.U,"p":myo.p,
 "rank_biserial_a_over_b":myo.rank_biserial_a_over_b,
 "bootstrap_ci95_low":np.nan,"bootstrap_ci95_high":np.nan,
 "FDR_within_gate_all15":myo.FDR_all15,
 "FDR_within_gate_contrast5":myo.FDR_within_contrast5
})
new=pd.concat([old,add],ignore_index=True)
new.to_csv(ST/"tables/TableS9_IMC_gate_sensitivity.tsv",sep="\t",index=False)
new.to_csv(SRCS/"FigureS9_IMC_threeway_gate_tests_v0_1.tsv",sep="\t",index=False)

# ---------- Figure 5 ----------
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"axes.labelsize":7.5,
"xtick.labelsize":6.3,"ytick.labelsize":6.3,"legend.fontsize":6.0,
"axes.linewidth":0.7,"svg.fonttype":"none","pdf.fonttype":42,"savefig.transparent":False})
C={"Normal":"#4C78A8","DCIS":"#B07AA1","Tumor":"#ECA82C","up":"#B65A50","down":"#4C78A8","gray":"#6B7280"}
va=pd.read_csv(SRC5/"Figure5A_independent_normal_reference_validation_v0_1.tsv",sep="\t")
bg=pd.read_csv(SRC5/"Figure5B_matched_background_directional_consistency_v0_1.tsv",sep="\t")
imcp=pd.read_csv(SRC5/"Figure5C_IMC_murallike_patient_marker_medians_v0_1.tsv",sep="\t")
imct=pd.read_csv(SRC5/"Figure5C_IMC_murallike_marker_tests_v0_1.tsv",sep="\t")
cpt=pd.read_csv(SRC5/"Figure5D_paired_CPTAC_stromal_proteins_v0_1.tsv",sep="\t")
s100=pd.read_csv(SRC5/"Figure5F_prespecified_S100A4_ERBB2_evaluation_v0_1.tsv",sep="\t")

def clean(ax):
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
def letter(ax,s):
    ax.text(-.12,1.08,s,transform=ax.transAxes,fontweight="bold",fontsize=10,va="top")
def save(fig,stem):
    for ext in ["svg","pdf","png"]:
        fig.savefig(str(stem)+"."+ext,dpi=300 if ext=="png" else None,bbox_inches="tight")
    plt.close(fig)
def f5a(ax):
    q=va[va.tumour_group.isin(["all_primary","ER+","TNBC"])].copy()
    lab={"all_primary":"All primary","ER+":"ER+","TNBC":"TNBC"};x=np.arange(len(q));w=.28
    ax.bar(x-w/2,q.direction_concordant,width=w,color="#9CA3AF",label="Direction")
    ax.bar(x+w/2,q.FDR05_concordant,width=w,color="#4F8A70",label="Direction + FDR<0.05")
    ax.set_xticks(x,[lab[x] for x in q.tumour_group]);ax.set_ylim(0,36);ax.set_ylabel("Interactions (of 34)");clean(ax);ax.legend(frameon=False)
def f5b(ax):
    q=bg.copy();lab={"all_noncore":"All non-prioritized","matched_noncore":"Matched background","core":"Prioritized"}
    y=np.arange(len(q));vals=100*q.consistency_rate
    ax.barh(y,vals,color=["#BFC5CC","#8FA3B8","#4F8A70"],height=.62)
    ax.set_yticks(y,[lab[x] for x in q.set]);ax.set_xlim(75,102);ax.set_xlabel("Directional consistency (%)");clean(ax)
    for yi,v in enumerate(vals):ax.text(v-.4,yi,f"{v:.1f}%",ha="right",va="center",fontsize=6,color="white" if v>90 else "#111827")
def f5c(ax):
    q=imcp[imcp.eligible_ge20].copy()
    groups=["Normal","DCIS","Tumor"]; pos=np.arange(3)
    data=[q.loc[q.condition==g,"SMA"].dropna().to_numpy(float) for g in groups]
    bp=ax.boxplot(data,positions=pos,widths=.58,patch_artist=True,showfliers=False,
                  medianprops={"color":"#111827","linewidth":1.0},
                  whiskerprops={"linewidth":.7},capprops={"linewidth":.7})
    for patch,g in zip(bp["boxes"],groups):
        patch.set_facecolor(C[g]);patch.set_alpha(.55);patch.set_edgecolor(C[g])
    rng=np.random.default_rng(1813)
    for x,v,g in zip(pos,data,groups):
        vv=v if len(v)<=140 else rng.choice(v,140,replace=False)
        ax.scatter(x+rng.normal(0,.055,len(vv)),vv,s=5,alpha=.27,color=C[g],edgecolor="none",zorder=2)
        ax.text(x,-.12,f"n={len(v)}",ha="center",va="top",fontsize=5,transform=ax.get_xaxis_transform())
    ax.set_xticks(pos,["Normal","DCIS","Invasive"])
    ax.set_ylabel("Patient median α-SMA IMC intensity");clean(ax)
    tt=imct[imct.marker.eq("SMA")].set_index("contrast")
    ax.text(.02,.98,
      "Myofibroblast-defined compartment\n"+
      f"FDR: N–DCIS {tt.loc['DCIS_vs_Normal','FDR_all15']:.3g}; "+
      f"N–Inv {tt.loc['Tumor_vs_Normal','FDR_all15']:.3g}; "+
      f"DCIS–Inv {tt.loc['DCIS_vs_Tumor','FDR_all15']:.3g}",
      transform=ax.transAxes,ha="left",va="top",fontsize=5.0)
def f5d(ax):
    q=cpt.copy();y=np.arange(len(q))
    ax.barh(y,q.CPTAC_paired_protein_delta,color="#7B6F9B",height=.6)
    ax.set_yticks(y,q.gene);ax.set_xlabel("Paired tumour − normal protein difference");clean(ax);ax.invert_yaxis()
    for yi,(v,f) in enumerate(zip(q.CPTAC_paired_protein_delta,q.CPTAC_paired_protein_FDR)):
        ax.text(v+.04,yi,f"FDR {f:.2g}",va="center",fontsize=5.4)
def f5e(ax):
    wanted=["all","ERBB2_pos=True","PAM50=HER2","PAM50=Basal","PAM50=Luminal A","PAM50=Luminal B"]
    q=s100[s100.group.isin(wanted)].copy()
    q["label"]=q.group.replace({"all":"All","ERBB2_pos=True":"ERBB2+","PAM50=HER2":"PAM50 HER2","PAM50=Basal":"PAM50 Basal","PAM50=Luminal A":"Luminal A","PAM50=Luminal B":"Luminal B"})
    y=np.arange(len(q));ax.scatter(q.median_log2OR,y,s=28,color="#6B7280");ax.axvline(0,color="#64748B",lw=.8,ls="--")
    ax.set_yticks(y,q.label);ax.set_xlabel("Median spatial log2 odds ratio");clean(ax);ax.invert_yaxis();ax.set_xlim(-.68,.03)
    for yi,(x,n) in enumerate(zip(q.median_log2OR,q.n_or)):
        ax.annotate(f"n={int(n)}",(x,yi),xytext=(6,0),textcoords="offset points",ha="left",va="center",fontsize=5.3)

fig=plt.figure(figsize=(7.2,6.8));gs=GridSpec(3,2,figure=fig,height_ratios=[1,1,1],hspace=.74,wspace=.72)
axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1]),fig.add_subplot(gs[2,:])]
for ax,nm,fn in zip(axes,"ABCDE",[f5a,f5b,f5c,f5d,f5e]): fn(ax);letter(ax,nm)
save(fig,OUTM/"Figure5")

# ---------- Supplementary Figure 13: phosphosite module ----------
site=pd.read_csv(SRCS/"FigureS13_CPTAC_phosphoproteomic_module_sites_v0_1.tsv",sep="\t")
summ=pd.read_csv(SRCS/"FigureS13_CPTAC_phosphoproteomic_module_summary_v0_1.tsv",sep="\t")
q=site[site.n_protein_adjusted>=5].copy()
fig=plt.figure(figsize=(7.2,3.35));gs=GridSpec(1,2,figure=fig,wspace=.48)
ax=fig.add_subplot(gs[0,0])
sig=q.protein_adjusted_module_family_FDR<.05
ax.scatter(q.raw_median_delta,q.protein_adjusted_median_delta,s=16,alpha=.55,color="#9CA3AF",label="Testable site")
ax.scatter(q.loc[sig,"raw_median_delta"],q.loc[sig,"protein_adjusted_median_delta"],s=22,alpha=.75,color="#7B6F9B",label="Adjusted FDR<0.05")
lo=float(np.nanmin([q.raw_median_delta.min(),q.protein_adjusted_median_delta.min()])); hi=float(np.nanmax([q.raw_median_delta.max(),q.protein_adjusted_median_delta.max()]))
ax.plot([lo,hi],[lo,hi],ls="--",lw=.7,color="#64748B")
ax.axhline(0,lw=.55,color="#CBD5E1");ax.axvline(0,lw=.55,color="#CBD5E1")
ax.set_xlabel("Raw paired phosphosite Δ");ax.set_ylabel("Protein-adjusted phosphosite Δ");clean(ax);letter(ax,"A")
# label representative prespecified sites
labels={"NP_001341945.1:s641":"PDGFRB S641","NP_001339623.1:y570":"PTK2 Y570",
        "NP_001074324.1:s533":"PXN S533","NP_003364.1:s721":"VCL S721",
        "NP_006280.3:s446":"TLN1 S446","NP_005397.1:s1341":"ROCK1 S1341",
        "NP_002464.1:s1943":"MYH9 S1943","NP_005498.1:s41":"CFL1 S41",
        "NP_001269026.1:s251":"YAP1 S251","NP_003205.2:s148":"TEAD3 S148"}
for r in q[q.phosphosite.isin(labels)].itertuples():
    ax.annotate(labels[r.phosphosite],(r.raw_median_delta,r.protein_adjusted_median_delta),
                xytext=(3,3),textcoords="offset points",fontsize=4.5)
ax.legend(frameon=False,loc="lower right")

ax=fig.add_subplot(gs[0,1])
order=["ECM_receptor_context","focal_adhesion","Rho_actomyosin","cofilin_axis","YAP_TAZ_mechanotransduction"]
s=summ.set_index("module").reindex(order).reset_index()
y=np.arange(len(s))
ax.barh(y,s.genes_with_testable_adjusted_site,color="#BFC5CC",height=.62,label="Genes with testable adjusted site")
ax.barh(y,s.genes_with_adjusted_module_FDR05_site,color="#4F8A70",height=.38,label="Genes with ≥1 adjusted FDR<0.05 site")
lab=["ECM/receptor","Focal adhesion","Rho/actomyosin","Cofilin axis","YAP/TAZ"]
ax.set_yticks(y,lab);ax.set_xlabel("Genes in predefined module");clean(ax);ax.invert_yaxis();letter(ax,"B");ax.legend(frameon=False,fontsize=5.4)
save(fig,OUTS/"FigureS13")

gate={
 "status":"M18_M17A_SOURCE_AND_FIGURE_PROPAGATION_COMPLETE",
 "Figure5C_definition":"deposited Myofibroblasts + Myofibroblasts PDPN+; alpha-SMA only",
 "historical_vessel_contact_IMC_retained":"Supplementary Table 9 / Figure S9 with semantic relabeling required in text",
 "TableS18_rows":int(len(site)),
 "FigureS13_testable_sites":int(q.shape[0]),
 "FigureS13_adjusted_FDR05_sites":int(sig.sum()),
}
(AUD/"M18_M17A_source_figure_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
