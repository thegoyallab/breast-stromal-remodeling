#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
SRC=ST/"source_data/main/Figure5"
OUT=ST/"figures/main"
OUT.mkdir(parents=True,exist_ok=True)

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"axes.labelsize":7.5,
"xtick.labelsize":6.3,"ytick.labelsize":6.3,"legend.fontsize":6.0,
"axes.linewidth":0.7,"svg.fonttype":"none","pdf.fonttype":42,"savefig.transparent":False})
C={"Normal":"#4C78A8","DCIS":"#B07AA1","Tumor":"#ECA82C","up":"#B65A50","down":"#4C78A8","gray":"#6B7280"}

va=pd.read_csv(SRC/"Figure5A_independent_normal_reference_validation_v0_1.tsv",sep="\t")
bg=pd.read_csv(SRC/"Figure5B_matched_background_directional_consistency_v0_1.tsv",sep="\t")
imc=pd.read_csv(SRC/"Figure5C_IMC_contractile_perivascular_validation_v0_1.tsv",sep="\t")
imcp=pd.read_csv(SRC/"Figure5C_IMC_patient_marker_medians_threeway_v0_1.tsv",sep="\t")
cpt=pd.read_csv(SRC/"Figure5D_paired_CPTAC_stromal_proteins_v0_1.tsv",sep="\t")
s100=pd.read_csv(SRC/"Figure5F_prespecified_S100A4_ERBB2_evaluation_v0_1.tsv",sep="\t")

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
    q=imcp[(imcp.gate=="exclude_is_epithelial") & imcp.eligible_ge20].copy()
    groups=["Normal","DCIS","Tumor"];markers=["SMA","Caveolin-1"]
    bases=[0,4]; positions=[];data=[];colors=[];labels=[]
    for bi,m in zip(bases,markers):
        for gi,g in enumerate(groups):
            v=q.loc[q.condition==g,m].dropna().to_numpy(float)
            positions.append(bi+gi);data.append(v);colors.append(C[g]);labels.append(g)
    bp=ax.boxplot(data,positions=positions,widths=.64,patch_artist=True,showfliers=False,
                  medianprops={"color":"#111827","linewidth":1.0},
                  whiskerprops={"linewidth":.7},capprops={"linewidth":.7})
    for patch,col in zip(bp["boxes"],colors):
        patch.set_facecolor(col);patch.set_alpha(.55);patch.set_edgecolor(col)
    rng=np.random.default_rng(1513)
    for pos,v,col in zip(positions,data,colors):
        # show all DCIS/normal; deterministic 120-point sample for large invasive group
        vv=v if len(v)<=120 else rng.choice(v,120,replace=False)
        ax.scatter(pos+rng.normal(0,.055,len(vv)),vv,s=5,alpha=.28,color=col,edgecolor="none",zorder=2)
        ax.text(pos,-.13,f"n={len(v)}",ha="center",va="top",fontsize=4.8,transform=ax.get_xaxis_transform())
    ax.set_xticks([0,1,2,4,5,6],["N","DCIS","Inv.","N","DCIS","Inv."])
    ax.text(1,-.25,"α-SMA",ha="center",va="top",fontsize=6.4,transform=ax.get_xaxis_transform())
    ax.text(5,-.25,"Caveolin-1",ha="center",va="top",fontsize=6.4,transform=ax.get_xaxis_transform())
    ax.set_ylabel("Patient median IMC intensity");clean(ax)
    # expanded 15-test family annotations
    test=imc[(imc.gate=="exclude_is_epithelial") & imc.marker.isin(markers)].copy()
    sm=test[test.marker=="SMA"].set_index("contrast").FDR_within_gate_all15
    cv=test[test.marker=="Caveolin-1"].set_index("contrast").FDR_within_gate_all15
    ax.text(.01,.98,f"α-SMA FDR: N–DCIS {sm['DCIS_vs_Normal']:.3g}; N–Inv {sm['Tumor_vs_Normal']:.3g}; DCIS–Inv {sm['DCIS_vs_Tumor']:.3g}",
            transform=ax.transAxes,ha="left",va="top",fontsize=4.6)
    ax.text(.01,.90,f"CAV1 FDR: N–DCIS {cv['DCIS_vs_Normal']:.3g}; N–Inv {cv['Tumor_vs_Normal']:.3g}; DCIS–Inv {cv['DCIS_vs_Tumor']:.3g}",
            transform=ax.transAxes,ha="left",va="top",fontsize=4.6)

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

funcs=[f5a,f5b,f5c,f5d,f5e]
fig=plt.figure(figsize=(7.2,6.8));gs=GridSpec(3,2,figure=fig,height_ratios=[1,1,1],hspace=.74,wspace=.72)
axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1]),fig.add_subplot(gs[2,:])]
for ax,nm,fn in zip(axes,"ABCDE",funcs):fn(ax);letter(ax,nm)
save(fig,OUT/"Figure5")
print("rendered",OUT/"Figure5")
