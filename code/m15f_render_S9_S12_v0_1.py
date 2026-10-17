#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"
SRC=ST/"source_data/supplementary"; OUT=ST/"figures/supplementary"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"axes.labelsize":7.5,"xtick.labelsize":6.3,"ytick.labelsize":6.3,"legend.fontsize":6.0,"axes.linewidth":.7,"pdf.fonttype":42,"svg.fonttype":"none","savefig.transparent":False})
def clean(ax): ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
def save(fig,name):
 for ext in ["pdf","png","svg"]:fig.savefig(OUT/f"{name}.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
 plt.close(fig)

# ---------------- S9 ----------------
over=pd.read_csv(SRC/"FigureS9_IMC_overlap_v0_1.tsv",sep="\t")
test=pd.read_csv(ST/"tables/TableS9_IMC_gate_sensitivity.tsv",sep="\t")
prot=pd.read_csv(SRC/"FigureS9_CPTAC_stromal_protein_validation_v0_1.tsv",sep="\t")
fig,axs=plt.subplots(1,3,figsize=(11.2,3.8),gridspec_kw={"wspace":.52})
# a
ax=axs[0]; order=["Normal","DCIS","Tumor"];q=over.set_index("condition").loc[order].reset_index();x=np.arange(3)
ax.bar(x,q.overlap_percent,color="#93A4B8");ax.set_xticks(x,["Normal","DCIS","Invasive"]);ax.set_ylabel("Epithelial overlap in published\nperivascular gate (%)");clean(ax)
for i,r in q.iterrows():ax.text(i,r.overlap_percent+.25,f"{r.overlap_percent:.1f}%\n{int(r.n_cells):,} cells\n{int(r.n_patients)} patients",ha="center",va="bottom",fontsize=5.4)
ax.set_ylim(0,max(q.overlap_percent)*1.35);ax.text(-.16,1.04,"a",transform=ax.transAxes,fontweight="bold",fontsize=11)
# b: invasive-normal headline with expanded 15-test correction
ax=axs[1];gates=["original","exclude_is_epithelial","exclude_epi_and_marker_high"];glab={"original":"Original","exclude_is_epithelial":"Exclude epithelial","exclude_epi_and_marker_high":"Stricter CK/panCK"}
q=test[(test.contrast=="Tumor_vs_Normal") & test.marker.isin(["SMA","Caveolin-1"])].copy()
ys=[];labs=[]
for gi,g in enumerate(gates):
 for mi,m in enumerate(["SMA","Caveolin-1"]):
  r=q[(q.gate==g)&(q.marker==m)].iloc[0];pos=gi*2.4+mi*.75
  ax.scatter(r.delta_a_minus_b,pos,s=42,color="#B65A50" if m=="SMA" else "#4C78A8",edgecolor="white",linewidth=.4,zorder=3)
  ax.text(1.02,pos,f"FDR15 {r.FDR_within_gate_all15:.3g}; n={int(r.n_b)}/{int(r.n_a)}",va="center",fontsize=5.4)
  ys.append(pos);labs.append(f"{glab[g]} | {'α-SMA' if m=='SMA' else 'Caveolin-1'}")
ax.axvline(0,color="#64748B",ls="--",lw=.8);ax.set_yticks(ys,labs);ax.invert_yaxis();ax.set_xlabel("Invasive − normal patient median");clean(ax);ax.set_xlim(-2.55,2.35)
ax.text(-.16,1.04,"b",transform=ax.transAxes,fontweight="bold",fontsize=11)
# c
ax=axs[2];p=prot.sort_values("paired_median_delta").copy();yy=np.arange(len(p))
ax.scatter(p.paired_median_delta,yy,s=55,color="#7B6F9B",edgecolor="white",linewidth=.5);ax.axvline(0,color="#64748B",ls="--",lw=.8)
ax.set_yticks(yy,p.gene);ax.set_xlabel("Paired tumour − normal protein abundance");clean(ax)
for i,r in p.reset_index(drop=True).iterrows():ax.text(r.paired_median_delta+.08,i,f"FDR {r.paired_FDR:.3g}; n={int(r.n_paired)}",va="center",fontsize=5.7)
ax.text(-.16,1.04,"c",transform=ax.transAxes,fontweight="bold",fontsize=11)
fig.text(.5,.005,"Panel b displays invasive–normal effects using the expanded 15 marker-by-contrast FDR family; DCIS–normal and DCIS–invasive results are reported in Supplementary Table 9.",ha="center",fontsize=6.2)
save(fig,"FigureS9")

# ---------------- S12 ----------------
A=pd.read_csv(SRC/"FigureS12A_scanvi_heldout_state_accuracy_v0_1.tsv",sep="\t")
B=pd.read_csv(SRC/"FigureS12B_crossmapper_role_concordance_v0_1.tsv",sep="\t")
C=pd.read_csv(SRC/"FigureS12C_crossmapper_pivot_effects_v0_1.tsv",sep="\t")
fig,axs=plt.subplots(1,3,figsize=(13.4,4.2),gridspec_kw={"width_ratios":[1.15,.9,1.45],"wspace":.48})
# a
ax=axs[0];order=A.sort_values(["lineage","true_state"]).copy();y=np.arange(len(order))
cols={"Fibroblast":"#B07AA1","Endothelial":"#59A14F","Perivascular":"#F28E2B"}
ax.barh(y,100*order.state_accuracy,color=[cols[x] for x in order.lineage]);ax.set_yticks(y,order.true_state);ax.invert_yaxis();ax.set_xlim(60,100)
ax.set_xlabel("Held-out state accuracy (%)");clean(ax);ax.text(-.16,1.04,"a",transform=ax.transAxes,fontweight="bold",fontsize=11)
# b, only normal and invasive primary
ax=axs[1];bb=B[B.role.isin(["reference_normal","core_primary"])].copy();lins=["Fibroblast","Endothelial","Perivascular"];x=np.arange(3);w=.34
n=bb[bb.role=="reference_normal"].set_index("lineage").loc[lins]
p=bb[bb.role=="core_primary"].set_index("lineage").loc[lins]
ax.bar(x-w/2,100*n.exact_state_concordance,w,label="Normal",color="#4C78A8")
ax.bar(x+w/2,100*p.exact_state_concordance,w,label="Invasive",color="#ECA82C")
ax.set_xticks(x,["Fibroblast","Endothelial","Perivascular"],rotation=22,ha="right");ax.set_ylabel("Exact cross-mapper concordance (%)");ax.set_ylim(35,90);clean(ax);ax.legend(frameon=False)
for xx,v,nn in zip(x-w/2,100*n.exact_state_concordance,n.n_cells):ax.text(xx,v+1,f"n={int(nn):,}",ha="center",fontsize=5)
for xx,v,nn in zip(x+w/2,100*p.exact_state_concordance,p.n_cells):ax.text(xx,v+1,f"n={int(nn):,}",ha="center",fontsize=5)
ax.text(-.16,1.04,"b",transform=ax.transAxes,fontweight="bold",fontsize=11)
# c
ax=axs[2]
for lin,q in C.groupby("lineage"):
 ax.scatter(q.pivot_effect,q.delta_pivot,s=36,label=lin,color=cols[lin],alpha=.85)
 for r in q.itertuples():ax.text(r.pivot_effect+.08,r.delta_pivot+.03,r.state,fontsize=5.1)
lo=min(C.pivot_effect.min(),C.delta_pivot.min())-.4;hi=max(C.pivot_effect.max(),C.delta_pivot.max())+.4
ax.plot([lo,hi],[lo,hi],ls="--",lw=.8,color="#64748B");ax.axhline(0,lw=.5,color="#CBD5E1");ax.axvline(0,lw=.5,color="#CBD5E1")
ax.set(xlim=(lo,hi),ylim=(lo,hi));ax.set_xlabel("Historical pivot effect");ax.set_ylabel("scANVI pivot effect");clean(ax);ax.legend(frameon=False,fontsize=5.5)
ax.text(-.16,1.04,"c",transform=ax.transAxes,fontweight="bold",fontsize=11)
fig.text(.5,.01,"Panel b separates DCIS from invasive primary disease; DCIS contains only 25 fibroblast, 46 endothelial and 11 perivascular mapped cells and is not included in the invasive concordance bars.",ha="center",fontsize=6.0)
fig.subplots_adjust(bottom=.20)
save(fig,"FigureS12")
print("PASS S9/S12")
