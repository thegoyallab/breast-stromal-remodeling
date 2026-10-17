#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM15F_FinalScientificFreeze_v0_1"; SRC=ST/"source_data/supplementary";OUT=ST/"figures/supplementary"
p=pd.read_csv(SRC/"FigureS4_patient_confidence_v0_1.tsv",sep="\t")
p=p[p.stage06a_role.isin(["reference_normal","core_primary"])].copy()
cal=pd.read_csv(SRC/"FigureS4_degradation_calibration_v0_1.tsv",sep="\t");cal=cal[cal.condition=="tumor_depth_matched"]
null=pd.read_csv(SRC/"FigureS4_full_state_degradation_null_v0_1.tsv",sep="\t")
order=["FB1","FB2","FB3","FB4","VEV","VEA","VEAT","VEC","LE1","LE2","PV1","PV2","PV3","PV4","PV5"]
null=null.set_index("state").loc[order].reset_index()
lins=["Fibroblast","Endothelial","Perivascular"];C={"N":"#4C78A8","T":"#ECA82C"}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"pdf.fonttype":42,"svg.fonttype":"none"})
fig,axs=plt.subplots(2,2,figsize=(10.8,7.2),gridspec_kw={"wspace":.32,"hspace":.40})
for ax,metric,ylabel,letter in [(axs[0,0],"median_max_corr","Patient median maximum correlation","a"),(axs[0,1],"median_margin","Patient median assignment margin","b")]:
 pos=[];labs=[]
 for i,lin in enumerate(lins):
  for j,(role,col) in enumerate([("reference_normal",C["N"]),("core_primary",C["T"])]):
   v=p[(p.lineage==lin)&(p.stage06a_role==role)][metric].to_numpy();x=i*2.2+j*.72
   bp=ax.boxplot([v],positions=[x],widths=.48,patch_artist=True,showfliers=False,medianprops={"color":"#111827","linewidth":1.2},whiskerprops={"color":"#64748B"},capprops={"color":"#64748B"})
   bp["boxes"][0].set_facecolor(col);bp["boxes"][0].set_alpha(.6);bp["boxes"][0].set_edgecolor(col)
   jitter=np.linspace(-.14,.14,len(v)) if len(v)>1 else np.array([0.]);ax.scatter(x+jitter,v,s=8,color=col,alpha=.65,zorder=3)
  pos.append(i*2.2+.36);labs.append(lin)
 ax.set_xticks(pos,labs,rotation=20,ha="right");ax.set_ylabel(ylabel);ax.spines[["top","right"]].set_visible(False);ax.text(-.13,1.04,letter,transform=ax.transAxes,fontweight="bold",fontsize=11)
axs[0,0].scatter([],[],color=C["N"],label="Normal");axs[0,0].scatter([],[],color=C["T"],label="Invasive primary");axs[0,0].legend(frameon=False,loc="lower left",fontsize=6.5)
ax=axs[1,0];q=cal.set_index("lineage").loc[lins];x=np.arange(3);w=.24
ax.bar(x-w,q.fraction_observed_maxcorr_drop_reproduced*100,width=w,label="Max-corr drop reproduced")
ax.bar(x,q.fraction_observed_margin_drop_reproduced*100,width=w,label="Margin drop reproduced")
ax.bar(x+w,q.mean_label_stability*100,width=w,label="Label stability")
ax.set_xticks(x,lins,rotation=20,ha="right");ax.set_ylabel("Percent");ax.set_ylim(0,105);ax.spines[["top","right"]].set_visible(False);ax.legend(frameon=False,fontsize=6,loc="upper left");ax.text(-.13,1.04,"c",transform=ax.transAxes,fontweight="bold",fontsize=11)
ax=axs[1,1];x=np.arange(len(null));w=.36
ax.bar(x-w/2,null.observed_primary_minus_normal,w,label="Observed invasive − normal")
ax.bar(x+w/2,null.null_degraded_minus_normal,w,label="Depth-matched degraded normal − normal")
ax.axhline(0,lw=.7,color="#64748B");ax.set_xticks(x);ax.set_xticklabels(null.state,rotation=45,ha="right");ax.set_ylabel("Change in mapped state fraction");ax.legend(frameon=False,fontsize=5.7,loc="upper left")
j=order.index("VEAT");ratio=float(null.loc[null.state=="VEAT","null_to_observed_ratio"].iloc[0]);yv=float(null.loc[null.state=="VEAT","null_degraded_minus_normal"].iloc[0])
ax.annotate(f"VEAT null/observed = {ratio:.2f}",xy=(j+w/2,yv),xytext=(j+1.1,.12),arrowprops=dict(arrowstyle="->",lw=.7),fontsize=6.2)
ax.spines[["top","right"]].set_visible(False);ax.text(-.13,1.04,"d",transform=ax.transAxes,fontweight="bold",fontsize=11)
fig.text(.5,.008,"Panels a–b exclude DCIS from invasive-primary summaries. Panel c: tumour-depth-matched UMI-thinning calibration with the corrected invasive denominator. Panel d: complete 15-state degradation null.",ha="center",fontsize=6.4)
fig.subplots_adjust(bottom=.12)
for ext in ["pdf","png","svg"]:fig.savefig(OUT/f"FigureS4.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
plt.close(fig)
print("PASS S4")
