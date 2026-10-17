from pathlib import Path
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
SRC=ROOT/'StageM15F_FinalScientificFreeze_v0_1/source_data/supplementary'
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/figures/supplementary'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,'legend.fontsize':7,'axes.linewidth':0.8,'svg.fonttype':'none','pdf.fonttype':42})
metrics=['FB2','FB3','FB4','VEA','VEC','LE1','LE2','PV1','PV2','PV3','PV5']
sens=pd.read_csv(SRC/'FigureS3_cell_threshold_and_cohort_sensitivity.tsv',sep='\t')
mapq=pd.read_csv(SRC/'FigureS3_mapping_quality_exclusion_sensitivity.tsv',sep='\t')
fig,axs=plt.subplots(1,2,figsize=(7.2,3.8),gridspec_kw={'width_ratios':[1.10,1.0],'wspace':.62})
ax=axs[0]
cols=[('all_cohorts',1,'All | 1'),('all_cohorts',20,'All | 20'),('all_cohorts',30,'All | 30'),('GSE161529_only',10,'GSE161529 | 10'),('GSE161529_only',20,'GSE161529 | 20'),('GSE161529_only',30,'GSE161529 | 30')]
arr=np.full((len(metrics),len(cols)),np.nan)
z=sens[sens.contrast.eq('Normal_vs_Primary')]
for i,m in enumerate(metrics):
  for j,(sc,f,lab) in enumerate(cols):
    q=z[(z.metric==m)&(z.scope==sc)&(z.min_cells==f)]
    if len(q): arr[i,j]=q.iloc[0].delta_B_minus_A
im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-.65,vmax=.65)
ax.set_yticks(range(len(metrics)),metrics); ax.set_xticks(range(len(cols)),[c[2] for c in cols],rotation=43,ha='right')
ax.set_title('Cell-threshold and cohort sensitivity',fontsize=8.5)
cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.03); cb.set_label('Primary − normal state fraction')
ax.text(-.14,1.04,'a',transform=ax.transAxes,fontweight='bold',fontsize=11)

ax=axs[1]
q=mapq[mapq.min_cells.eq(1)].copy()
y=np.arange(len(metrics))
allv=[]; exv=[]
for m in metrics:
  qa=q[(q.metric==m)&(~q.exclude_accuracy_lt0p90)]
  qe=q[(q.metric==m)&(q.exclude_accuracy_lt0p90)]
  allv.append(qa.delta.iloc[0] if len(qa) else np.nan); exv.append(qe.delta.iloc[0] if len(qe) else np.nan)
for yi,a,e in zip(y,allv,exv):
  ax.plot([a,e],[yi,yi],color='#CBD5E1',lw=1.2,zorder=1)
ax.scatter(allv,y,s=30,color='#4C78A8',label='All patients',zorder=2)
ax.scatter(exv,y,s=30,facecolors='white',edgecolors='#D95F5F',linewidths=1.2,label='Exclude mapping accuracy <0.90',zorder=3)
ax.axvline(0,color='#64748B',ls='--',lw=.8)
ax.set_yticks(y,metrics); ax.invert_yaxis(); ax.set_xlabel('Primary − normal state fraction'); ax.set_title('Mapping-quality exclusion sensitivity',fontsize=8.5)
ax.legend(frameon=False,loc='lower right')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.text(-.14,1.04,'b',transform=ax.transAxes,fontweight='bold',fontsize=11)
for ext in ['png','pdf','svg']:
  kwargs={'dpi':300} if ext=='png' else {}
  fig.savefig(OUT/f'FigureS3.{ext}',bbox_inches='tight',**kwargs)
plt.close(fig)
