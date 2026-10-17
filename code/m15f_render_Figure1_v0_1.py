#!/usr/bin/env python3
from pathlib import Path
import os
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib.gridspec import GridSpec

ROOT=Path(os.environ.get('BC_PROJECT_ROOT','.'))
SRC=ROOT/'StageM15F_FinalScientificFreeze_v0_1/source_data/main/Figure1'
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/figures/main'
OUT.mkdir(parents=True,exist_ok=True)

f1b=pd.read_csv(SRC/'Figure1B_normal_reference_scale_v0_1.tsv',sep='\t')
f1c=pd.read_csv(SRC/'Figure1C_normal_breast_lineage_support_v0_1.tsv',sep='\t')
f1d=pd.read_csv(SRC/'Figure1D_Reed_fine_state_donor_heldout_validation_v0_1.tsv',sep='\t')

def clean(ax):
 ax.spines[['top','right']].set_visible(False)
 ax.tick_params(labelsize=7)
def box(ax,xy,w,h,text,fc='#F8FAFC',fs=7):
 p=FancyBboxPatch(xy,w,h,boxstyle='round,pad=.015,rounding_size=.02',fc=fc,ec='#64748B',lw=.8)
 ax.add_patch(p);ax.text(xy[0]+w/2,xy[1]+h/2,text,ha='center',va='center',fontsize=fs)
def arrow(ax,a,b):
 ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=8,lw=.8,color='#64748B'))
def letter(ax,s):
 ax.text(-.08,1.05,s,transform=ax.transAxes,fontweight='bold',fontsize=11,va='top')

def A(ax):
 ax.set(xlim=(0,1),ylim=(0,1)); ax.axis('off')
 box(ax,(.03,.55),.29,.29,'Original human\ngenerating\nresources\n123',fc='#EEF4FB',fs=6.1)
 box(ax,(.38,.55),.24,.29,'Master\nregistry records\n137',fc='#F5F0FA',fs=6.1)
 box(ax,(.68,.53),.29,.33,'Overlap / reuse\nadjudication\nretained explicitly',fc='#FEF4E8',fs=6.0)
 arrow(ax,(.32,.695),(.38,.695));arrow(ax,(.62,.695),(.68,.695))
 ax.text(.50,.23,'Patient denominators are reported for each biological analysis;\nno resource-wide patient total is used as an inferential result.',ha='center',va='center',fontsize=6.2)
def B(ax):
 q=f1b[f1b.layer.isin(['CELLxGENE_primary_breast_normal','HBCA_anchor'])].copy()
 names={'CELLxGENE_primary_breast_normal':'Normal breast','HBCA_anchor':'Reed atlas anchor'}
 q['label']=q.layer.map(names); q['m']=q.n_cells/1e6
 y=np.arange(len(q)); ax.barh(y,q.m,height=.62)
 ax.set_yticks(y,q.label);ax.set_xscale('log');ax.set_xlabel('Cells (millions; log scale)',fontsize=7);clean(ax)
 for yi,v in zip(y,q.m): ax.text(v*1.08,yi,f'{v:.2f} M',va='center',fontsize=6.3)
 ax.invert_yaxis()
def C(ax):
 q=f1c[f1c.broad_lineage!='Other'].copy().sort_values('n_cells')
 ax.barh(q.broad_lineage,q.n_cells/1000,height=.65)
 ax.set_xlabel('Normal-breast cells (thousands)',fontsize=7);clean(ax)
 for y,v in enumerate(q.n_cells/1000):ax.text(v+8,y,f'{v:,.0f}',va='center',fontsize=6)
def D(ax):
 q=f1d.copy();y=np.arange(len(q));vals=100*q.LODO_accuracy
 ax.barh(y,vals,height=.62);ax.set_yticks(y,q.lineage);ax.set_xlim(90,100.5)
 ax.set_xlabel('Leave-one-donor-out accuracy (%)',fontsize=7);clean(ax)
 for yi,v in zip(y,vals):ax.text(v-.15,yi,f'{v:.2f}%',ha='right',va='center',fontsize=6.3)
 ax.invert_yaxis()
def E(ax):
 ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
 rows=[('Reference validation','Donor'),('State composition','Patient'),('Within-state expression','Patient × state'),('Ligand–receptor expression','Patient × state pair'),('Primary–LN comparison','Matched patient'),('Protein imaging','Patient')]
 y=.88
 for a,b in rows:
  ax.text(.03,y,a,ha='left',va='center',fontsize=6.5);box(ax,(.59,y-.055),.36,.11,b,fs=6.3);y-=.145
 ax.text(.03,.02,'Inference follows biological replicates, not cells or spots.',fontsize=6.2)

funcs=[A,B,C,D,E]
fig=plt.figure(figsize=(7.2,7.0));gs=GridSpec(3,2,figure=fig,height_ratios=[1,1,1.1],hspace=.55,wspace=.55)
axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1]),fig.add_subplot(gs[2,:])]
for ax,nm,fn in zip(axes,'ABCDE',funcs):fn(ax);letter(ax,nm)
for ext in ['png','pdf','svg']:
 fig.savefig(OUT/f'Figure1.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)
