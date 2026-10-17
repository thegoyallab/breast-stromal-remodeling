#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
BR=ROOT/'StageM11_PostClaudeScientificRevision_v0_1'
M02=ROOT/'StageM02_FigureSourceEvidencePackage_v0_2'
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/figures/work'
for x in ['Figure2']: (OUT/x).mkdir(parents=True,exist_ok=True)

plt.rcParams.update({
 'font.family':'DejaVu Sans','font.size':7.0,'axes.labelsize':7.4,'xtick.labelsize':6.2,'ytick.labelsize':6.2,
 'legend.fontsize':5.8,'axes.linewidth':0.7,'xtick.major.width':0.6,'ytick.major.width':0.6,
 'xtick.major.size':3,'ytick.major.size':3,'svg.fonttype':'none','pdf.fonttype':42,'savefig.transparent':False
})
C={'Normal':'#4C78A8','Pre-neoplastic':'#72B7B2','Primary':'#ECA82C','LN':'#D95F5F',
   'Fibroblast':'#B07AA1','Endothelial':'#59A14F','Perivascular':'#F28E2B',
   'up':'#B65A50','down':'#4C78A8','gray':'#6B7280','light':'#E5E7EB'}

def clean(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
def letter(ax,s):
    ax.text(-.12,1.07,s,transform=ax.transAxes,fontweight='bold',fontsize=10,va='top')
def save(fig,stem):
    for ext in ['svg','pdf','png']:
        fig.savefig(str(stem)+'.'+ext,dpi=300 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)
def save_panel(fig,stem):
    for ext in ['svg','png']:
        fig.savefig(str(stem)+'.'+ext,dpi=300 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)
def box(ax,xy,w,h,text,fc='#F8FAFC',ec='#64748B',fs=7):
    p=FancyBboxPatch(xy,w,h,boxstyle='round,pad=0.02,rounding_size=0.02',fc=fc,ec=ec,lw=.8)
    ax.add_patch(p); ax.text(xy[0]+w/2,xy[1]+h/2,text,ha='center',va='center',fontsize=fs)
def arrow(ax,a,b,rad=0):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=8,lw=.8,color='#64748B',connectionstyle=f'arc3,rad={rad}'))

# ---------------- FIGURE 2 ----------------
raw=pd.read_csv(ROOT/'StageM15F_FinalScientificFreeze_v0_1/source_data/main/Figure2/patient_state_fraction_summary.tsv',sep='\t')
piv=pd.read_csv(ROOT/'StageM13_PostIndependentReview_v0_1/results/M13A_compositional_pivot_state_tests_v0_1.tsv',sep='\t')
glob=pd.read_csv(ROOT/'StageM13_PostIndependentReview_v0_1/results/M13A_compositional_global_PERMANOVA_v0_1.tsv',sep='\t')
role_order=['reference_normal','reference_preneoplastic','core_primary','linked_progression_LN']

def f2a(ax):
    ax.axis('off'); ax.set(xlim=(0,1),ylim=(0,1))
    box(ax,(.02,.61),.18,.18,'Normal\nbreast',fc='#EDF4FA',fs=7)
    box(ax,(.27,.58),.26,.24,'BRCA1-associated\nprophylactic /\npre-neoplastic',fc='#ECF7F5',fs=6.1)
    box(ax,(.61,.61),.18,.18,'Primary\ntumour',fc='#FFF4DD',fs=7)
    box(ax,(.72,.18),.24,.18,'Tumour-involved\nlinked LN',fc='#FCEDED',fs=6.5)
    arrow(ax,(.70,.61),(.82,.36),rad=.06)
    ax.text(.43,.91,'Cross-sectional contexts — not a serial trajectory',ha='center',fontsize=5.5,color='#475569')

def f2b(ax):
    lins=['Fibroblast','Endothelial','Perivascular']; offs=[-.22,0,.22]; x=np.arange(4)
    for j,L in enumerate(lins):
        q=raw[raw.lineage==L]
        for i,r in enumerate(role_order):
            v=q.loc[q.role==r,'median_max_corr'].dropna().values
            if not len(v): continue
            pos=x[i]+offs[j]; rng=np.random.default_rng(3+j*10+i)
            ax.scatter(pos+rng.normal(0,.035,len(v)),v,s=8,alpha=.45,color=C[L],edgecolor='none')
            med=np.median(v); q1,q3=np.quantile(v,[.25,.75])
            ax.plot([pos,pos],[q1,q3],color=C[L],lw=2.2); ax.scatter([pos],[med],s=20,color=C[L],edgecolor='white',lw=.4,zorder=3)
    ax.set_xticks(x,['Normal','Pre-neoplastic','Primary','Linked LN'],rotation=22,ha='right')
    ax.set_ylabel('Similarity to normal reference'); ax.set_ylim(.25,.78); clean(ax)
    ax.legend([Line2D([0],[0],marker='o',ls='',color=C[L]) for L in lins],lins,frameon=False,ncol=3,loc='lower left',bbox_to_anchor=(0,1.01))

def f2_state(ax,lineage,states):
    q=raw[(raw.lineage==lineage)&(raw.n_cells>=20)&raw.role.isin(['reference_normal','core_primary'])].copy()
    n0=int((q.role=='reference_normal').sum()); n1=int((q.role=='core_primary').sum())
    status_map={
      ('Fibroblast','FB1'):'mapper-robust',('Fibroblast','FB2'):'mapper-robust',
      ('Fibroblast','FB3'):'mapper-sensitive',('Fibroblast','FB4'):'mapper-sensitive',
      ('Endothelial','VEV'):'mapper-sensitive',('Endothelial','VEA'):'mapper-sensitive',
      ('Endothelial','VEAT'):'mapper-sensitive',('Endothelial','VEC'):'mapper-sensitive',
      ('Endothelial','LE1'):'mapper-robust',('Endothelial','LE2'):'mapper-robust',
      ('Perivascular','PV1'):'mapper-robust',('Perivascular','PV2'):'mapper-robust',
      ('Perivascular','PV3'):'mapper-robust',('Perivascular','PV4'):'method-sensitive',
      ('Perivascular','PV5'):'secondary†'}
    y=np.arange(len(states))
    for i,st in enumerate(states):
        a=q.loc[q.role=='reference_normal',st].to_numpy(float)
        b=q.loc[q.role=='core_primary',st].to_numpy(float)
        n=float(np.median(a)); p=float(np.median(b))
        pp=piv[(piv.scope=='all_cohorts')&(piv.lineage==lineage)&(piv.state==st)].iloc[0]
        f=float(pp.FDR_within_lineage_scope)
        ax.plot([n,p],[i,i],color='#CBD5E1',lw=2)
        ax.scatter(n,i,s=28,color=C['Normal'],zorder=3); ax.scatter(p,i,s=28,color=C['Primary'],zorder=3)
        ft='<0.001' if f<.001 else f'{f:.3f}'
        ax.text(.69,i,f"{status_map[(lineage,st)]}\nFDR {ft}",va='center',fontsize=4.4,color='#475569')
    ax.set_yticks(y,states); ax.invert_yaxis(); ax.set_xlim(-.02,.98)
    ax.set_xlabel('Median patient state fraction'); clean(ax)
    gp=glob[(glob.scope=='all_cohorts')&(glob.lineage==lineage)].iloc[0]
    ax.text(0,1.02,f'Global P={gp.PERMANOVA_p:.1e} in both mappers; n={n0}/{n1}',transform=ax.transAxes,fontsize=5.2)
    if lineage=='Fibroblast':
        note='FB1+FB3 vs FB2+FB4 axis retained by both mappers'
    elif lineage=='Endothelial':
        note='LE1/LE2 robust; vascular fine-state allocation mapper-sensitive'
    else:
        note='† PV5 not retained after strict reference acceptance'
    ax.text(0,-.18,note,transform=ax.transAxes,fontsize=4.7,color='#475569')
    ax.legend([Line2D([0],[0],marker='o',ls='',color=C['Normal']),Line2D([0],[0],marker='o',ls='',color=C['Primary'])],
              ['Normal','Primary'],frameon=False,loc='upper center',bbox_to_anchor=(.52,-.27),ncol=2)

def f2c(ax): f2_state(ax,'Fibroblast',['FB1','FB2','FB3','FB4'])
def f2d(ax): f2_state(ax,'Endothelial',['VEV','VEA','VEAT','VEC','LE1','LE2'])
def f2e(ax): f2_state(ax,'Perivascular',['PV1','PV2','PV3','PV4','PV5'])
def f2f(ax):
    vals=[]
    for L in ['Fibroblast','Endothelial','Perivascular']:
        q=raw[(raw.lineage==L)&raw.role.isin(['core_primary','linked_progression_LN'])]
        pv=q.pivot_table(index='patient_uid',columns='role',values='median_max_corr',aggfunc='first').dropna()
        for dd in (pv.linked_progression_LN-pv.core_primary): vals.append((L,dd))
    for i,L in enumerate(['Fibroblast','Endothelial','Perivascular']):
        v=np.array([dd for ll,dd in vals if ll==L]); rng=np.random.default_rng(20+i)
        ax.scatter(np.full(len(v),i)+rng.normal(0,.055,len(v)),v,s=17,color=C[L],alpha=.75,edgecolor='white',lw=.3)
        ax.plot([i-.16,i+.16],[np.median(v),np.median(v)],color='#111827',lw=1.6)
        ax.text(i,ax.get_ylim()[1] if i else 0,'',fontsize=1)
    ax.axhline(0,color='#64748B',lw=.8,ls='--'); ax.set_xticks(range(3),['Fibroblast','Endothelial','Perivascular'],rotation=15)
    ax.set_ylabel('Linked LN − primary\nreference similarity'); clean(ax)
    ax.text(.02,.02,'Matched pairs: 11 / 11 / 7',transform=ax.transAxes,fontsize=5.5,color='#475569')

f2func=[f2a,f2b,f2c,f2d,f2e,f2f]
for nm,fn in zip('ABCDEF',f2func):
    fig,ax=plt.subplots(figsize=(3.25,2.45)); fn(ax); letter(ax,nm); save_panel(fig,OUT/'Figure2'/f'Panel_{nm}')
fig=plt.figure(figsize=(7.35,7.65)); gs=GridSpec(3,2,figure=fig,hspace=.78,wspace=.66)
axes=[fig.add_subplot(gs[i,j]) for i in range(3) for j in range(2)]
for ax,nm,fn in zip(axes,'ABCDEF',f2func): fn(ax); letter(ax,nm)
save(fig,OUT/'Figure2/Figure2')

