#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1')
OUT=ROOT/'StageM15F_FinalScientificFreeze_v0_1/results/lr_corrected/stage09'
RES=OUT
DE_RES=ROOT/'Stage09_StromalPrograms_v0_1/results'
cts=['Epithelial','Fibroblast','Endothelial','Perivascular']

# DE lookup: best FDR row per symbol
DE={}
for ct in cts:
    x=pd.read_csv(DE_RES/f'DE_GSE161529_{ct}_Primary_vs_Normal_edgeR_v0_1.tsv',sep='\t')
    x=x.sort_values(['gene_symbol','FDR','PValue']).drop_duplicates('gene_symbol')
    DE[ct]=x.set_index('gene_symbol')[['logFC','FDR']].to_dict('index')

m=pd.read_csv(RES/'LR_integrated_LIANA_pseudobulk_primary_vs_normal_v0_1.tsv.gz',sep='\t')
rows=[]
for i,r in m.iterrows():
    sign=np.sign(r.delta_median_log2score)
    lig=r.ligand_complex.split('_')
    rec=r.receptor_complex.split('_')
    details=[]
    support=0; opposite=0; total_sig=0
    for side,ct,genes in [('L',r.source,lig),('R',r.target,rec)]:
        for g in genes:
            d=DE.get(ct,{}).get(g)
            if d is None:
                details.append(f'{side}:{g}:NA')
                continue
            fc=float(d['logFC']); fdr=float(d['FDR'])
            details.append(f'{side}:{g}:{fc:.3f}:{fdr:.3g}')
            if fdr<0.05:
                total_sig+=1
                if sign!=0 and np.sign(fc)==sign: support+=1
                elif sign!=0 and np.sign(fc)==-sign: opposite+=1
    rows.append((support,opposite,total_sig,';'.join(details)))
m[['de_support_same_direction','de_opposite_significant','de_significant_subunits','de_details']]=pd.DataFrame(rows,index=m.index)

# stringent mechanistic gate
m['mechanistic_gate']=(
    (m['FDR_global']<1e-5) &
    (m['FDR_route']<1e-4) &
    (m['delta_median_log2score'].abs()>=0.75) &
    (m['magnitude_rank_primary'].fillna(1)<=0.10) &
    (m['primary_cohorts_evaluable']>=2) &
    (m['primary_cohort_direction_consistency']>=0.75) &
    (m['direction_concordant']==True) &
    (m['de_support_same_direction']>=1) &
    (m['de_opposite_significant']==0)
)
m['priority_score']=(
    m['delta_median_log2score'].abs() *
    (-np.log10(m['FDR_global'].clip(lower=1e-300))) *
    (1-m['magnitude_rank_primary'].fillna(1))
)
m=m.sort_values(['mechanistic_gate','priority_score'],ascending=[False,False])
m.to_csv(RES/'LR_integrated_prioritized_v0_1.tsv.gz',sep='\t',index=False,compression='gzip')
short=m[m.mechanistic_gate].copy()
short.to_csv(RES/'LR_mechanistic_shortlist_v0_1.tsv',sep='\t',index=False)

# top 100 for figure/manual audit, balanced by direction
up=short[short.delta_median_log2score>0].head(60)
dn=short[short.delta_median_log2score<0].head(40)
top=pd.concat([up,dn]).sort_values('priority_score',ascending=False)
top.to_csv(RES/'LR_top100_manual_audit_v0_1.tsv',sep='\t',index=False)

summary=(short.groupby(['source','target'],observed=True)
         .agg(candidates=('mechanistic_gate','size'),
              increased=('delta_median_log2score',lambda x:int((x>0).sum())),
              decreased=('delta_median_log2score',lambda x:int((x<0).sum())),
              median_abs_effect=('delta_median_log2score',lambda x:float(np.median(np.abs(x)))))
         .reset_index())
summary.to_csv(RES/'LR_mechanistic_shortlist_route_summary_v0_1.tsv',sep='\t',index=False)
print('MECHANISTIC SHORTLIST',len(short))
print(summary.to_string(index=False))
print('\nTOP 40')
cols=['source','target','ligand_complex','receptor_complex','delta_median_log2score',
      'FDR_global','magnitude_rank_primary','primary_cohort_direction_consistency',
      'de_details','priority_score']
print(short[cols].head(40).to_string(index=False))
