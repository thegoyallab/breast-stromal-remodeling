import pandas as pd, os
base='/mnt/stk35data/BreastCancer_Census_v1/StageM20B_MethodConcordance_v0_1/results'
jobs=['FB1plusFB3','PV3','VEAT','VEATplusVEC']
targets=['COL1A1','COL1A2','COL3A1','COL4A1','COL4A2','COL5A1','COL5A2','FN1','ACTA2','MYL9','TAGLN','ITGA11','POSTN','JAG1','JAG2','NOTCH3','NOTCH4','ITGA1','PECAM1','VWF','CDH5']
rows=[]
for j in jobs:
    p=os.path.join(base,f'M20B_{j}_edgeR_voom_gene_concordance.tsv')
    d=pd.read_csv(p,sep='\t')
    x=d[d.gene_symbol.isin(targets)].copy()
    x.insert(0,'contrast',j)
    x['claim_gene']=True
    rows.append(x)
out=pd.concat(rows,ignore_index=True)
out.to_csv(os.path.join(base,'M20B_manuscript_target_gene_concordance.tsv'),sep='\t',index=False)
summ=out.groupby('contrast').agg(
    n_target_genes=('gene_symbol','size'),
    direction_concordant=('direction_concordant','sum'),
    edgeR_FDR05=('edgeR_FDR',lambda s:(s<0.05).sum()),
    voom_FDR05=('voom_FDR',lambda s:(s<0.05).sum()),
    both_FDR05=('both_FDR05','sum')
)
print(summ.to_string())
print('\nPotential target-gene discordances:')
print(out[(~out.direction_concordant) | ((out.edgeR_FDR<0.05)!=(out.voom_FDR<0.05))][['contrast','gene_symbol','edgeR_log2FC','edgeR_FDR','voom_log2FC','voom_FDR']].to_string(index=False))
