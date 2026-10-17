import pandas as pd, os
base='/mnt/stk35data/BreastCancer_Census_v1/StageM20B_MethodConcordance_v0_1/results'
jobs=['FB1plusFB3','PV3','VEAT','VEATplusVEC']
frames=[]
for j in jobs:
    d=pd.read_csv(os.path.join(base,f'M20B_{j}_edgeR_voom_gene_concordance.tsv'),sep='\t')
    d.insert(0,'contrast',j)
    frames.append(d)
out=pd.concat(frames,ignore_index=True)
out.to_csv('/mnt/stk35data/BreastCancer_Census_v1/StageM20B_MethodConcordance_v0_1/candidate/TableS19_method_concordance_full.tsv',sep='\t',index=False)
print(out.shape)
