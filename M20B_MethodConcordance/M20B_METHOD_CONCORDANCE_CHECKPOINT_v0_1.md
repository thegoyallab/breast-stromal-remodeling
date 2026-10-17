# Stage M20B — edgeR/limma-voom method-concordance validation

Date: 13 October 2026
Parent scientific freeze: M18
Purpose: methodological validation only; no new biological hypothesis.

## Design
Primary method remains edgeR quasi-likelihood with TMM normalization and the exact frozen patient-level pseudobulk count matrices.
Sensitivity method: limma-voom, using the identical filtered count matrices, sample inclusion criteria and Normal-versus-Primary design.
Principal contrasts tested:
- scANVI FB1+FB3 aggregate fibroblast axis, floor >=20 cells/patient
- scANVI PV3, floor >=20
- scANVI VEAT, floor >=20
- scANVI VEAT+VEC, floor >=20

## Global concordance results
FB1+FB3: Spearman log2FC rho=0.9926; Pearson r=0.9889; union-FDR-significant direction concordance=99.97%; FDR<0.05 Jaccard=0.934.
PV3: Spearman rho=0.9932; Pearson r=0.9941; union-FDR-significant direction concordance=100%; Jaccard=0.925.
VEAT: Spearman rho=0.9892; Pearson r=0.9715; union-FDR-significant direction concordance=100%; Jaccard=0.807.
VEAT+VEC: Spearman rho=0.9933; Pearson r=0.9907; union-FDR-significant direction concordance=100%; Jaccard=0.906.

## Manuscript-highlighted genes
For the prespecified manuscript-relevant gene set, all biologically promoted genes retained the same direction between edgeR and voom.

Observed method differences were limited to:
- FB1+FB3 COL4A1: nonsignificant by both methods and sign-unstable; not a promoted fibroblast claim.
- VEAT+VEC TAGLN: nonsignificant by both methods and sign-unstable; not a promoted endothelial claim.
- VEAT+VEC ACTA2: positive by both methods; edgeR FDR=0.0472 versus voom FDR=0.0991. This is a threshold-level significance difference and should not be promoted as a method-robust endothelial claim.

All principal promoted ECM/contractile fibroblast and PV3 genes, and the core endothelial JAG2/NOTCH4/COL4A1/COL4A2/ITGA1/PECAM1/VWF/CDH5 program, are directionally concordant. The exact per-gene table is in results/M20B_manuscript_target_gene_concordance.tsv.

## Gate decision
PASS for edgeR-versus-limma-voom concordance.

Interpretation:
- edgeR remains the primary method.
- limma-voom independently supports the principal effect directions and nearly all significance calls.
- No central manuscript conclusion depends on an edgeR-only sign reversal.
- The VEAT+VEC ACTA2 threshold difference is explicitly classified as method-sensitive and is not needed for the core endothelial conclusion.

## DESeq2
DESeq2 is not present in the existing isolated R environment. It is not required to establish the primary two-method concordance gate. If desired, DESeq2 can be added later in a separate isolated environment; it should not alter or replace the frozen edgeR primary analysis.

## Scientific-freeze rule
This stage validates the analytical framework. It does not reopen biological discovery and does not modify M18/M19 unless a concise methods/sensitivity statement is consciously added at final editorial integration.
