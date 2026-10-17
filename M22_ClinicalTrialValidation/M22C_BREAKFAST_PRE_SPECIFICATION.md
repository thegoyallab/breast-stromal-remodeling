# Stage M22C — BREAKFAST paired-treatment stromal validation

Date: 16 October 2026
Status: FROZEN BEFORE BREAKFAST OUTCOME ANALYSIS
Trial: NCT04248998 / BREAKFAST
Resource: GSE279219 / PRJNA1171406
Parent: M22/M22B clinical-trial stromal validation

## Rationale and scope

BREAKFAST provides a second independent treatment-state cohort with public paired tumour single-cell data. The current public molecular crosswalk resolves 10 complete baseline (T0) / early-on-treatment (T1, approximately 14–20 days) patient pairs. Sample-level trial arm and pCR linkage are unresolved, so this cohort is used only for paired treatment-state validation; no response claim is permitted.

The analysis is outcome-independent and uses only already-frozen stromal endpoints from the Nature Communications manuscript.

## Input objects

- Primary-QC common-gene H5AD:
  AtlasExpansion_Readiness_20261014_v0_9/GSE279219_PRIMARY_QC_COMMON35371_UNINTEGRATED.h5ad
- Broad-lineage annotation:
  AtlasExpansion_Lineage_20261014_v0_10/GSE279219_broad_lineage_annotations.tsv.gz
- Frozen T0/T1 crosswalk:
  AtlasExpansion_BREAKFAST_20261015_v0_31/BREAKFAST_GSE279219_patient_timepoint_crosswalk_v0_1.tsv

Only cells assigned by the already-frozen provisional_joint_lineage as Fibroblast, Endothelial or Perivascular are used. Ambiguous cells are not reclassified.

## Eligibility gate established before outcomes

All 10 complete T0/T1 patients have >=20 cells at both timepoints in each of:
- Fibroblast
- Endothelial
- Perivascular

Therefore the >=20-cell threshold is retained without relaxation. >=10 remains sensitivity only.

## Endpoint family A — broad-lineage abundance

For each sample and target lineage:
log((n_lineage + 0.5) / (n_all_other non-ambiguous retained broad-lineage cells + 0.5))

Raw fractions are retained descriptively.
Three paired T0/T1 tests; BH correction across the three lineages.

## Endpoint family B — fixed within-lineage expression programs

Use exactly the M22 pre-specified gene sets:
- fibroblast matrix-contractile: COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, FN1, ACTA2, ITGA11, POSTN, JAG1
- mural/perivascular matrix-contractile: COL1A1, COL1A2, COL3A1, ACTA2, JAG1
- endothelial NOTCH-ECM: JAG2, COL4A1, COL4A2, NOTCH4, ITGA1

Patient x timepoint x lineage raw-count pseudobulk; log2(CPM+1); within-lineage gene standardization across evaluable pseudobulks; unweighted mean score; no fitted weights or post-hoc gene replacement.

Three paired T0/T1 tests; BH correction across lineages.

## Endpoint family C — frozen Reed fine-state balances

Mapping uses the exact historical Reed correlation geometry:
- donor-state Reed centroids;
- frozen fine_state_markers.tsv;
- log1p(count × 10,000 / total cell UMI);
- Pearson correlation;
- nearest state by maximum correlation.

Primary state balances:
- Fibroblast: log((FB1 + FB3 + 0.5)/(FB2 + FB4 + 0.5))
- Perivascular: log((PV3 + 0.5)/(PV1 + PV2 + PV4 + PV5 + 0.5))
- Endothelial: log((LE2 + 0.5)/(LE1 + 0.5))

Three paired T0/T1 tests; BH correction across lineages.

## Fine-state applicability gate

The fine-state endpoint for a lineage is interpretable only if the median cell-level max correlation across the BREAKFAST paired cohort is at least the 5th percentile of the original Wave-A core-primary distribution for that lineage, AND the median margin is at least the corresponding Wave-A 5th percentile.

Frozen Wave-A core-primary gates:
- Fibroblast: max_corr >= 0.3168938385; margin >= 0.00269698504
- Endothelial: max_corr >= 0.3484501570; margin >= 0.0020778835
- Perivascular: max_corr >= 0.3765755535; margin >= 0.003566771765

If a lineage fails this applicability gate, its fine-state P value is retained for audit but is not used as biological validation.

## Statistics

Patient is the biological replicate.
Paired T0/T1 comparisons use two-sided Wilcoxon signed-rank tests.
BH correction is performed separately for endpoint families A, B and C.
Report paired n, median change, raw medians/fractions, P and FDR.

## Robustness

Any promoted result must:
- retain direction at the >=10-cell sensitivity threshold;
- retain median direction under leave-one-patient-out analysis;
- not depend on redefining a lineage, state, program or endpoint.

## Interpretation boundary

BREAKFAST can support treatment-associated stromal remodeling but cannot currently support arm-specific or response-specific claims because sample-level trial arm and response linkage are unresolved.

No result from M22C changes the original normal-versus-tumour discovery claim.
