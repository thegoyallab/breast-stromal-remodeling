# Stage M22B — Prospective Reed-state longitudinal validation

Date: 16 October 2026
Status: FROZEN BEFORE FELINE REED-STATE MAPPING
Parent: M22 pre-specified GSE211434 validation
Reference framework: frozen Reed normal-breast stromal nearest-state correlation mapper used in the Nature Communications manuscript

## Objective

Test whether treatment changes the pre-existing tumour-associated stromal state architecture in serial FELINE biopsies. No FELINE-derived clustering, marker selection or state discovery is permitted.

## Mapping method

Use the exact frozen historical mapper geometry:
- Reed donor-by-state pseudobulk reference centroids;
- frozen fine_state_markers.tsv marker universe;
- per-cell log1p(count × 10,000 / total UMI);
- Pearson correlation across lineage-specific marker space;
- nearest state = highest reference-state correlation;
- mapping confidence = maximum correlation and best-minus-second-best margin.

Source lineage crosswalk:
- FELINE Fibroblasts -> Reed Fibroblast
- FELINE Pericytes -> Reed Perivascular
- FELINE Endothelial cells -> Reed Endothelial

No source response, treatment arm or timepoint is used during mapping.

## Primary state endpoints

1. Fibroblast remodeling balance
   log((FB1 + FB3 + 0.5) / (FB2 + FB4 + 0.5))
   Rationale: the current manuscript retains FB1+FB3 versus FB2+FB4 as the mapper-robust fibroblast axis.

2. Mural remodeling balance
   log((PV3 + 0.5) / (PV1 + PV2 + PV4 + PV5 + 0.5))
   Rationale: PV3 gain is mapper-robust; PV1 and PV2 losses are also robust. PV4 and PV5 remain in the complement so that the endpoint is a complete one-versus-rest balance and does not promote those states individually.

3. Lymphatic endothelial remodeling balance
   log((LE2 + 0.5) / (LE1 + 0.5))
   Rationale: LE1 loss and LE2 gain are the fine-state endothelial effects retained by both mapping frameworks. VEA and VEV/VEAT/VEC partitions are mapper-sensitive and are not primary treatment endpoints.

Raw mapped-state fractions are reported for interpretation but are not treated as independent primary outcomes.

## Eligibility

Primary: source lineage has >=20 cells at both timepoints in the pair.
Sensitivity: >=10 cells at both timepoints.

The same source-lineage eligibility determined before mapping is retained. No minimum mapped-state count is added after inspecting results.

## Primary contrasts

For each of the three fixed state balances:
- D0 -> D14 paired change;
- D0 -> D180 paired change.

Response association:
- compare patient-level D14-D0 change between source-defined Response and Non-response;
- compare patient-level D180-D0 change between Response and Non-response.

Patient is the biological replicate.

## Secondary contrasts

Treatment-arm differences in the same patient-level changes:
ribociclib + letrozole versus letrozole alone.

Descriptive supporting quantities:
- FB1+FB3 fraction;
- FB2+FB4 fraction;
- PV1, PV2 and PV3 fractions;
- LE1 and LE2 fractions;
- median max correlation and median margin by lineage/timepoint.

These descriptive quantities cannot replace the pre-specified balance tests.

## Statistics

- paired D0/timepoint comparisons: two-sided Wilcoxon signed-rank test;
- response and treatment-arm differences in within-patient change: two-sided Mann-Whitney U;
- Benjamini-Hochberg correction separately across six paired tests, six response-change tests and six treatment-arm-change tests;
- report exact evaluable patient n and effect sizes;
- >=10-cell analyses are sensitivity only.

## Robustness gates

Any promoted result must:
- retain the same direction at the >=10-cell sensitivity threshold;
- retain its median direction in every leave-one-patient-out reanalysis;
- not depend on changing the fixed state endpoint;
- be interpreted alongside mapping-confidence distributions.

## Interpretation

Positive state-balance result: longitudinal support that therapy modifies an already-defined tumour-associated stromal architecture.
Null result: evidence that the cross-sectional tumour-associated state architecture is not detectably reversed or intensified by this FELINE treatment context.
Response-associated result: only if within-patient balance change differs between source-defined sensitive and resistant groups after family-wise FDR control.

No result in M22B may redefine the original tumour-versus-normal discovery claims.
