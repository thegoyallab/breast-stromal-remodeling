# Stage M22 — Pre-specified clinical-trial stromal validation

Date: 16 October 2026
Status: PRE-SPECIFIED BEFORE EXPRESSION-OUTCOME ANALYSIS
Parent manuscript: AtlasExpansion_PreSubmissionIntegrity_20261016_v0_97
Primary validation cohort: FELINE / NCT02712723 / GSE211434

## Purpose
Test whether stromal programs defined before inspection of this clinical-trial expression dataset show longitudinal treatment-associated behavior in an independent, serially sampled breast-cancer cohort. This is a bounded validation of the existing Nature Communications stromal model, not a new atlas analysis.

## Cohort
GSE211434 is the independent FELINE validation cohort.
- 27 ER-positive breast-cancer patients.
- all 27 have day 0, day 14 and day 180 serial biopsies in the source metadata.
- 16 received ribociclib + letrozole; 11 received letrozole alone.
- 15 are source-defined sensitive; 12 resistant.
- source patient, timepoint, arm and response labels will be used without reclassification.

## Pre-specified cellular compartments
Use source cell-type annotations.
1. Fibroblasts.
2. Pericytes/mural cells.
3. Endothelial cells.

No immune-only cohort will substitute for a missing stromal compartment. No new cell type will be promoted after outcome inspection.

## Pre-specified expression programs
These gene sets are derived from the already-frozen manuscript findings.

Fibroblast matrix-contractile program:
COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, FN1, ACTA2, ITGA11, POSTN, JAG1

Mural/pericyte matrix-contractile program:
COL1A1, COL1A2, COL3A1, ACTA2, JAG1

Endothelial NOTCH-ECM program:
JAG2, COL4A1, COL4A2, NOTCH4, ITGA1

Secondary endothelial identity/context markers, reported separately and not used to redefine the primary program:
PECAM1, VWF, CDH5

## Pre-specified abundance endpoints
For every biopsy with sufficient cells:
- fibroblast fraction of all retained source-annotated cells;
- pericyte fraction;
- endothelial fraction.

## Pre-specified expression endpoint construction
Within each patient x day x lineage:
- aggregate raw counts to patient-lineage pseudobulk;
- retain exact cell counts and library sizes;
- require >=20 cells for the primary endpoint;
- repeat at >=10 cells as a sparse-lineage sensitivity only;
- calculate normalized expression from pseudobulk, not cell-level pseudo-replication;
- calculate each fixed program score from the available member genes without changing membership based on observed direction.

If a named gene is absent from the deposited matrix, record it as unavailable; do not replace it post hoc.

## Pre-specified contrasts
Primary:
1. paired day 0 -> day 14 program change;
2. paired day 0 -> day 180 program change;
3. sensitive versus resistant difference in those within-patient changes.

Secondary:
4. letrozole-alone versus ribociclib+letrozole difference in within-patient change;
5. broad-lineage abundance change over the same intervals.

The day-14 contrast tests early pharmacodynamic remodeling; day-180 tests persistent/end-of-treatment remodeling.

## Statistical unit and tests
Patient is the biological replicate.
- paired longitudinal contrasts: two-sided paired tests or patient-blocked pseudobulk models as appropriate to the endpoint;
- response/arm comparisons: compare patient-level change values, not cells;
- report effect sizes, evaluable patient counts and exact/nominal P values;
- Benjamini-Hochberg correction within the pre-specified stromal-program family;
- no significance claim from cell counts alone.

A patient lacking an evaluable lineage at either member of a paired contrast is excluded only from that lineage/contrast and remains eligible elsewhere.

## Interpretation gates
A program will be described as longitudinally supported only if:
- direction is defined from the pre-specified score;
- the paired patient-level effect is FDR-supported or shows a consistent effect with an explicitly reported non-significant result;
- result is not driven solely by one patient on leave-one-patient-out sensitivity;
- minimum-cell-threshold sensitivity does not reverse the principal direction.

Response association requires a difference in within-patient change between source-defined sensitive and resistant groups; cross-sectional responder/nonresponder differences are not sufficient.

## Manuscript-integration rule
This analysis will not change the existing discovery claims.
- Positive result: add as independent treatment-aware validation, preferably Supplementary with one concise main-text sentence if compelling.
- Null result: report as a boundary/negative validation if informative.
- Discordant result: report as context dependence; do not redefine the original program.
- No post hoc gene-set optimization is permitted.

## Data integrity
Original GEO files are retained read-only with SHA256 hashes. Derived objects are written to this M22 directory. The v0.97 submission package is not modified during analysis.
