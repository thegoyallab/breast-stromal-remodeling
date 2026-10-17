# Stage M22 — Clinical-trial stromal validation checkpoint

Date: 16 October 2026
Parent manuscript package: AtlasExpansion_PreSubmissionIntegrity_20261016_v0_97
Status: ANALYSIS COMPLETE — MIXED/BOUNDARY RESULT; MANUSCRIPT PACKAGE NOT MODIFIED

## Objective

Test whether the stromal axes already defined in the Nature Communications manuscript behave dynamically during treatment in independent longitudinal clinical-trial cohorts, without reopening discovery or optimizing signatures after outcome inspection.

Two cohorts passed the technical stromal-eligibility screen before outcome analysis:
1. FELINE / NCT02712723 / GSE211434.
2. BREAKFAST / NCT04248998 / GSE279219.

Other treatment resources were held or excluded for this bounded analysis because of immune enrichment, absence of a directly deposited/frozen stromal annotation layer, controlled-access patient linkage, or unresolved response/sample linkage. The full decision ledger is M22_TRIAL_STROMAL_ELIGIBILITY.tsv.

## FELINE: response-linked three-timepoint validation

GSE211434 contains 251,421 deposited cells from 27 ER-positive breast-cancer patients with D0, D14 and D180 biopsies. Source labels define 15 sensitive and 12 resistant patients; 16 received ribociclib + letrozole and 11 letrozole alone.

Stromal representation was strong at the prespecified >=20-cell threshold:
- Fibroblasts: 25 paired D0/D14 and 25 paired D0/D180 patients.
- Pericytes: 20 paired D0/D14 and 18 paired D0/D180.
- Endothelial cells: 24 paired D0/D14 and 24 paired D0/D180.

### Broad-lineage abundance

A pre-specified one-versus-rest log-ratio analysis identified one FDR-supported longitudinal change:
- endothelial D0 -> D180: n=24, median raw fraction 5.28% -> 7.72%, median log-ratio change +0.492, P=0.00178, FDR=0.0107.

This direction was retained:
- at the >=10-cell sensitivity threshold (n=25, P=0.000912, FDR=0.00547);
- in every leave-one-patient-out median-direction analysis.

Fibroblasts showed an early D14 increase at nominal P=0.0342 but not after the six-test family correction (FDR=0.103). Pericyte abundance changes were not FDR-supported.

No broad-lineage change differed significantly between source-defined sensitive and resistant patients after FDR correction, and no treatment-arm difference in longitudinal change was FDR-supported.

### Fixed transcriptional programs

All 15 pre-specified genes were present and the count matrix matched the annotation CellID universe exactly (251,421/251,421).

None of the fixed fibroblast matrix-contractile, mural matrix-contractile or endothelial NOTCH-ECM program scores changed significantly at D14 or D180 after FDR correction. No response-associated or treatment-arm-associated program change was FDR-supported.

Therefore FELINE does not validate these fixed expression programs as generic pharmacodynamic markers.

### Reed fine-state transfer

A frozen-marker Reed mapping pilot used all 1,322 historical mapper genes and did not use treatment or response labels. However, mapping similarity was far outside the original Wave-A primary-tumour domain:
- Fibroblast median max correlation: -0.0043 versus Wave-A core-primary 5th-percentile benchmark 0.3169.
- Endothelial: 0.1542 versus 0.3485.
- Perivascular: 0.1420 versus 0.3766.

The deposited nCount_RNA field also showed a small approximately 3-4% discrepancy from direct count-table column sums in a 10-cell spot audit. This is not large enough to explain the correlation gap, but it reinforces the decision not to rescue or reinterpret the mapping post hoc.

Accordingly, apparent FELINE fine-state P values are retained for audit only and are NOT accepted as biological validation.

## BREAKFAST: independent early-treatment paired validation

GSE279219 contributes 10 complete public T0/T1 pairs, where T1 is approximately 14-20 days after treatment initiation. Sample-level trial arm and response linkage remain unresolved, so only paired treatment-state testing was permitted.

The already-frozen expansion object contains 253,150 primary-QC cells. All 10 paired patients had >=20 cells at both timepoints in each target lineage. The analysis used 63,769 non-ambiguous target stromal cells:
- 47,532 fibroblasts;
- 10,406 endothelial cells;
- 5,831 perivascular cells.

### Reference applicability

Unlike FELINE, BREAKFAST passed the predeclared Wave-A-domain applicability gate for every lineage:
- Fibroblast median max correlation 0.546, median margin 0.0297.
- Endothelial 0.523, margin 0.0257.
- Perivascular 0.528, margin 0.0226.

All values exceeded the corresponding Wave-A core-primary 5th-percentile gates.

### Paired treatment-state results

No endpoint family showed an FDR-supported T0 -> T1 change.

Broad abundance:
- Fibroblast: median balance change +0.439, P=0.131, FDR=0.196.
- Endothelial: -0.173, P=0.375, FDR=0.375.
- Perivascular: -0.641, P=0.105, FDR=0.196.

Fixed programs:
- Fibroblast matrix-contractile: P=0.193, FDR=0.349.
- Mural matrix-contractile: P=0.770, FDR=0.770.
- Endothelial NOTCH-ECM: P=0.232, FDR=0.349.

Frozen Reed state balances:
- FB1+FB3 versus FB2+FB4: P=0.695, FDR=1.0.
- PV3 versus rest: P=0.492, FDR=1.0.
- LE2 versus LE1: P=1.0, FDR=1.0.

Because the mapper passed the applicability gate, these are interpretable negative early-treatment results rather than mapping failures.

## Scientific conclusion

The clinical-trial layer does NOT support a claim that the cross-sectional tumour-associated stromal architecture is a universal treatment-response biomarker.

Instead, it establishes a narrower and defensible boundary:
- early treatment in BREAKFAST TNBC did not measurably shift the pre-specified stromal abundance, expression-program or Reed-state endpoints despite good reference applicability;
- prolonged FELINE endocrine/CDK4/6 treatment was associated with a robust increase in broad endothelial representation, but not with an FDR-supported fixed transcriptional program or response-associated stromal change;
- FELINE fine-state mapping is out of the validated reference domain and is not interpreted.

This argues that tumour-associated stromal remodeling is treatment-context- and timescale-dependent and should not be presented as uniformly reversible with therapy.

## Manuscript recommendation

Use M22 as a bounded longitudinal stress test, not as a new headline biological result.

Recommended integration:
1. one concise Results paragraph;
2. one Supplementary table summarizing cohort eligibility, applicability and paired effects;
3. optionally one compact Supplementary figure showing the FELINE endothelial paired change plus a forest/summary view of all pre-specified trial endpoints;
4. one Discussion sentence stating that two longitudinal cohorts do not support a generic pharmacodynamic interpretation of the cross-sectional stromal axes.

Do NOT:
- change the title;
- add clinical-trial claims to the abstract unless the editorial strategy specifically requires it;
- describe FELINE fine-state results as validated;
- claim response prediction;
- optimize gene sets or state endpoints after seeing these results.

The v0.97 package remains unchanged pending integration review.
