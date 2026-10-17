# Stage M17A — vascular spatial, clinical and phosphoproteomic extension
**Date:** 13 October 2026
**Branch:** `StageM17A_VascularSpatialClinicalPhospho_v0_1`
**Status:** EXTENSION COMPLETE — spatial/prognostic gates negative; phosphoproteomic gate positive; IMC semantic correction required.

## 1. Source-annotation correction discovered during spatial analysis
The Danenberg single-cell annotation file defines:
- `is_perivascular`: **whether the cell is in contact with a vessel**, defined from vessel masks.
- `cellPhenotype = Endothelial`: deposited endothelial cell identity.

Therefore `is_perivascular` is a **spatial vessel-contact flag**, not a pericyte/perivascular/mural lineage label. It includes endothelial and multiple other phenotypes.

Historical analyses that called `is_perivascular==1` a "perivascular gate" are valid only if the term is understood as a **vessel-contacting spatial compartment**. They must not be presented as direct protein measurements of transcriptomically defined perivascular/mural cells.

## 2. Outcome-blind spatial hypothesis test
A lineage-independent mural-like proxy was locked before survival analysis:
- Myofibroblasts
- Myofibroblasts PDPN+

Primary endpoint:
- negative log2 enrichment of myofibroblast vessel-mask contact relative to the within-image expectation.

Independent secondary endpoint:
- negative log2 enrichment of myofibroblast–endothelial graph edges relative to the within-image random-mixing expectation.

Eligibility:
- >=20 mural-like cells for the primary endpoint;
- >=20 mural-like and >=20 endothelial cells for the graph endpoint.

Eligible patient-condition counts:
- primary metric: 25 normal, 12 DCIS, 631 invasive tumour;
- graph metric: 25 normal, 11 DCIS, 515 invasive tumour.

The two tumour-level metrics correlate moderately (Spearman rho=0.580, P=1.22e-47), showing that they capture related spatial structure. However, the tissue-context hypothesis did **not** replicate across metrics:
- primary normal vs tumour: delta decoupling +0.261, P=0.0550, FDR across six locked tests=0.304;
- graph normal vs tumour: delta -0.0246, P=0.875, FDR=0.875.

**Spatial promotion gate: FAIL.**
The large-scale Danenberg data do not provide robust evidence for a general tumour-associated mural-like/endothelial "detachment" phenotype under these prespecified definitions.

Do not claim pericyte detachment from these data.

## 3. METABRIC prognosis
The locked spatial scores were joined to the Danenberg clinical table only after the spatial phenotype was frozen.

Clinical overlap:
- 682 Danenberg invasive-tumour patients;
- discovery: 517 patients / 159 breast-cancer deaths;
- validation: 165 patients / 60 deaths.

Primary endpoint:
- breast-cancer-specific survival;
- continuous spatial-decoupling score scaled from the discovery cohort;
- Cox model adjusted for grade, tumour-size category, lymph-node burden and ER status.
- HER2 was a prespecified sensitivity covariate because of missing values.

Primary spatial score:
- discovery adjusted HR per 1 SD = 1.104 (95% CI 0.939–1.298), P=0.231;
- validation adjusted HR = 0.991 (95% CI 0.769–1.277), P=0.943.

Secondary graph score:
- discovery HR=0.952, P=0.576;
- validation HR=0.828, P=0.197.

Prespecified myofibroblast contractile score:
- discovery HR=1.250, P=0.146;
- validation HR=1.461, P=0.171.

No score violated the simple Schoenfeld-residual/time diagnostic at P<0.05.

**Clinical validation gate: FAIL.**
There is no supported prognostic claim for these spatial or contractile IMC scores.

## 4. CPTAC paired phosphoproteomics
A predefined module was tested in the 18 paired CPTAC prospective breast tumour–normal cases.

Modules:
1. ECM/receptor context: COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, ITGA1, ITGB1, PDGFRB.
2. Focal adhesion: PTK2, PXN, VCL, TLN1, SRC.
3. Rho/actomyosin: RHOA, ROCK1, ROCK2, MYL9, MYH9, ACTA2.
4. Cofilin axis: LIMK1, LIMK2, CFL1.
5. YAP/TAZ mechanotransduction: YAP1, WWTR1, TEAD3.

For each phosphosite, the raw tumour–normal change was retained, and where matched total-protein data existed a **protein-adjusted phosphosite delta** was calculated:
> phosphosite tumour–normal delta − corresponding total-protein tumour–normal delta.

This is a relative site-abundance measure, not phosphorylation stoichiometry.

Results:
- 45 module phosphosites had >=5 case-matched protein-adjusted pairs;
- 37/45 were FDR<0.05 within the predefined adjusted-site family;
- 12 genes contained >=1 significant protein-adjusted phosphosite;
- all five predefined module categories contained at least one significant adjusted site.

Representative results:
- PDGFRB S641: raw delta +2.242, global phosphosite FDR=0.00322; protein-adjusted delta +3.602, module FDR=0.00129.
- PTK2 Y570: raw +2.875, global FDR=6.96e-5; adjusted +2.894, module FDR=8.6e-5.
- PXN S533: raw +1.141, global FDR=0.00285; adjusted +1.673, module FDR=0.000286.
- VCL S721: raw +0.526, global FDR=0.00345; adjusted +2.026, module FDR=8.6e-5.
- TLN1 S446: raw +1.005, global FDR=0.0177; adjusted +1.903, module FDR=0.000549.
- ROCK1 S1341: raw +3.084, global FDR=0.0131; adjusted +3.407, module FDR=0.0121.
- MYH9 S1943: raw +1.884, global FDR=6.96e-5; adjusted +1.896, module FDR=8.6e-5.
- CFL1 S41: raw +1.108, global FDR=0.000546; adjusted +1.659, module FDR=0.000549.
- YAP1 S251: raw +2.376, global FDR=6.96e-5; adjusted +3.207, module FDR=8.6e-5.
- TEAD3 S148: raw +1.490, global FDR=0.000820; adjusted +1.243, module FDR=0.00214.

**Molecular gate: PASS.**
The data support coordinated phosphosite remodeling across PDGFR/focal-adhesion/actomyosin/YAP-associated machinery in paired breast tumours, including changes that exceed corresponding total-protein shifts.

Interpretation lock:
- Allowed: "coordinated phosphosite remodeling" or "protein-adjusted phosphosite remodeling".
- Not allowed without site-specific functional evidence: "pathway activation", "FAK activation", "ROCK activation", "YAP activation" or direct proof of collagen–integrin signal transmission.
- CPTAC is bulk tissue and cannot localize these phosphorylation events specifically to mural/perivascular cells.

## 5. Independent IMC mural-like marker reanalysis
Because `is_perivascular` is a vessel-contact flag rather than a cell lineage, the core markers were rerun using only deposited myofibroblast phenotypes.

Eligible >=20-cell patients:
- normal 25;
- DCIS 12;
- invasive tumour 632.

### α-SMA
- invasive vs normal: median difference +0.681, P=1.12e-6, FDR15=1.68e-5.
- DCIS vs normal: +1.034, P=6.98e-4, FDR15=0.00523.
- DCIS vs invasive: +0.353, P=0.110, FDR15=0.207.

Thus the **contractile α-SMA phenotype is independently supported in a myofibroblast-defined compartment** and is detectable in DCIS and invasive cancer.

### Caveolin-1
- invasive vs normal: difference -0.099, P=0.0715, FDR15=0.153.
- DCIS vs normal: -0.068, P=0.339, FDR15=0.508.

Therefore the strong historical caveolin-1 decrease pertains to the **vessel-contacting spatial compartment**, not demonstrably to the myofibroblast-defined mural-like population.

## 6. Manuscript implications
### Required correction
Replace any statement equating the Danenberg `is_perivascular` field with a perivascular/mural cell lineage.

Recommended terminology:
> "vessel-contacting compartment (the Danenberg `is_perivascular` spatial annotation)"

The myofibroblast reanalysis can then provide the lineage-independent contractile validation:
> "An independent analysis restricted to deposited myofibroblast phenotypes confirmed increased α-SMA in invasive tumours and DCIS relative to normal tissue."

Caveolin-1 should be described separately:
> "Caveolin-1 was reduced in the broader vessel-contacting compartment, but this decrease was not significant when analysis was restricted to myofibroblast phenotypes."

### Do not add
- pericyte detachment as a positive result;
- METABRIC prognostic claims for the spatial/contractile scores.

### Optional high-value addition
Add the CPTAC phosphoproteomic module as a supplementary mechanistic layer. It directly addresses the manuscript limitation that ligand–receptor analysis measures expression potential rather than downstream signaling state, while remaining explicitly bulk-tissue and non-causal.

## Overall M17A decision
**Spatial detachment: FAIL.**
**Prognostic relevance: FAIL.**
**CPTAC phosphoproteomic remodeling: PASS.**
**IMC α-SMA contractile validation after semantic correction: PASS.**
**IMC caveolin-1 as mural-cell-intrinsic validation: FAIL; retain only as vessel-contact compartment evidence.**

The negative spatial and survival results should be retained in the audit trail but should not expand the manuscript. The phosphoproteomic result and the corrected α-SMA validation are suitable for controlled propagation into the next manuscript freeze.
