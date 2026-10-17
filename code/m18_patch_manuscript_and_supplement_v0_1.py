#!/home/pankaj/venvs/bc_scanpy/bin/python
from pathlib import Path
import shutil,json

ROOT=Path("/mnt/stk35data/BreastCancer_Census_v1")
ST=ROOT/"StageM18_IntegratedRevision_v0_1"
MAN0=ST/"manuscript/NatureCommunications_Manuscript_M15F_WORKING.md"
MAN=ST/"manuscript/NatureCommunications_Manuscript_M18_WORKING.md"
SUP0=ST/"supplementary/Supplementary_Information_M15F_WORKING.md"
SUP=ST/"supplementary/Supplementary_Information_M18_WORKING.md"
shutil.copy2(MAN0,MAN); shutil.copy2(SUP0,SUP)

lines=MAN.read_text().splitlines()
def replace_prefix(prefix,new):
    n=0
    for i,x in enumerate(lines):
        if x.startswith(prefix):
            lines[i]=new; n+=1
    if n!=1: raise RuntimeError((prefix,n))
def replace_heading(old,new):
    n=0
    for i,x in enumerate(lines):
        if x==old: lines[i]=new;n+=1
    if n!=1: raise RuntimeError((old,n))

for i,x in enumerate(lines):
    if x.startswith("Tumour-associated single-cell differences conflate"):
        lines[i]=x.replace("Imaging mass cytometry, paired proteomics and a negative pre-specified interaction test provide orthogonal checks.",
                           "Imaging mass cytometry, paired proteomics, paired phosphoproteomics and a negative pre-specified interaction test provide orthogonal checks.")
        break
else: raise RuntimeError("abstract not found")

replace_prefix(
"Independent protein imaging supported the contractile mural phenotype",
"""Independent protein imaging supported a contractile mural-like phenotype and showed that it was already detectable in DCIS. In the Danenberg source annotation, is_perivascular denotes spatial contact with a vessel mask rather than a mural-cell lineage. We therefore repeated the lineage-oriented analysis using only cells with the deposited phenotypes Myofibroblasts or Myofibroblasts PDPN+. Within this myofibroblast-defined compartment, patient-level median α-smooth-muscle actin (α-SMA) was higher in invasive tumours than in normal breast (median difference=0.681, FDR=1.68×10⁻⁵) and was also higher in DCIS (difference=1.03, FDR=0.00523) (Fig. 5c). The direct cross-sectional DCIS-versus-invasive contrast was not FDR-supported (difference=0.353, FDR=0.207). Caveolin-1 was reduced in the broader vessel-contacting compartment used in the original spatial annotation, but this decrease was not significant after restriction to myofibroblast phenotypes (invasive versus normal FDR=0.153; DCIS versus normal FDR=0.508; Supplementary Table 9). Thus, the imaging data independently support a contractile α-SMA-high stromal phenotype in DCIS and invasive disease, whereas caveolin-1 loss is retained only as evidence from the broader vessel-contacting spatial compartment rather than as a mural-cell-intrinsic change."""
)

replace_prefix(
"Paired mass-spectrometry proteomics provided independent support for the matrix-remodeling fibroblast program.",
"""Paired mass-spectrometry proteomics provided independent support for the matrix-remodeling fibroblast program. Relative to matched normal tissue, tumour protein abundance increased for FAP (median paired difference=1.01, FDR=7.29×10⁻⁴), COL11A1 (1.86, FDR=3.53×10⁻⁴) and ITGA11 (0.804, FDR=0.00141) (Fig. 5d). The concordance between transcript-level and paired protein-level changes supports the biological relevance of this stromal matrix program while preserving the distinction between total-protein abundance and site-specific phosphorylation.

In the same prospective CPTAC resource, we next tested a predefined phosphoproteomic module spanning ECM/receptor context, focal-adhesion, Rho/actomyosin, cofilin and YAP/TAZ-associated proteins. For each phosphosite with matched total-protein measurements, the paired tumour-minus-normal phosphosite change was additionally expressed relative to the corresponding paired total-protein change for the same gene. Forty-five sites had at least five such matched pairs, and 37 were FDR<0.05 within this predefined protein-adjusted site family, distributed across 12 genes (Supplementary Fig. 13; Supplementary Table 18). Representative protein-adjusted changes included PDGFRB S641 (Δ=3.60, FDR=0.00129), PTK2 Y570 (Δ=2.89, FDR=8.6×10⁻⁵), PXN S533 (Δ=1.67, FDR=2.86×10⁻⁴), MYH9 S1943 (Δ=1.90, FDR=8.6×10⁻⁵), CFL1 S41 (Δ=1.66, FDR=5.49×10⁻⁴) and YAP1 S251 (Δ=3.21, FDR=8.6×10⁻⁵). Because these measurements are bulk-tissue and site-specific functional consequences differ by residue, we interpret them as coordinated phosphosite remodeling rather than direct evidence of pathway activation or stromal-cell-specific signalling flux."""
)

replace_prefix(
"The mural compartment showed an equally strong but less commonly emphasized pattern.",
"""The mural compartment showed an equally strong but less commonly emphasized pattern. Pericyte perturbation can influence hypoxia, epithelial–mesenchymal transition and metastatic behaviour in breast-cancer models [29], and recent human work implicates breast pericytes in tumour-cell proliferation and altered tumour localization [10]. We observed a robust loss of PV1/PV2-like affinities together with expansion of a contractile PV3-like program in invasive tumours, accompanied by increased ACTA2 and fibrillar-collagen expression within PV3-like cells. Imaging mass cytometry independently confirmed increased α-SMA in a myofibroblast-defined compartment in invasive tumours and DCIS. The Danenberg is_perivascular annotation was found to denote vessel contact rather than mural lineage identity; accordingly, the strong caveolin-1 decrease observed in that spatially defined compartment is no longer interpreted as mural-cell-intrinsic because it was not significant after restriction to myofibroblast phenotypes. The direct DCIS-versus-invasive α-SMA comparison was not FDR-supported, so these imaging data establish cross-sectional presence rather than temporal ordering. PV5-like abundance was more sensitive to reference-confidence filtering and is therefore retained as a secondary rather than defining result."""
)

for i,x in enumerate(lines):
    if x.startswith("Several limitations remain."):
        lines[i]=x.replace(
            "Finally, ligand–receptor inference remains transcript based and inherits uncertainty from fine-state assignment, while the available orthogonal protein datasets measure only a subset of the proposed mechanisms.",
            "Finally, ligand–receptor inference remains transcript based and inherits uncertainty from fine-state assignment. The paired phosphoproteomic analysis provides bulk-tissue evidence of site-level remodeling across relevant adhesion and mechanotransduction proteins, but it cannot localize those phosphorylation events to the stromal populations that generated the single-cell interaction hypotheses or by itself establish pathway activation."
        )
        break
else: raise RuntimeError("limitations paragraph not found")

replace_prefix(
"Protein-level stromal validation used the publicly available METABRIC breast-cancer imaging mass cytometry resources",
"""Protein-level stromal validation used the publicly available METABRIC breast-cancer imaging mass cytometry resources [6,7]. Source-data semantics were taken directly from the Danenberg annotations. The field is_perivascular denotes whether a cell is in contact with a vessel mask and was therefore treated as a vessel-contacting spatial annotation, not as a mural/pericyte lineage label. For lineage-oriented contractile validation, we instead defined a mural-like protein-imaging compartment using only the deposited phenotypes Myofibroblasts and Myofibroblasts PDPN+. Marker measurements were summarized by patient and condition after requiring at least 20 mural-like cells per patient-condition group. In parallel, the original vessel-contacting annotation was retained as a spatial-compartment sensitivity analysis, including epithelial exclusion and a stricter exclusion of epithelial plus CK5/pan-cytokeratin-high cells."""
)
replace_prefix(
"For each gating definition, invasive-tumour-versus-normal",
"""For both the mural-like phenotype analysis and vessel-contacting sensitivity analyses, invasive-tumour-versus-normal, DCIS-versus-normal and DCIS-versus-invasive patient distributions were compared using two-sided Mann–Whitney U tests. Five markers were tested for each of the three contrasts, and P values were adjusted by the Benjamini–Hochberg method across the resulting 15 marker-by-contrast tests within each definition. The direct DCIS-versus-invasive comparison was treated as cross-sectional and was not used to infer temporal progression. α-SMA was interpreted as evidence of a contractile protein phenotype rather than as a marker specific for an individual nearest-normal mural state; caveolin-1 changes from the source is_perivascular analysis were interpreted only at the level of the vessel-contacting spatial compartment."""
)

replace_heading("## Paired proteomic validation","## Paired proteomic and phosphoproteomic validation")
replace_prefix(
"Mass-spectrometry proteomic validation used the CPTAC/PDC prospective breast-cancer proteome",
"""Mass-spectrometry proteomic validation used the CPTAC/PDC prospective breast-cancer proteome and phosphoproteome [19]. Protein abundance measurements were collapsed to patient level, and patients with both primary-tumour and normal tissue measurements were analysed as paired observations. Eighteen matched patients contributed to the FAP, COL11A1 and ITGA11 comparisons shown in Fig. 5d. Tumour-minus-normal abundance differences were tested using two-sided Wilcoxon signed-rank tests, and P values were adjusted by the Benjamini–Hochberg method across quantified proteins.

A predefined phosphoproteomic module was specified before examining module-level significance and comprised ECM/receptor-context genes (COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, ITGA1, ITGB1 and PDGFRB), focal-adhesion genes (PTK2, PXN, VCL, TLN1 and SRC), Rho/actomyosin genes (RHOA, ROCK1, ROCK2, MYL9, MYH9 and ACTA2), a cofilin-axis set (LIMK1, LIMK2 and CFL1) and YAP/TAZ-associated genes (YAP1, WWTR1 and TEAD3). Phosphosite tumour-minus-normal differences were tested in paired cases using two-sided Wilcoxon signed-rank tests. Where the corresponding total-protein measurement was available in the same patient, a protein-adjusted site change was calculated as the phosphosite tumour-minus-normal difference minus the same-gene total-protein tumour-minus-normal difference. Sites with at least five paired protein-adjusted observations were included in the predefined adjusted-site family, and these P values were corrected by the Benjamini–Hochberg method. The original experiment-wide phosphosite FDR values calculated across all quantified CPTAC phosphosite rows were retained alongside the module-family results. Protein-adjusted differences were treated as relative site-abundance changes, not phosphorylation stoichiometry or evidence that an individual site is activating or inhibitory."""
)

for i,x in enumerate(lines):
    if x.startswith("Custom analysis code central to the conclusions"):
        lines[i]=x.replace("independent validation and figure rendering.",
                    "independent validation, IMC source-annotation reanalysis, paired phosphoproteomic sensitivity analysis and figure rendering.")
        break

replace_prefix(
"**a,** Independent validation of prioritized ligand–receptor directions using a 55-donor Reed normal-breast reference",
"""**a,** Independent validation of prioritized ligand–receptor directions using a 55-donor Reed normal-breast reference and an external primary-tumour cohort. Patient- or donor-level scores were compared using two-sided Mann–Whitney U tests with Benjamini–Hochberg correction within each tumour group. All 34 evaluable interactions were directionally concordant and FDR<0.05 across all primary tumours and in the oestrogen-receptor-positive subgroup; 33 of 34 were directionally concordant and 32 of 34 were FDR<0.05 in triple-negative breast cancer. The HER2-only subgroup did not meet the prespecified cell threshold and was not evaluated. **b,** Directional consistency in the prioritized interactions compared with all non-prioritized interactions and an abundance- and route-matched background set, illustrating that direction alone is insufficient as a validation criterion. **c,** Patient-level α-SMA imaging-mass-cytometry measurements after restricting cells to the deposited Myofibroblasts and Myofibroblasts PDPN+ phenotypes. α-SMA was higher in invasive tumours and DCIS than in normal breast, whereas the direct DCIS-versus-invasive contrast was not FDR-supported. Patients with fewer than 20 mural-like cells were excluded; two-sided Mann–Whitney U tests were corrected across the 15 marker-by-contrast tests in this definition. The Danenberg is_perivascular field denotes vessel contact and is analysed separately as a spatial-compartment sensitivity in Supplementary Figure 9 and Supplementary Table 9. **d,** Paired CPTAC mass-spectrometry measurements showing increased FAP, COL11A1 and ITGA11 protein abundance in tumour relative to matched normal tissue. Paired differences were tested using two-sided Wilcoxon signed-rank tests and corrected by the Benjamini–Hochberg method across assayed proteins. **e,** Orthogonal evaluation of a pre-specified S100A4–ERBB2 interaction prediction. Patient-level spatial association measures did not support the predicted positive relationship, including within ERBB2-positive and PAM50-HER2 tumours. This interaction was therefore excluded from the direct interaction model."""
)

MAN.write_text("\n".join(lines)+"\n")

slines=SUP.read_text().splitlines()
def sreplace_prefix(prefix,new):
    n=0
    for i,x in enumerate(slines):
        if x.startswith(prefix):
            slines[i]=new;n+=1
    if n!=1: raise RuntimeError(("supp",prefix,n))

sreplace_prefix(
"**a,** Fraction of cells in the published perivascular gate",
"""**a,** Fraction of cells carrying the Danenberg is_perivascular flag that also have explicit epithelial annotation in normal, DCIS and invasive-tumour samples. In the source annotation, is_perivascular means contact with a vessel mask and is treated here as a vessel-contacting spatial compartment, not a mural-cell lineage. **b,** Tumour-versus-normal patient-level α-SMA and caveolin-1 median differences under the original vessel-contacting annotation, after epithelial exclusion and after the stricter cytokeratin/panCK exclusion. DCIS-versus-normal and direct DCIS-versus-invasive tests for the same definitions and markers are reported in Supplementary Table 9, with Benjamini–Hochberg correction across the 15 marker-by-contrast tests within each definition. **c,** Paired tumour–normal CPTAC protein differences for FAP, COL11A1 and ITGA11 with 18 matched pairs and FDR values shown. A separate myofibroblast-defined α-SMA analysis is presented in main Figure 5c and Supplementary Table 9."""
)

insert_idx=None
for i,x in enumerate(slines):
    if x.startswith("# Supplementary Tables"):
        insert_idx=i
        break
if insert_idx is None: raise RuntimeError("supp tables header missing")
fig13=[
"## Supplementary Figure 13 | Paired CPTAC phosphoproteomic remodeling across adhesion and mechanotransduction modules",
"**a,** Raw paired tumour-minus-normal phosphosite differences plotted against same-gene total-protein-adjusted phosphosite differences for the 45 sites with at least five matched observations. Filled points denote FDR<0.05 within the predefined protein-adjusted site family; representative sites are labelled. The diagonal indicates no change after accounting for the parent-protein difference. **b,** Number of genes per predefined module with at least one testable protein-adjusted phosphosite and with at least one FDR<0.05 adjusted site. These results indicate coordinated site-level remodeling but are not interpreted as site-specific activation or inhibition and are not cell-type localized because CPTAC measurements are bulk tissue.",
"",
]
slines[insert_idx:insert_idx]=fig13

sreplace_prefix(
"**Supplementary Table 9.** Imaging mass cytometry perivascular-gate sensitivity",
"""**Supplementary Table 9.** Imaging mass cytometry sensitivity across the Danenberg vessel-contacting spatial annotation and an independently defined mural-like myofibroblast phenotype, including invasive-versus-normal, DCIS-versus-normal and direct DCIS-versus-invasive contrasts."""
)
for i,x in enumerate(slines):
    if x.startswith("**Supplementary Table 17.**"):
        slines.insert(i+1,"**Supplementary Table 18.** Predefined paired CPTAC phosphoproteomic module, including raw phosphosite tumour-minus-normal differences, same-gene total-protein-adjusted differences, experiment-wide phosphosite FDR values and predefined module-family adjusted FDR values.")
        break
else: raise RuntimeError("table17 missing")
SUP.write_text("\n".join(slines)+"\n")

for f in (ROOT/"StageM17A_VascularSpatialClinicalPhospho_v0_1/code").glob("m17a*.py"):
    shutil.copy2(f,ST/"code"/f.name)

gate={
 "status":"M18_MANUSCRIPT_SUPPLEMENT_PATCH_COMPLETE",
 "manuscript":str(MAN),
 "supplement":str(SUP),
 "key_interpretive_changes":[
  "Danenberg is_perivascular relabelled as vessel-contact spatial annotation",
  "main Figure5C uses myofibroblast-defined alpha-SMA validation",
  "caveolin-1 no longer described as mural-cell-intrinsic",
  "CPTAC phosphoproteomic module added as supplementary mechanistic evidence",
  "negative spatial/prognostic M17A analyses remain audit-only"
 ]
}
(ST/"audit/M18_manuscript_patch_gate_v0_1.json").write_text(json.dumps(gate,indent=2)+"\n")
print(json.dumps(gate,indent=2))
