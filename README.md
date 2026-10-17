# Reference-anchored stromal remodeling in invasive breast cancer — analysis code v0.98

This repository contains the manuscript-specific analysis and audit code associated with:

**Reference-anchored decomposition of stromal redistribution and reprogramming in invasive breast cancer**

Version: **v0.98**  
Release date: **17 October 2026**  
License: **BSD 3-Clause**

## Scope

The repository supports the central analyses used to distinguish:

1. redistribution among normal-reference stromal states; and
2. transcriptional reprogramming within corresponding reference affinities.

It includes code for patient-level stromal composition, fine-state mapping/confidence and degradation calibration, patient-state pseudobulk construction, differential expression, pathway analysis, ligand–receptor prioritization, independent validation, imaging-mass-cytometry reanalysis, paired phosphoproteomic sensitivity analysis, method-concordance testing, and the pre-specified longitudinal clinical-trial stress tests in FELINE and BREAKFAST.

## Repository structure

- `code/` — core manuscript analysis and figure-rendering scripts.
- `M20B_MethodConcordance/` — edgeR versus limma-voom concordance analyses.
- `M22_ClinicalTrialValidation/` — pre-specified longitudinal clinical-trial analyses, audit tables and summary outputs.
- `environment/` — frozen Python/R/system version information.
- `documentation/` — workflow maps and scientific checkpoints.
- `audit/` — final scientific gate and decision-ledger records.
- `PATH_REMAP.md` — instructions for adapting the original workstation paths.
- `tools/remap_project_root.py` — helper for creating a local path-remapped working copy.
- `DATA_ACCESS.md` — data-access and redistribution notes.
- `OMITTED_PUBLIC_FILES.txt` — reviewer-only derived files intentionally not redistributed.

## Frozen-code principle

The scripts are preserved from the manuscript analysis environment. Many record the original project root:

`/mnt/stk35data/BreastCancer_Census_v1`

This path is provenance, not a requirement. Use the supplied path-remapping helper or edit the root constants to point to your local data tree. The helper modifies a *working copy*; the tagged release remains unchanged.

## Historical fine-state mapping limitation

The original standalone executable used for the historical correlation-based fine-state mapping was not retained. Frozen cell-level assignments and downstream tables are the authoritative historical outputs. Where a deterministic reconstruction/audit implementation is supplied, it is explicitly an audit/reconstruction implementation and must not be represented as the original historical executable.

## Data

Primary datasets are **not redistributed** in this repository. They are previously published public or controlled-access datasets and must be obtained from their original repositories under the applicable terms. Accessions, cohort definitions and provenance are documented in the manuscript, Supplementary Tables and resource catalogue.

Two large cell-level derived mapping exports from the reviewer package are also intentionally omitted from the public release. They can be regenerated from the public source datasets using the supplied scripts; see `OMITTED_PUBLIC_FILES.txt`.

## Reproducibility

The manuscript uses the patient/donor as the biological replicate whenever identity is resolvable. Cells are not treated as independent biological replicates for inferential testing.

Environment records are provided in `environment/`. The release also contains a SHA256 manifest so users can verify file integrity.

## Recommended reproduction order

The complete historical pipeline spans several frozen stages. For review or reproduction:

1. obtain the cited source datasets;
2. reproduce the normal-reference and Wave-A input objects described in the manuscript;
3. follow `documentation/M18_ACTIVE_ANALYSIS_MAP.md` for the principal analysis order;
4. run the method-concordance module in `M20B_MethodConcordance/`;
5. run the longitudinal validation module in `M22_ClinicalTrialValidation/` when reproducing the clinical-trial sensitivity analyses;
6. compare generated outputs with the retained summary/audit tables and SHA256 manifest.

## Public release identity

Git tag: `v0.98`  
Git commit: the exact tagged commit can be obtained with `git rev-parse v0.98`.  
Permanent archive DOI: **to be added after Zenodo deposition**.

## Citation

See `CITATION.cff`. Once the permanent archive DOI is minted, cite both the manuscript and the archived software release.

## Contact

Pankaj Goyal  
Department of Biotechnology, Central University of Rajasthan  
pankaj_bio@curaj.ac.in
