# Path remapping

The frozen scripts preserve the workstation root used for the manuscript:

`/mnt/stk35data/BreastCancer_Census_v1`

This is intentional provenance. To run the code elsewhere, either:

1. edit the top-level `ROOT`/`BASE` constants in the scripts you intend to run; or
2. create a disposable working copy and run:

```bash
python tools/remap_project_root.py /path/to/your/BreastCancer_Census_v1
```

The helper replaces only the literal historical project-root prefix in text/code files. It does not download data, change analysis logic, or modify the tagged release if you run it in a copied working tree.

For strict provenance comparisons, keep an untouched copy of the tagged release and perform remapping only in a separate working directory.
