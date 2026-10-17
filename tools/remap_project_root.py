#!/usr/bin/env python3
from pathlib import Path
import sys

OLD = "/mnt/stk35data/BreastCancer_Census_v1"

if len(sys.argv) != 2:
    raise SystemExit("Usage: python tools/remap_project_root.py /path/to/BreastCancer_Census_v1")

new = str(Path(sys.argv[1]).expanduser().resolve())
root = Path(__file__).resolve().parents[1]
suffixes = {".py", ".R", ".r", ".md", ".txt", ".json", ".tsv", ".yaml", ".yml", ".cff"}

changed = []
for p in root.rglob("*"):
    if not p.is_file() or ".git" in p.parts:
        continue
    if p.suffix not in suffixes and p.name != "CITATION.cff":
        continue
    try:
        s = p.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        continue
    if OLD in s:
        p.write_text(s.replace(OLD, new), encoding="utf-8")
        changed.append(str(p.relative_to(root)))

print(f"Remapped {len(changed)} files from {OLD} to {new}")
for x in changed:
    print(x)
