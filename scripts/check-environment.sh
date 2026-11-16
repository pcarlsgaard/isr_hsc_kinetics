#!/usr/bin/env bash
set -euo pipefail
for tool in python Rscript claude gh samtools bcftools bedtools tabix fastqc; do
    command -v "$tool" >/dev/null || { echo "Missing: $tool" >&2; exit 1; }
done
.venv/bin/python -c 'import Bio, pandas, numpy, scipy, matplotlib, seaborn, requests; print("Python scientific imports OK")'
Rscript -e 'cat(R.version.string, "\n")'
claude --version
samtools --version | head -n 1 || test "${PIPESTATUS[0]}" -eq 141
