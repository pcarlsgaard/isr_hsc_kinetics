# Project instructions

Build reproducible bioinformatics analyses of public data.

- Run analysis and dependency installation inside the devcontainer.
- Start by reading docs/analysis-brief.md, README.md, and metadata/datasets.tsv. The brief defines the proposed scientific question and datasets; verify source contents and metadata before analysis. Do not invent unavailable measurements.
- Record source accessions, exact URLs, access dates, licenses, organism, genome build, checksums, and sample annotations before analysis.
- Keep raw inputs immutable; write derived data under data/processed and outputs under results/.
- Keep data, credentials, local Claude settings, and large generated artifacts out of Git.
- Use scripts or workflows for analysis; notebooks should call reusable analysis code where practical.
- Record seeds, parameters, software versions, QC decisions, and statistical assumptions. Keep exploratory results separate from confirmatory claims.
- Resolve and commit dependency locks before running scientific analyses; capture R dependencies with renv when R analysis starts.
- Never fabricate observations, accessions, results, citations, or successful checks.
- Use normal Claude permission prompts. Do not bypass permissions or alter host configuration.

## First execution stage

Before full analysis, produce docs/data-audit.md describing processed files, object classes and package versions, available count layers, species/gene identifiers, cell annotations, clone/recipient/sample mappings, kinetic parameter definitions and uncertainty, and independent biological units. Identify which requested contrasts are estimable. Download only the small source-code/metadata bundle first; inspect file sizes before downloading expression objects. This stage does not require raw FASTQs.

Define primary signatures and contrasts before inspecting associations with kinetics. Preserve original counts for pseudobulk; do not perform count-based differential expression on integrated or scaled assays. Treat transcriptomic scores as proxies and test in both directions. With five primary systems, emphasize clone-level effect sizes, plots, and leave-one-system-out sensitivity; do not fit a heavily adjusted model to five observations. Separate within-clone cell-state sensitivity checks from between-clone inference. Record gene-set versions, mouse mappings, overlap, and expressed-gene coverage. Do not redistribute upstream data or source code without checking the applicable license.
