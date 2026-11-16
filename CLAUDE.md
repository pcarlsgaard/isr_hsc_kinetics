# Project instructions

Build reproducible bioinformatics analyses of public data.

- Run analysis and dependency installation inside the devcontainer.
- Start by reading README.md and metadata/datasets.tsv. The biological question and dataset selection are still pending; do not invent them.
- Record source accessions, exact URLs, access dates, licenses, organism, genome build, checksums, and sample annotations before analysis.
- Keep raw inputs immutable; write derived data under data/processed and outputs under results/.
- Keep data, credentials, local Claude settings, and large generated artifacts out of Git.
- Use scripts or workflows for analysis; notebooks should call reusable analysis code where practical.
- Record seeds, parameters, software versions, QC decisions, and statistical assumptions. Keep exploratory results separate from confirmatory claims.
- Resolve and commit dependency locks before running scientific analyses; capture R dependencies with renv when R analysis starts.
- Never fabricate observations, accessions, results, citations, or successful checks.
- Use normal Claude permission prompts. Do not bypass permissions or alter host configuration.
