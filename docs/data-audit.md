# Data audit

Audit date 2026-11-16. Provenance, URLs, licenses and checksums are in `metadata/datasets.tsv`. Audit logs are regenerated under `results/audit/` by `scripts/audit_seurat_object.R`, `scripts/milsom_score_hspc.R` and `scripts/larry_clone_scores.py`.

## Primary: Milsom lab clonal systems (Zenodo 19486868)

Paper: Nat Cell Biol, DOI 10.1038/s41556-026-01958-0, published 2026-05-25. A correction (10.1038/s41556-026-02103-7, 2026-10-16) exists and was not reviewed. License of the deposit: DKFZ Academic Non-Commercial Research License; nothing from it is redistributed here.

| File | Size | Status |
|---|---|---|
| `manuscript.zip` | 5.0 MB | Inspected: figure scripts, kinetic fits, FACS chimerism tables, stemness gene sets |
| `scRNAseqData/HSPCs.rds` | 690 MB | Audited and analysed |
| `scRNAseqData/scBM_atlas.rds` | 1.24 GB | Not loaded (7 GB container memory; not needed for HSC-level questions) |
| `scRNAseqData/slingshot.rds` | 516 MB | Not loaded |

**Object.** `HSPCs.rds` is a Seurat v3.2.2 object (read with SeuratObject 5.4.0 by slot access; the Seurat package is not installed). One assay, `RNA`: `counts` 31,053 genes x 29,730 cells, integer-valued, original UMI counts; `data` log-normalised; `scale.data` 3,000 genes. Reductions `iNMF_raw`, `integrated`, `umap` come from LIGER integration and were not used. Gene identifiers are mouse symbols (no Ensembl IDs), 13 `mt-` genes, 80 cytoplasmic ribosomal protein genes by the pre-specified pattern. Genome build GRCm38, Cell Ranger 3.1 (paper Methods).

**Cell annotations.** `Celltypes_Louvain` (8 HSPC states including `HSCs` and `MPPs`), `Celltype` (13 states including an `HSCs end` subset), `Phase`, `S.Score`, `G2M.Score`, `percent.mt`, `percent.ribo`, `Mice`, `Run`, `Batch`. All cells are cKit-enriched bone marrow at week 20.

**Independent biological units.** Seven mice: five single-HSC-derived systems and two polyclonal controls (30 LT-HSCs each). One mouse per system; no replicate animals per clone.

| System (object ID) | Kinetic order | Run | Cells | Annotated HSCs | `HSCs end` |
|---|---|---|---|---|---|
| M9.6 (`M96`) | 1, slowest | 2 | 3,616 | 825 | 824 |
| M6.1 (`M61`) | 2 | 1 | 3,854 | 147 | 119 |
| M3.2 (`M32`) | 3 | 2 | 6,954 | 243 | 188 |
| M7.3 (`M73`) | 4 | 1 | 3,101 | 26 | 14 |
| M5.1 (`M51`) | 5, fastest | 1 | 1,373 | 25 | 16 |
| Control `M14` | – | 1 | 4,792 | 205 | 185 |
| Control `M16` | – | 2 | 6,040 | 371 | 405 |

- Sequencing run is partly confounded with kinetics: the slowest system shares a run only with M3.2 and one control.
- The object has no donor-reporter flag or transgene feature, so clone origin of individual cells cannot be verified from the data. The paper describes the sequenced cells as clonally derived progeny. M5.1 has 0% phenotypic HSC chimerism by flow yet 25 transcriptionally annotated HSCs, so "HSC" here is a transcriptional label.

**Kinetic parameters.** Per lineage, peripheral-blood chimerism over weeks 4–20 is fitted with `y(t) = R(t−t0) / (1 + ((t−t0)/θ)^n)` (`scTransplant_functions.R`). Derived: `t0` (delay), `tyMax`, `tGrowth`, `tDecline`, `tHalfReg`, `AUC`. Uncertainty: only per-fit RMSE is provided; no intervals on parameters. For the three slowest sequenced systems `tDecline` and `tHalfReg` take identical capped values (39.06 and 42.74 weeks) and several fits sit on parameter bounds, so these are censored. The kinetic variable is therefore used as the authors' rank order, with `tyMax` as an uncapped check.

| System | tyMax (wk) | HSC share, primary (%) | HSC share, secondary (%) | Mean PB AUC | Blood output pattern |
|---|---|---|---|---|---|
| M9.6 | 19.8 | 98.8 | 69.8 | 1.17 | platelets ~73%, other lineages ≤ ~9% |
| M6.1 | 19.4 | 59.9 | 49.1 | 3.80 | platelet, myeloid, erythroid; no lymphoid |
| M3.2 | 18.3 | 3.8 | 0.4 | 7.70 | high multilineage |
| M7.3 | 16.4 | 0.9 | 0 | 3.96 | multilineage |
| M5.1 | 10.4 | 0 | 0 | 1.98 | early, transient |

Kinetic rank, HSC share and lineage restriction are almost perfectly collinear across the five systems.

**Estimable contrasts.** Between-system ordering of HSC-level signature means (n = 5, exact permutation; minimum two-sided p = 0.017); within-run ordered pairs (4); comparison with the same-run control; within-system cell-level correlations. Not estimable: any separation of kinetics from HSC share or lineage restriction; curve shape; adjusted models; clone-replicate variance.

## Validation: LARRY barcoding (GSE134242; Tcf15_HSCs repository)

- GEO holds 90 per-library raw count tables (368 MB tar). They are not cell-aligned to the published annotations, so they were not used.
- The repository (commit `9747d71`, no license) provides per-experiment clone tables, barcode lists, sort-gate flags and clusters. The authors' filtered, normalised matrices (527 MB, Dropbox link in the README) are cell-aligned to those files: 16,976 x 25,289; 10,375 x 38,458; 4,845 x 38,458 (cells x genes, mouse symbols). Values are non-integer, so only rank-based scoring is valid; no count-based differential expression was run.
- Clone tables give, per clone: total cells, HSC and non-HSC counts, `relout` (non-HSC frequency divided by HSC frequency, i.e. output per HSC), lineage fractions, `inactive`/`active` (low/high output) calls and Mk bias. Output per HSC is therefore available as published.
- Persistence over time and reconstitution kinetics are not measured in these tables.
- The published clone tables cannot be reproduced from the deposited barcode lists: cells carry several concatenated barcodes and per-clone counts do not match under exact or containment matching. Following the original notebook, a cell is assigned to a clone only if its barcode string equals the clone key. This keeps 2,052 + 1,331 + 183 sorted-HSC cells in 51 + 75 + 30 clones; 84 clones have at least 5 scored HSCs.
- The Sample 1 barcode list is four entries shorter than its matrix (16,972 vs 16,976). Alignment from the first cell is assumed, as the notebook does implicitly; this cannot be verified.
- Independent units are clones nested in three experiments (recipient/experiment batch). Experiment 3 contributes only 8 clones at the 5-cell threshold.

## Gene sets

MSigDB mouse 2025.1.Mm (native symbols, no ortholog mapping) plus curated lists in `metadata/gene_sets_curated.tsv`. Literature citations attached to the curated ISR, 5'TOP and immediate-early lists were written from memory and are not yet checked against the primary papers.

| Set | Defined | In Milsom object | In LARRY exp. 1 / 2 |
|---|---|---|---|
| ISR adaptive | 52 (incl. old/new tRNA-synthetase aliases) | 39 | 38 / 39 |
| ISR terminal | 16 | 15 | 15 / 15 |
| Ribosome biogenesis, non-RP | 285 | 273 | 256 / 263 |
| Cytoplasmic RP genes | regex | 80 | 81 / 81 |
| mTORC1 hallmark | 199 | 188 | 189 / 190 |

Overlaps that matter: 10 of 39 ISR-adaptive genes are also mTORC1 hallmark genes (*Asns, Psat1, Phgdh, Psph, Shmt2, Mthfd2, Slc7a5, Slc7a11, Slc1a5, Cth*); four ISR-terminal genes (*Txnip, Nupr1, Trib3, Atf3*) sit in stemness or immediate-early sets. Full matrix: `results/milsom/gene_set_overlap.tsv`.

## Software

R 4.3.3, SeuratObject 5.4.0, UCell 2.6.2, Matrix, data.table (`renv.lock`). Python 3.12 with numpy, pandas, scipy, matplotlib (`requirements.lock`). Seeds: 20261116. UCell `maxRank` 1500 in both datasets (reimplemented in Python for LARRY with the same statistic).
