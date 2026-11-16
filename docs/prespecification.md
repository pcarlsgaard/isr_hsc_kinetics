# Pre-specification and testability statement

Written 2026-11-16, after inspecting the kinetic/FACS tables in the Zenodo `manuscript.zip` bundle and **before** opening any expression object. Everything below was fixed without knowledge of gene expression in either dataset.

## Units and variables (Milsom / Zenodo 19486868)

- Independent biological units: five single-HSC-derived systems (author IDs M9.6, M6.1, M3.2, M7.3, M5.1) and two polyclonal controls (M14, M16). Cells are not replicates.
- **Kinetic rank** (primary exposure): the authors' slow-to-fast order 9.6 < 6.1 < 3.2 < 7.3 < 5.1. The fitted `tDecline`/`tHalfReg` are capped at identical values for the three slowest systems, so they are not used as continuous variables. `tyMax` (uncapped; 19.8, 19.4, 18.3, 16.4, 10.4 weeks) is the continuous sensitivity variable.
- **Self-renewal proxy**: clone share of phenotypic HSCs after primary transplant (98.8, 59.9, 3.8, 0.9, 0 %), with the secondary-transplant value as support.
- **Output**: mean fitted peripheral-blood AUC (1.17, 3.80, 7.70, 3.96, 1.98 in kinetic order). Output per HSC is not computable for systems with ~0 HSC share and is not used.
- Known structure, seen before expression: kinetic rank, HSC share and lineage restriction are almost perfectly confounded across the five systems; output is non-monotone in kinetic rank (peak at M3.2).

## Cell population

Cells annotated as HSCs by the original authors, in the HSPC object if it carries that label. A system contributes to HSC-level comparisons only with at least 20 annotated HSCs; systems below that are reported as not estimable. If cells in clonal systems are not exclusively clone-derived, that is reported as a limitation, not corrected.

## Signatures

Definitions are in `metadata/gene_sets_curated.tsv` and `metadata/gene_sets_msigdb.tsv` (MSigDB mouse 2025.1.Mm, native mouse symbols, no ortholog mapping). Sets are restricted to genes present in the object; coverage and pairwise overlap are reported.

Primary (five, Bonferroni/BH family): `ISR_adaptive`, `ISR_terminal`, `RiBi_nonRP`, `RP_cytoplasmic`, `mTORC1_hallmark`.
Secondary (exploratory): UPR sets with ISR genes removed, MYC targets V1/V2, rRNA processing, translation initiation/elongation, non-RP 5'TOP, GO ISR signaling components.
Controls: authors' stemness/quiescence sets shipped in the bundle (MolO, dormant-vs-active HSC, Rodriguez-Fraticelli low-output), proliferation (E2F, G2M hallmarks), dissociation immediate-early genes.

## Analyses

1. Per-cell UCell scores on the original count layer; per-system mean and median among HSCs.
2. Pseudobulk: sum raw counts over HSCs per system, log2 CPM, per-gene values for the primary sets.
3. Between-system inference (n = 5): Spearman correlation of per-system score with kinetic rank, reported with the exact permutation p-value (minimum attainable two-sided p = 1/60 ≈ 0.017), leave-one-system-out range, and the polyclonal controls plotted as reference. No regression adjustment across systems.
4. Effect size: difference between the two slowest (M9.6, M6.1) and the remaining systems, in units of the pooled within-system cell-level SD.
5. Within-system sensitivity (not between-clone inference): cell-level Spearman correlation of each primary score with stemness, proliferation, IEG and library size, per system; and system differences recomputed within strata of the stemness score.
6. Both directions are tested for every signature.

## Predictions

- H1: `ISR_adaptive` higher and `RiBi_nonRP`, `RP_cytoplasmic`, `mTORC1_hallmark` lower in slower systems; `ISR_terminal` not elevated.
- H2: the slowest, HSC-dominated, low-output system (M9.6) additionally shows elevated `ISR_terminal`, or lies beyond the range of the other slow system and the controls on the restraint signatures. Signatures peak or trough non-monotonically against output.

## What these data can and cannot resolve (stated in advance)

- Can: whether translation-related transcriptional programs differ between systems in the predicted direction, their size relative to within-system variation, and whether differences survive stratification by stemness state.
- Cannot: separate kinetics from HSC share or lineage restriction (confounded by design at n = 5); establish curve shape (H2) beyond a single extreme clone; measure protein synthesis or eIF2α phosphorylation; establish causality or direction. A perfect monotone trend across five systems is the strongest possible result and still corresponds to p ≈ 0.017 for one signature before multiplicity.
- LARRY (GSE134242) is used only as a between-clone consistency check on output and persistence, with its own pre-specified contrast (low-output vs high-output HSC clones as annotated by the original authors; tertile test of shape), and cannot test kinetics.

## Amendment 1 (2026-11-16, after the object metadata audit, before any expression score was computed)

The HSPC object shows a run structure not known when the plan was written: Run 1 holds M6.1, M7.3, M5.1 and control M14; Run 2 holds M9.6, M3.2 and control M16. Kinetic rank is therefore partly confounded with sequencing run. Added, without changing the primary analysis:

- Every between-system summary is also reported as a difference from the same-run polyclonal control (control-referenced score).
- Within-run ordering is reported separately: Run 1 gives three ordered pairs (M6.1 slower than M7.3 slower than M5.1), Run 2 gives one (M9.6 slower than M3.2). Concordance over these four within-run pairs is the batch-robust read-out; it has no usable p-value.
- HSC counts are 825, 147, 243, 26 and 25 for M9.6, M6.1, M3.2, M7.3 and M5.1; the two fastest systems are just above the inclusion threshold and their means are reported with bootstrap intervals over cells.
- The deposited object carries no donor-reporter flag, so whether each annotated HSC is clone-derived cannot be checked from the data. The paper describes the sequenced cells as clonally derived bone-marrow progeny.

LARRY convention fixed at the same time, before LARRY scores were inspected: a cell is assigned to a clone only when its barcode string exactly equals the clone key, as in the original notebook (cells carrying several barcodes are excluded because the published clone tables cannot be reproduced from the deposited barcode lists). Clones need at least 5 scored sorted-HSC cells (sensitivity: 3 and 10). Analyses are stratified by experiment.
