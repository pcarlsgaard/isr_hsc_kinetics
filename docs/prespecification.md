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

## Amendment 2 (2026-11-16, exploratory second round, written before any of it was run)

Status: exploratory. The first-round results in `docs/report.md` had been seen when this was written, so nothing below is confirmatory. All comparisons use annotated/sorted HSCs only.

### Change of stance on adjustment

Stemness, dormancy and growth-programme scores may be consequences of ISR-regulated translation and not confounders of it. Headline numbers in this round are therefore **unadjusted**. Stemness-adjusted values are shown beside them and read as a different question ("is there anything beyond what the stemness score already captures"), not as a correction. Gene overlap between the stemness sets and the ISR/translation sets is reported, because shared genes make adjustment partly circular.

### Additional gene sets (`metadata/gene_sets_round2.tsv`, built by `scripts/build_gene_sets_round2.py`)

ATF4 target transcripts can be raised by two routes with opposite meaning for the hypothesis: eIF2α phosphorylation (ISR; global translation down) and mTORC1/MYC growth signalling (translation up). The first-round `ISR_adaptive` set cannot tell them apart, and the hallmark mTORC1 set contains ATF4 targets. No regression of one score on another is used to separate them. Instead three non-overlapping read-outs are scored and reported side by side:

- `ISR_stress_specific`: Ddit3, Trib3, Chac1, Ppp1r15a, Atf5, Sesn2, Ddit4, Nupr1, Atf3, Gadd45a, Stc2. Targets associated with the eIF2α-P route and with ISR feedback on mTORC1 (Sesn2, Ddit4). Selection is from memory of Park 2017 (Cell Rep) and Torrence 2021 (eLife); citations and membership unverified. `ISR_stress_specific_noIEG` drops members of the dissociation immediate-early set.
- `ATF4_biosynth`: first-round `ISR_adaptive` minus Sesn2, Eif4ebp1, Cebpg, Stc2, Vegfa, Atf4 (amino-acid synthesis and transport, tRNA charging, one-carbon). Route-ambiguous.
- `Growth_noATF4`: hallmark mTORC1 minus every gene in `ISR_adaptive`, `ISR_terminal` or `ISR_stress_specific`. `RiBi_nonRP` is kept as the second growth read-out.

Readings fixed in advance: restraint by ISR predicts `ISR_stress_specific` higher together with `Growth_noATF4`/`RiBi_nonRP` lower in slower or lower-output HSCs. Growth-route ATF4 predicts `ATF4_biosynth` moving with the growth read-outs and `ISR_stress_specific` flat or following the immediate-early set.

### Milsom: HSCs compared between clones

- Cell-level standardised differences (slower minus faster, pooled cell SD, bootstrap interval over cells) for the within-run clone pairs M9.6–M3.2 (run 2) and M6.1–M7.3, M6.1–M5.1, M7.3–M5.1 (run 1), and for each clone against its same-run control. Cross-run clone pairs are not compared directly. Intervals describe sampling of cells from two particular clones; each clone is still one biological unit.
- Sets: the three above, `ISR_stress_specific_noIEG`, `RiBi_nonRP`, `RP_cytoplasmic`, first-round `ISR_adaptive`, `ISR_terminal`, `mTORC1_hallmark`, `Stem_MolO`, `Stem_dHSC`, `IEG_dissociation`, and log10 total UMIs per cell.
- Genome-wide: HSC pseudobulk (raw counts summed, log2 CPM + 1) fold change for two independent slower-versus-faster contrasts, A = M9.6 vs M3.2 and B = M6.1 vs (M7.3 + M5.1 pooled), over genes with at least 10 counts in all four pseudobulks. Spearman correlation of A with B is the read-out of a shared kinetics-associated programme. No p-value (genes are not independent). Mean fold change of each set in A and B is tabulated.

### LARRY

- Low- versus high-output clones (authors' labels, at least 5 scored HSCs, stratified by experiment, same permutation test and seed as round one) for the new sets, unadjusted; rank correlation with output per HSC with and without stemness adjustment shown together.
- "HSC-rich, low-output" clones: within each experiment, clones above the median HSC frequency and below the median output per HSC, against all other clones. Descriptive; H2 predicts the highest `ISR_stress_specific` in this group.
- Genome-wide: per-clone mean of log1p normalised expression in HSCs, per-gene within-experiment standardised difference low minus high output; rank position of each set.
- Cross-dataset: Spearman correlation of that per-gene difference with Milsom contrasts A, B and M9.6 versus its control.

Dropped from the earlier proposal: stemness-matched cell pairing (superseded by the change of stance; the clone-level adjusted value is retained for reference) and trajectory analysis.

## Amendment 3 (2026-11-16, exploratory third round, written before any of these sets was scored in HSCs)

Status: exploratory, as for Amendment 2. Purpose: replace the hand-picked round-2 ISR lists with sets taken from published data, and widen the biosynthesis read-outs. Gene lists were inspected (membership only) before this was written; no HSC score had been computed for them.

Correction to Amendment 2: the source data of Torrence et al. 2021 (eLife 10:e63326, Figure 1-source data 1) list Ddit3, Chac1, Ddit4, Sesn2 and Atf5 among the 61 ATF4 targets induced by **both** the stress route and mTORC1. Five of the eleven `ISR_stress_specific` genes were therefore not route-specific, and the round-2 result for that set is superseded by the sets below.

### Sets (`metadata/gene_sets_round3.tsv`, built by `scripts/build_gene_sets_round3.py`)

From Torrence 2021 (mouse embryonic fibroblasts; symbols as published, with current aminoacyl-tRNA synthetase and Ero1a aliases added):
- `T21_ATF4_shared`: the 61 genes flagged "Shared ISR and mTORC1 gene".
- `T21_ATF4_stress_only`: ATF4-dependent tunicamycin-induced genes, not flagged shared, not in the insulin-sensitive list, tunicamycin log2 fold change at least 1. Tunicamycin acts through PERK, so this is one ISR branch in one cell type; many members are not expressed in HSCs. No expression filter is applied (undetected genes add a constant under the rank statistic).
- `T21_mTORC1_nonATF4`: insulin-induced, rapamycin-sensitive genes that are not in the ATF4-dependent insulin list.

From human Reactome (MSigDB 2025.1.Hs), mapped through MGI one-to-one homology classes, cytoplasmic ribosomal protein genes removed: `R_ATF4_ER_stress`, `R_PERK_gene_expression`, `R_GCN2_aa_deficiency`, `R_HRI_heme_deficiency`. These mix pathway components with targets.

Biosynthesis sets from mouse MSigDB 2025.1.Mm, unmodified unless stated: amino-acid biosynthesis, serine-family biosynthesis, one-carbon metabolism, purine and pyrimidine nucleotide biosynthesis (GO BP); RNA polymerase I transcription, mitochondrial translation, cholesterol biosynthesis (Reactome); large and small ribosomal subunit biogenesis, fatty-acid biosynthesis, cytoplasmic translation without ribosomal protein genes (GO BP); oxidative phosphorylation and glycolysis (hallmark).

### Analyses and readings

Same contrasts, statistics, seed and thresholds as Amendment 2: Milsom within-run clone pairs and clone-versus-control (cell-level standardised differences); LARRY low- versus high-output clones, rank correlation with output per HSC, and HSC-rich low-output clones. Unadjusted values are the headline.

- ISR restraint predicts `T21_ATF4_stress_only` higher in slower or lower-output HSCs while `T21_mTORC1_nonATF4` is lower.
- Growth-route ATF4 predicts `T21_ATF4_shared` moving with `T21_mTORC1_nonATF4` and `T21_ATF4_stress_only` flat.
- For biosynthesis: a difference confined to ribosome and translation sets points to protein synthesis specifically; a similar difference across nucleotide, lipid and energy sets points to general biosynthetic quiet.

## Amendment 4 (2026-11-16, exploratory, written before running)

Question added: is ISR-route transcript activity a property of the HSC cell type, as opposed to a difference between HSC clones? This is the brief's "HSC-to-MPP" comparison (section 2) applied to the round-3 sets.

- Data: Milsom HSPC object, all seven mice (five clonal, two polyclonal controls). Round-3 sets and the round-1 `ISR_adaptive`, `ISR_terminal`, `Stem_MolO`, `IEG_dissociation` scored per cell (UCell on raw counts).
- Read-out: per mouse, mean score in annotated HSCs minus mean in annotated MPPs, divided by the pooled within-mouse cell SD. Mice are the units (n = 7, paired by design, so run and animal cancel). Reported: per-mouse values, median, number of mice with a positive sign, Wilcoxon signed-rank p (minimum attainable two-sided p = 0.016). Repeated for HSCs versus all non-HSC progenitors pooled.
- Reading: cell-type ISR activity predicts `T21_ATF4_stress_only` higher in HSCs than MPPs in most mice while `T21_mTORC1_nonATF4` is lower. Clonal systems are not normal haematopoiesis; the two controls are shown separately.
