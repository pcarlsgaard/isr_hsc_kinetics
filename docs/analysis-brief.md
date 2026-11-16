# Project: Translational regulation and HSC reconstitution kinetics

## Scope and context

This project is a computational reanalysis of published, publicly deposited single-cell transcriptomic data from mouse hematopoietic stem and progenitor cells. It is basic research into how stem cells balance self-renewal against blood production, a question relevant to both normal hematopoiesis and to understanding how that balance is disturbed in age-related clonal hematopoiesis and myeloid neoplasia.

Constraints that apply throughout:

- All inputs are processed objects from public repositories, analyzed under their stated licenses. No new experiments are performed and no laboratory protocols are produced by this work.
- All readouts are transcriptional. Pathway scores are proxies for pathway activity, never measurements of protein synthesis rate, self-renewal capacity, or disease state. Every claim in the report must carry the proxy qualifier.
- Where disease relevance is discussed, it is interpretive framing for basic mechanism, not a clinical, diagnostic, or prognostic claim.
- Sample sizes here are small (five clonal systems in the primary dataset). Conclusions are exploratory and hypothesis-generating by construction.

## Scientific question

Investigate whether hematopoietic stem cells (HSCs) with slow reconstitution kinetics and high self-renewal exhibit distinctive programs of integrated stress response (ISR), ribosome biogenesis, and translational regulation; and whether the same programs relate to how much blood output a clone generates and to its lineage composition.

The motivating idea has two parts, which are tested separately because they make different predictions.

**H1 — translational restraint supports stemness.** HSCs with slow kinetics and high self-renewal maintain a state of tightly regulated protein synthesis, potentially involving elevated ATF4/ISR signaling, reduced ribosome biogenesis, and altered translation initiation or elongation. The current datasets can address this.

**H2 — the same program taken too far trades output for self-renewal.** Beyond some level, translational restraint or ISR activation may accompany a state that retains or expands HSCs while contributing little mature blood. The available healthy datasets can only partially address this; see Analysis section 5 and the limitations deliverable.

Both are exploratory and bidirectional. H1 predicts a monotone association with kinetics; H2 predicts a non-monotone relationship with output, in which intermediate restraint is most productive. A simple linear correlation cannot distinguish them, and prior work found that moving HSC protein synthesis either up or down impaired function (Signer et al., Nature 2014), so test for both shapes rather than assuming one.

Do not build into the premise that slow reconstitution and high self-renewal are the same property. Whether they coincide in these five systems is a question for the audit, to be answered from persistence or secondary-transplant data if available, and reported as unresolved otherwise.

## Datasets

**Primary — Milsom et al., Nature Cell Biology 2026**

- Paper: https://doi.org/10.1038/s41556-026-01958-0
- Data/code: https://zenodo.org/records/19486868
- Analyze the HSC/HSPC compartments from the five single-HSC-derived hematopoietic systems, ordered by experimentally measured reconstitution kinetics, plus appropriate polyclonal controls.
- Retrieve the actual kinetic parameters and original clone annotations. Do not equate differentiation pseudotime with engraftment kinetics.
- Position relative to the source study: the deposit's title attributes skewed lineage output to extrinsic regulation, while the study treats reconstitution kinetics as a clone-intrinsic, heritable property. This project targets the kinetics axis. The claim under test is that intrinsic kinetic differences are linked to ISR, protein synthesis, or ribosome biogenesis programs; lineage skew is a secondary readout that may be extrinsically driven and is not assumed to have an intrinsic transcriptional explanation. Confirm the authors' intrinsic-kinetics evidence (e.g. secondary-transplant heritability) from the paper and source data before relying on it.

**Validation — Rodríguez-Fraticelli et al., Nature 2020**

- GEO: GSE134242
- Code: https://github.com/rodriguez-fraticelli/Tcf15_HSCs
- Use LARRY clonal annotations to examine whether transcriptomic programs correlate with long-term HSC persistence, self-renewal, and lineage output.
- Recognize that endpoint lineage bias is not a direct measurement of reconstitution kinetics.

Use processed data where possible. Avoid raw FASTQ processing unless essential.

## Analysis

### 1. Define translation-related programs

Obtain curated, nonredundant gene sets from MSigDB, GO, Reactome, and primary literature:

- ISR and ATF4 downstream transcriptional targets. Split these into two subsets, because H1 and H2 make different predictions about them:
  - **Adaptive arm:** ATF4 targets supporting amino-acid supply, tRNA charging, one-carbon metabolism, and redox balance.
  - **Terminal arm:** DDIT3/CHOP and its downstream targets, including PPP1R15A/GADD34 and pro-apoptotic effectors.
  - H1 predicts slow, high-self-renewal HSCs sit high on the adaptive arm without elevation of the terminal arm. Joint elevation of both arms would be the signature H2 expects of a state pushed past the productive range. Report the two arms separately in every figure and table; do not merge them into a single ISR score.
- Ribosome biogenesis, rRNA processing, ribosome assembly.
- Cytoplasmic translation and translation initiation/elongation.
- mTORC1, 5′TOP-mediated translation regulation, EIF4E/4EBP signaling.
- Proteostasis, unfolded protein response, and protein quality control.
- Optional: HRI/EIF2AK1-specific and GCN2/EIF2AK4-associated programs, if sufficiently well characterized.

Keep canonical ISR targets separate from generic ER stress/UPR signatures. Do not use ATF4 mRNA abundance alone as a proxy for ISR activation, or ribosomal mRNA abundance as a direct measurement of protein synthesis rates.

### 2. Relate pathway activity to HSC kinetics

Within phenotypically/transcriptionally defined HSCs:

- Score pathways using UCell, AUCell, or an equivalent robust method.
- Compare pathway scores across slow-to-fast reconstitution systems.
- Assess continuous associations with measured kinetic parameters where available.
- Examine whether signatures change systematically along HSC-to-MPP differentiation trajectories.
- Identify genes and pathways that distinguish slow from fast HSC systems.
- Record, alongside the kinetic parameters, whatever output measures the dataset provides: chimerism over time, per-lineage contribution, and secondary-transplant results. Output magnitude and lineage composition are separate axes and must not be collapsed into one another.

Perform sensitivity analyses controlling for cell-cycle status, quiescence, sequencing depth, mitochondrial/ribosomal transcript fractions, and HSC-state composition.

Avoid treating thousands of cells from five clonal systems as independent biological replicates. Use clone-level summaries or pseudobulk methods and describe the limits of statistical inference.

Do not fit curves or test for non-monotone relationships in the Milsom systems. Five observations cannot support it. Report clone-level effect sizes, per-system plots, and leave-one-system-out sensitivity, and note the shape question as unresolved at this sample size.

### 3. Validate using LARRY barcoding

In GSE134242:

- Reconstruct clone identities, clone sizes, HSC persistence, and available lineage outputs using the original analysis conventions.
- Define two distinct per-clone quantities and analyze them separately:
  - **Output magnitude:** mature-cell contribution relative to the clone's HSC contribution, i.e. blood produced per stem cell. This is the H2-relevant axis.
  - **Lineage composition:** Mk/myeloid versus lymphoid share of that output. This is the H1/section-4C axis.
- Make the low-output versus high-output clone contrast the primary LARRY comparison, conditional on the published clone-fate annotations supporting it. Low-output long-lived clones in this dataset are the closest available healthy analogue of the state H2 describes; the original study associated such clones with Tcf15. Confirm in the data audit that output and persistence are separably annotated before committing to this comparison.
- Prefer HSC-state expression profiles rather than comparing mature myeloid and lymphoid cells directly.
- Account for recipient/batch effects, clone size, and sampling uncertainty.
- Determine whether any observed associations persist after controlling for established quiescence/stemness signatures.
- Pre-specify one test of relationship shape before inspecting results: group clones into tertiles of each signature score and compare output magnitude and persistence across tertiles. Record the test and the grouping rule in the analysis script before running it. A monotone tertile pattern supports H1; a peak in the middle tertile is the H2 prediction.

If the available barcoded cells cannot support a particular comparison, document that limitation instead of substituting an inappropriate proxy. In particular, if output per HSC cannot be derived from the deposited annotations, say so and drop the H2 arm of the LARRY analysis rather than approximating it with clone size alone.

### 4. Synthesize the findings

Evaluate three potentially distinct relationships:

A. **ISR activation and slow kinetics:** Do slow HSCs express stronger ISR-associated transcriptional programs?

B. **Ribosome/translation regulation and slow kinetics:** Do slow HSCs exhibit altered ribosome biogenesis or translation machinery, independent of quiescence?

C. **Translation programs and apparent lineage bias (secondary):** The source study attributes skewed lineage output to extrinsic regulation. Ask only whether any association between these programs and Mk/myeloid output remains after accounting for kinetics, differentiation stage, and clonal competition. A null result here is compatible with the main hypothesis.

D. **Translation programs and output magnitude:** Is signature score related to blood output per HSC, and is that relationship monotone or peaked? Does the adaptive ISR arm behave differently from the terminal arm on this axis?

Distinguish results supported by direct functional readouts, clonal proxies, and simple transcriptional correlations.

State explicitly which of H1 and H2 each result speaks to, and which parts of H2 remain untested.

### 5. Optional contrast stage for H2

Only after sections 1–4 are complete, and only if the earlier results motivate it. H2 requires a state in which restraint or ISR is elevated beyond the healthy range, which five healthy clonal systems cannot provide.

- Add at most one focused public dataset. Candidates to evaluate during the audit: aged mouse HSCs (an expanded pool with reduced output per cell and a myeloid/platelet lean), or a Tet2- or Dnmt3a-mutant HSC comparison.
- Name the chosen accession in `metadata/datasets.tsv` with the same provenance fields as the primary datasets.
- Keep this stage out of the primary deliverables and label it exploratory. A single added dataset with different batch structure cannot establish a dose relationship; it can only show whether the adaptive and terminal arms separate in the predicted direction.
- Do not expand this into a multi-dataset search. If no suitable dataset is identified quickly, record H2 as untestable with current public data and stop.

## Deliverables

Produce a concise exploratory report, not a manuscript or comprehensive literature review.

1. **One-page executive summary:** Findings, effect sizes, biological interpretation, limitations, and confidence.
2. **Figure 1:** Pathway activity across slow-to-fast HSC systems.
3. **Figure 2:** Heatmap of key ISR, ribosome biogenesis, and translation-regulation genes across HSC states.
4. **Figure 3:** Relationship between translation signatures and LARRY clone fate/persistence, if adequately supported.
5. **Results table:** Pathway, direction of association, effect size, statistical evidence, confounding assessment, and consistency across datasets. Report the adaptive and terminal ISR arms as separate rows.
6. **Reproducible code:** R/Seurat preferred, with documented packages, gene sets, inputs, and scripts.
7. **Testability statement:** A short section written and committed *before* results are inspected, recording which parts of H1 and H2 the available data can and cannot address, and the pre-specified tests from sections 2 and 3.

Use appropriate multiple-testing correction, report sample sizes at the clone and cell levels, and distinguish exploratory patterns from statistically supported findings.

Carry the proxy qualifier into the report's wording. Transcriptional ISR and ribosome scores do not measure translation rate, so "restrained translation" remains an inferred state throughout. Do not describe any result as evidence about disease causation; the datasets are healthy-donor-derived mouse transplant experiments and the strongest available statement concerns association within normal hematopoiesis.

## Execution constraints

Prioritize speed, reproducibility, and biological interpretability. Start with processed objects. Reuse published annotations. Do not train a predictive model or conduct a broad multi-dataset search at this stage.

Complete the primary Milsom analysis first, then use LARRY as a focused validation. If data accessibility or sample sizes limit inference, finish with the available results and document the limitation.

Conclude by naming the three most informative experimental tests that would follow from the computational results, at the level of study design rather than protocol. Priorities: ATF4 translational reporter activity, O-propargyl-puromycin protein synthesis measurement, and HSC functional behavior under ISR perturbation.

For H2, the informative design is graded rather than binary: ISR perturbation across several levels, with blood output per HSC as the primary readout rather than HSC retention alone. A two-arm on/off experiment cannot distinguish "restraint supports stemness" from "restraint past a threshold costs output."

Consider how a future WT versus Tet2-deficient comparison could distinguish intrinsic translational regulation from altered differentiation kinetics.

These are proposed directions for the reader's own experimental planning, stated as designs and readouts. This document does not specify laboratory procedures.
