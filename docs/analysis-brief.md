# Project: Translational regulation and HSC reconstitution kinetics

## Scientific question

Investigate whether slowly reconstituting, high-self-renewal hematopoietic stem cells (HSCs) exhibit distinctive programs of integrated stress response (ISR), ribosome biogenesis, and translational regulation, and whether these programs correlate with apparent megakaryocyte/myeloid lineage bias.

**Working hypothesis:** Slow-kinetics HSCs maintain a translational state characterized by tightly regulated protein synthesis, potentially involving elevated ATF4/ISR signaling, reduced ribosome biogenesis, and altered translation initiation or elongation. These programs may be associated with long-term self-renewal and delayed differentiation.

Treat this as an exploratory, bidirectional hypothesis, not an assumed result.

## Datasets

**Primary — Milsom et al., Nature Cell Biology 2026**

- Paper: https://doi.org/10.1038/s41556-026-01958-0
- Data/code: https://zenodo.org/records/19486868
- Analyze the HSC/HSPC compartments from the five single-HSC-derived hematopoietic systems, ordered by experimentally measured reconstitution kinetics, plus appropriate polyclonal controls.
- Retrieve the actual kinetic parameters and original clone annotations. Do not equate differentiation pseudotime with engraftment kinetics.

**Validation — Rodríguez-Fraticelli et al., Nature 2020**

- GEO: GSE134242
- Code: https://github.com/rodriguez-fraticelli/Tcf15_HSCs
- Use LARRY clonal annotations to examine whether transcriptomic programs correlate with long-term HSC persistence, self-renewal, and lineage output.
- Recognize that endpoint lineage bias is not a direct measurement of reconstitution kinetics.

Use processed data where possible. Avoid raw FASTQ processing unless essential.

## Analysis

### 1. Define translation-related programs

Obtain curated, nonredundant gene sets from MSigDB, GO, Reactome, and primary literature:

- ISR and ATF4 downstream transcriptional targets.
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

Perform sensitivity analyses controlling for cell-cycle status, quiescence, sequencing depth, mitochondrial/ribosomal transcript fractions, and HSC-state composition.

Avoid treating thousands of cells from five clonal systems as independent biological replicates. Use clone-level summaries or pseudobulk methods and describe the limits of statistical inference.

### 3. Validate using LARRY barcoding

In GSE134242:

- Reconstruct clone identities, clone sizes, HSC persistence, and available lineage outputs using the original analysis conventions.
- Examine relationships between the translation-related signatures and clonal HSC persistence, self-renewal, and Mk/myeloid versus lymphoid output.
- Prefer HSC-state expression profiles rather than comparing mature myeloid and lymphoid cells directly.
- Account for recipient/batch effects, clone size, and sampling uncertainty.
- Determine whether any observed associations persist after controlling for established quiescence/stemness signatures.

If the available barcoded cells cannot support a particular comparison, document that limitation instead of substituting an inappropriate proxy.

### 4. Synthesize the findings

Evaluate three potentially distinct relationships:

A. **ISR activation and slow kinetics:** Do slow HSCs express stronger ISR-associated transcriptional programs?

B. **Ribosome/translation regulation and slow kinetics:** Do slow HSCs exhibit altered ribosome biogenesis or translation machinery, independent of quiescence?

C. **Translation programs and apparent lineage bias:** Are the same programs associated with long-term HSC persistence and Mk/myeloid output, or is the association explained by differentiation stage and clonal competition?

Distinguish results supported by direct functional readouts, clonal proxies, and simple transcriptional correlations.

## Deliverables

Produce a concise exploratory report, not a manuscript or comprehensive literature review.

1. **One-page executive summary:** Findings, effect sizes, biological interpretation, limitations, and confidence.
2. **Figure 1:** Pathway activity across slow-to-fast HSC systems.
3. **Figure 2:** Heatmap of key ISR, ribosome biogenesis, and translation-regulation genes across HSC states.
4. **Figure 3:** Relationship between translation signatures and LARRY clone fate/persistence, if adequately supported.
5. **Results table:** Pathway, direction of association, effect size, statistical evidence, confounding assessment, and consistency across datasets.
6. **Reproducible code:** R/Seurat preferred, with documented packages, gene sets, inputs, and scripts.

Use appropriate multiple-testing correction, report sample sizes at the clone and cell levels, and distinguish exploratory patterns from statistically supported findings.

## Execution constraints

Prioritize speed, reproducibility, and biological interpretability. Start with processed objects. Reuse published annotations. Do not train a predictive model or conduct a broad multi-dataset search at this stage.

Complete the primary Milsom analysis first, then use LARRY as a focused validation. If data accessibility or sample sizes limit inference, finish with the available results and document the limitation.

Conclude with the three most informative experimental tests, prioritizing ATF4 translational reporter activity, O-propargyl-puromycin protein synthesis measurements, and HSC functional behavior under ISR perturbation. Consider how a future WT versus Tet2-deficient comparison could distinguish intrinsic translational regulation from altered differentiation kinetics.
