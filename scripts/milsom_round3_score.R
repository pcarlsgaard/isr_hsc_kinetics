#!/usr/bin/env Rscript
# Round 3 (exploratory): UCell scores for the Amendment-3 sets on raw counts, annotated HSCs only.
suppressPackageStartupMessages({library(Matrix); library(SeuratObject); library(UCell); library(data.table)})
set.seed(20261116)
o <- readRDS("data/raw/zenodo_19486868/extracted/scRNAseqData/HSPCs.rds")
counts <- o@assays$RNA@counts; md <- o@meta.data; rm(o); invisible(gc())
keep <- rownames(md)[md$Celltypes_Louvain == "HSCs"]
counts <- counts[, keep]
gs <- fread("metadata/gene_sets_round3.tsv")
sets <- lapply(split(gs$gene, gs$set), function(g) intersect(unique(g), rownames(counts)))
cov <- data.table(set = names(sets), n_defined = sapply(split(gs$gene, gs$set), function(g) length(unique(g)))[names(sets)],
                  n_present = lengths(sets))
fwrite(cov, "data/processed/milsom/round3_gene_set_coverage.tsv", sep = "\t")
sc <- ScoreSignatures_UCell(counts, features = sets, maxRank = 1500, chunk.size = 200, ncores = 1, name = "")
fwrite(data.table(cell = keep, as.data.frame(sc[keep, ])), "data/processed/milsom/round3_hsc_cell_scores.tsv.gz", sep = "\t")
print(cov)
