#!/usr/bin/env Rscript
# Milsom HSPC object: UCell scores on raw counts for the pre-specified gene sets, plus
# pseudobulk raw-count sums per system x annotated cell type. Slot access only (legacy v3 object).
suppressPackageStartupMessages({library(Matrix); library(SeuratObject); library(UCell); library(data.table)})
set.seed(20261116)
inp <- "data/raw/zenodo_19486868/extracted/scRNAseqData/HSPCs.rds"
out <- "data/processed/milsom"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
o <- readRDS(inp)
counts <- o@assays$RNA@counts
md <- o@meta.data
stopifnot(identical(rownames(md), colnames(counts)))
rm(o); invisible(gc())
gs <- fread("data/processed/gene_sets.tsv")
rp <- grep("^(Rp[ls][0-9]+[a-z]?|Rplp[0-2]|Rpsa)$", rownames(counts), value = TRUE)
sets <- split(gs$gene, gs$set); sets$RP_cytoplasmic <- rp
cov <- rbindlist(lapply(names(sets), function(s) data.table(
  set = s, n_defined = length(unique(sets[[s]])), n_present = sum(unique(sets[[s]]) %in% rownames(counts)))))
sets <- lapply(sets, function(g) intersect(unique(g), rownames(counts)))
fwrite(cov, file.path(out, "gene_set_coverage.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(sets), function(s) data.table(set = s, gene = sets[[s]]))),
       file.path(out, "gene_sets_resolved.tsv"), sep = "\t")
sc <- ScoreSignatures_UCell(counts, features = sets, maxRank = 1500, chunk.size = 200, ncores = 1, name = "")
keep <- c("Mice", "Run", "Batch", "Celltypes_Louvain", "Celltype", "Phase", "S.Score", "G2M.Score",
          "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
cells <- data.table(cell = rownames(md), md[, keep], as.data.frame(sc[rownames(md), ]))
fwrite(cells, file.path(out, "hspc_cell_scores.tsv.gz"), sep = "\t")
grp <- factor(paste(md$Mice, md$Celltypes_Louvain, sep = "|"))
pb <- counts %*% sparse.model.matrix(~ 0 + grp)
colnames(pb) <- levels(grp)
pb <- as.matrix(pb)
fwrite(data.table(gene = rownames(pb), pb), file.path(out, "pseudobulk_counts_system_by_celltype.tsv.gz"), sep = "\t")
fwrite(as.data.table(table(group = grp)), file.path(out, "pseudobulk_ncells.tsv"), sep = "\t")
cat("done\n"); print(cov); print(sessionInfo())
