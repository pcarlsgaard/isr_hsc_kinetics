#!/usr/bin/env Rscript
# Structural audit of one .rds object (no analysis). Usage:
#   Rscript scripts/audit_seurat_object.R <path.rds> <out.txt>
# Run one object per process to bound memory.
suppressPackageStartupMessages({library(Matrix); library(SeuratObject)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
sink(args[2], split = TRUE)
cat("file:", args[1], "\nsize_bytes:", file.size(args[1]), "\n")
obj <- readRDS(args[1])
cat("class:", paste(class(obj), collapse = ", "), "\n")
if (identical(class(obj)[1], "Seurat")) {
  cat("object version:", as.character(obj@version), "\n")
  # slot access only: objects written by Seurat v3 fail S4 dispatch without the Seurat package
  cat("active assay:", obj@active.assay, " n cells (meta.data):", nrow(obj@meta.data), "\n")
  for (a in names(obj@assays)) {
    A <- obj@assays[[a]]
    cat("\n## assay", a, "class", class(A)[1], "\n")
    slots <- if (identical(class(A)[1], "Assay5")) names(A@layers) else c("counts", "data", "scale.data")
    for (s in slots) {
      m <- if (identical(class(A)[1], "Assay5")) A@layers[[s]] else methods::slot(A, s)
      cat(sprintf("  %-12s class=%s dim=%s", s, class(m)[1], paste(dim(m), collapse = "x")))
      if (length(m) > 0 && inherits(m, "sparseMatrix")) {
        x <- m@x
        cat(sprintf(" nnz=%d integer_valued=%s min=%g max=%g", length(x),
                    all(x[seq_len(min(length(x), 2e6))] %% 1 == 0), min(x), max(x)))
      }
      cat("\n")
    }
    rn <- if (identical(class(A)[1], "Assay5")) rownames(A@features) else rownames(A@counts)
    cat("  features (first 15):", head(rn, 15), "\n")
    cat("  n ENSMUSG-like:", sum(grepl("^ENSMUSG", rn)),
        " n mt-:", sum(grepl("^mt-", rn)),
        " n Rps/Rpl:", sum(grepl("^Rp[sl][0-9]", rn)), "\n")
  }
  cat("\nreductions:", names(obj@reductions), "\n")
  md <- obj@meta.data
  cat("\n## meta.data columns\n")
  for (cn in colnames(md)) {
    v <- md[[cn]]
    if (is.numeric(v)) {
      cat(sprintf("%-28s numeric  range=[%g, %g] NA=%d\n", cn, min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(is.na(v))))
    } else {
      u <- unique(as.character(v))
      cat(sprintf("%-28s %-8s n_unique=%d: %s\n", cn, class(v)[1], length(u), paste(head(u, 12), collapse = " | ")))
    }
  }
  cat("\n## cross-tabulations of low-cardinality columns against sample-like columns\n")
  lowc <- names(Filter(function(v) !is.numeric(v) && length(unique(v)) <= 40, md))
  samp <- intersect(c("Mice", "Mouse", "orig.ident", "Sample", "sample"), lowc)
  for (s in samp) for (cn in setdiff(lowc, s)) {
    cat("\n#", cn, "x", s, "\n"); print(table(md[[cn]], md[[s]], useNA = "ifany"))
  }
  if (length(obj@misc)) cat("\nmisc names:", names(obj@misc), "\n")
  if (length(obj@commands)) cat("commands:", names(obj@commands), "\n")
} else {
  cat("slots/names:\n"); print(utils::str(obj, max.level = 2, list.len = 40))
}
cat("\n"); print(sessionInfo())
sink()
