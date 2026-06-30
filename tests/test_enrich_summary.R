# Regression test for the enrich_kegg/enrich_go interactive-summary crash.
# The summary referenced pre-rename output columns (Description/ID/p.adjust/...)
# that no longer exist after the rename to snake_case, so $Description was NULL
# and nchar(NULL) -> "argument is of length zero". This verifies the renamed
# columns exist and the summary loop runs on real enriched results.

source("R/utils-orientation.R")
source("R/enrich_kegg.R")
source("R/enrich_go.R")

suppressPackageStartupMessages(library(dplyr))

gene_pool <- paste0("G", 1:20)
kegg_db <- data.frame(
  gene = rep(gene_pool, 3),
  kegg_id = rep(c("hsa00010", "hsa00020", "hsa00030"), each = 20),
  kegg_term = rep(c("Pathway A", "Pathway B", "Pathway C"), each = 20),
  stringsAsFactors = FALSE
)
go_db <- data.frame(
  gene = rep(gene_pool, 3),
  go_id = rep(c("GO:0000001", "GO:0000002", "GO:0000003"), each = 20),
  go_term = rep(c("Term A", "Term B", "Term C"), each = 20),
  stringsAsFactors = FALSE
)
input_genes <- paste0("G", 1:15)

res_k <- suppressWarnings(suppressMessages(
  enrich_kegg(input_genes, kegg_db, p_adjust = 1.0)))
res_g <- suppressWarnings(suppressMessages(
  enrich_go(input_genes, go_db, p_adjust = 1.0, max_gene_set = 30)))

stopifnot(nrow(res_k) >= 1, nrow(res_g) >= 1)

# Columns the interactive summary now references (must all resolve)
need_k <- c("id", "description", "p_adjust", "gene_count", "total_genes", "enrichment_score")
need_g <- c("id", "description", "p_adjust", "gene_count")
stopifnot(all(need_k %in% names(res_k)))
stopifnot(all(need_g %in% names(res_g)))
cat("✓ renamed summary columns present in both results\n")

# Run the exact summary loop on real results (the path that crashed)
out_k <- character(0)
top <- head(res_k, 5)
for (i in seq_len(nrow(top))) {
  nm <- top$description[i]
  if (!is.na(nm) && nchar(nm) > 50) nm <- paste0(substr(nm, 1, 47), "...")
  out_k <- c(out_k, sprintf("%d. %s: %s | p.adj=%.2e genes=%d/%d enrich=%.1f",
    i, top$id[i], nm, top$p_adjust[i], top$gene_count[i],
    top$total_genes[i], top$enrichment_score[i]))
}
cat("kegg summary line 1:", out_k[1], "\n")

out_g <- character(0)
topg <- head(res_g, 5)
for (i in seq_len(nrow(topg))) {
  out_g <- c(out_g, sprintf("%s: %s (p.adj=%.2e, genes=%d)",
    topg$id[i], substr(topg$description[i], 1, 50), topg$p_adjust[i], topg$gene_count[i]))
}
cat("go summary line 1  :", out_g[1], "\n")
stopifnot(length(out_k) == min(5, nrow(res_k)), length(out_g) == min(5, nrow(res_g)))

cat("\n✓ summary loops run without error (regression fixed)\n")
cat("=====================================\n")
cat("All enrich-summary regression tests passed! ✓\n")
cat("=====================================\n")
