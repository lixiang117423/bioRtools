# Test script for enrich_kegg/enrich_go batch mode: pass a find_degs_deseq2()-
# style data frame, auto-split by comparison, enrich per group, combine.

source("R/utils-orientation.R")
source("R/enrich_batch.R")
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

deg_df <- rbind(
  data.frame(gene = paste0("G", 1:15), comparison = "AvsC",
             regulation = c(rep("Up", 7), rep("Down", 7), "NS"), stringsAsFactors = FALSE),
  data.frame(gene = paste0("G", 6:20), comparison = "MvsC",
             regulation = c(rep("Down", 7), rep("Up", 7), "NS"), stringsAsFactors = FALSE),
  data.frame(gene = paste0("G", 1:10), comparison = "AMvsC",
             regulation = rep(c("Up", "Down"), 5), stringsAsFactors = FALSE)
)

quiet <- function(expr) suppressMessages(suppressWarnings(force(expr)))

# --- Test 1: kegg batch combines all comparisons -----------------------------
cat("Test 1: enrich_kegg batch mode\n")
res_k <- quiet(enrich_kegg(deg_df, kegg_db, p_adjust = 1.0))
stopifnot(is.data.frame(res_k), nrow(res_k) > 0)
stopifnot("comparison" %in% names(res_k))
stopifnot(setequal(unique(res_k$comparison), c("AvsC", "MvsC", "AMvsC")))
cat("✓ Test 1 passed —", nrow(res_k), "rows across",
    length(unique(res_k$comparison)), "comparisons\n\n")

# --- Test 2: batch per-comparison == single-mode (AvsC) ----------------------
cat("Test 2: batch matches single-mode per comparison\n")
avsc_genes <- unique(as.character(deg_df$gene[deg_df$comparison == "AvsC" &
                                              deg_df$regulation != "NS"]))
single <- quiet(enrich_kegg(avsc_genes, kegg_db, p_adjust = 1.0))
stopifnot(identical(sort(single$id), sort(res_k$id[res_k$comparison == "AvsC"])))
cat("✓ Test 2 passed — AvsC pathways identical to single-mode run\n\n")

# --- Test 3: keep_regulation restricts to Up only ----------------------------
cat("Test 3: keep_regulation = 'Up'\n")
res_up <- quiet(enrich_kegg(deg_df, kegg_db, p_adjust = 1.0, keep_regulation = "Up"))
stopifnot("comparison" %in% names(res_up), nrow(res_up) > 0)
up_genes_avsc <- unique(as.character(deg_df$gene[deg_df$comparison == "AvsC" &
                                                 deg_df$regulation == "Up"]))
single_up <- quiet(enrich_kegg(up_genes_avsc, kegg_db, p_adjust = 1.0))
stopifnot(identical(sort(single_up$id), sort(res_up$id[res_up$comparison == "AvsC"])))
cat("✓ Test 3 passed — Up-only genes used\n\n")

# --- Test 4: a comparison with no db overlap is skipped, not fatal -----------
cat("Test 4: errored comparison skipped gracefully\n")
deg_with_empty <- rbind(deg_df,
  data.frame(gene = paste0("X", 1:5), comparison = "Empty",
             regulation = rep("Up", 5), stringsAsFactors = FALSE))
res_e <- quiet(enrich_kegg(deg_with_empty, kegg_db, p_adjust = 1.0))
stopifnot(setequal(unique(res_e$comparison), c("AvsC", "MvsC", "AMvsC")))
cat("✓ Test 4 passed — 'Empty' comparison skipped; others intact\n\n")

# --- Test 5: go batch --------------------------------------------------------
cat("Test 5: enrich_go batch mode\n")
res_g <- quiet(enrich_go(deg_df, go_db, p_adjust = 1.0, max_gene_set = 30))
stopifnot(is.data.frame(res_g), nrow(res_g) > 0)
stopifnot("comparison" %in% names(res_g))
stopifnot(setequal(unique(res_g$comparison), c("AvsC", "MvsC", "AMvsC")))
cat("✓ Test 5 passed —", nrow(res_g), "GO rows across 3 comparisons\n\n")

# --- Test 6: vector mode unchanged (backward compat) -------------------------
cat("Test 6: vector input still works (backward compat)\n")
res_vec <- quiet(enrich_kegg(paste0("G", 1:15), kegg_db, p_adjust = 1.0))
stopifnot(is.data.frame(res_vec), nrow(res_vec) > 0)
stopifnot(!"comparison" %in% names(res_vec))
cat("✓ Test 6 passed — vector mode returns enriched df (no comparison col)\n\n")

cat("=====================================\n")
cat("All enrich batch-mode tests passed! ✓\n")
cat("=====================================\n")
