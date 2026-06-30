# Test script for the enrich_kegg / enrich_go db column rename (snake_case)
# with a backward-compat shim that still accepts legacy dot.case names.

source("R/utils-orientation.R")
source("R/enrich_kegg.R")
source("R/enrich_go.R")

suppressPackageStartupMessages({
  library(dplyr)
})

# A small annotation db with clear overlap so enrichment returns results
gene_pool <- paste0("G", 1:20)
kegg_db <- data.frame(
  gene      = rep(gene_pool, 3),
  kegg_id   = rep(c("hsa00010", "hsa00020", "hsa00030"), each = 20),
  kegg_term = rep(c("Pathway A", "Pathway B", "Pathway C"), each = 20),
  stringsAsFactors = FALSE
)
go_db <- data.frame(
  gene       = rep(gene_pool, 3),
  go_id      = rep(c("GO:0000001", "GO:0000002", "GO:0000003"), each = 20),
  go_term    = rep(c("Term A", "Term B", "Term C"), each = 20),
  stringsAsFactors = FALSE
)
input_genes <- paste0("G", 1:15)   # 15 of 20 genes appear in every term/pathway

# Legacy dot.case versions of the same dbs
kegg_db_legacy <- setNames(kegg_db, c("gene", "kegg.id", "kegg.term"))
go_db_legacy   <- setNames(go_db,   c("gene", "go.id", "go.term"))

# --- Test 1: reconcile_db_columns helper -------------------------------------
cat("Test 1: reconcile_db_columns helper\n")
out <- reconcile_db_columns(kegg_db_legacy,
                            c(kegg_id = "kegg.id", kegg_term = "kegg.term"), "kegg_db")
stopifnot(identical(names(out), c("gene", "kegg_id", "kegg_term")))
# snake_case input is left untouched (no warning)
out2 <- reconcile_db_columns(kegg_db,
                             c(kegg_id = "kegg.id", kegg_term = "kegg.term"), "kegg_db")
stopifnot(identical(names(out2), names(kegg_db)))
cat("✓ Test 1 passed\n\n")

# --- Test 2: enrich_kegg accepts snake_case and runs -------------------------
cat("Test 2: enrich_kegg (snake_case db)\n")
res_new <- suppressWarnings(suppressMessages(
  enrich_kegg(input_genes, kegg_db, p_adjust = 1.0)))
stopifnot(is.data.frame(res_new), nrow(res_new) >= 1)
cat("✓ Test 2 passed —", nrow(res_new), "pathways returned\n\n")

# --- Test 3: enrich_kegg accepts legacy dot.case via shim (same result) -------
cat("Test 3: enrich_kegg (legacy dot.case db)\n")
warn <- NULL
res_old <- suppressMessages(withCallingHandlers(
  enrich_kegg(input_genes, kegg_db_legacy, p_adjust = 1.0),
  warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") }))
stopifnot(grepl("deprecated", warn))
stopifnot(identical(res_new$id, res_old$id))
stopifnot(isTRUE(all.equal(res_new$pvalue, res_old$pvalue)))
cat("✓ Test 3 passed — deprecation warning emitted; results identical\n\n")

# --- Test 4: enrich_go accepts snake_case and runs ---------------------------
cat("Test 4: enrich_go (snake_case db)\n")
res_go <- suppressWarnings(suppressMessages(
  enrich_go(input_genes, go_db, p_adjust = 1.0, max_gene_set = 30)))
stopifnot(is.data.frame(res_go), nrow(res_go) >= 1)
cat("✓ Test 4 passed —", nrow(res_go), "terms returned\n\n")

# --- Test 5: enrich_go accepts legacy dot.case via shim -----------------------
cat("Test 5: enrich_go (legacy dot.case db)\n")
warn2 <- NULL
res_go_old <- suppressMessages(withCallingHandlers(
  enrich_go(input_genes, go_db_legacy, p_adjust = 1.0, max_gene_set = 30),
  warning = function(w) { warn2 <<- conditionMessage(w); invokeRestart("muffleWarning") }))
stopifnot(grepl("deprecated", warn2))
stopifnot(identical(res_go$id, res_go_old$id))
cat("✓ Test 5 passed — deprecation warning; results identical\n\n")

cat("=====================================\n")
cat("All enrich column-rename tests passed! ✓\n")
cat("=====================================\n")
