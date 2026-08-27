# Migrated from tests/test_enrich_columns.R (bare script).
# enrich_kegg / enrich_go db column rename (snake_case) with a
# backward-compat shim that still accepts legacy dot.case names.

make_enrich_columns_dbs <- function() {
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
  list(kegg_db = kegg_db, go_db = go_db)
}

test_that("reconcile_db_columns helper renames legacy dot.case columns", {
  dbs <- make_enrich_columns_dbs()
  kegg_db_legacy <- setNames(dbs$kegg_db, c("gene", "kegg.id", "kegg.term"))

  # the shim signals a deprecation warning by design; asserted below, muted here
  out <- suppressWarnings(bioRtools:::reconcile_db_columns(
    kegg_db_legacy,
    c(kegg_id = "kegg.id", kegg_term = "kegg.term"), "kegg_db"))
  expect_identical(names(out), c("gene", "kegg_id", "kegg_term"))

  # snake_case input is left untouched (no warning)
  out2 <- bioRtools:::reconcile_db_columns(
    dbs$kegg_db,
    c(kegg_id = "kegg.id", kegg_term = "kegg.term"), "kegg_db")
  expect_identical(names(out2), names(dbs$kegg_db))
})

test_that("enrich_kegg accepts snake_case db and runs", {
  dbs <- make_enrich_columns_dbs()
  input_genes <- paste0("G", 1:15)
  res_new <- suppressWarnings(suppressMessages(
    enrich_kegg(input_genes, dbs$kegg_db, p_adjust = 1.0)))

  expect_s3_class(res_new, "data.frame")
  expect_gte(nrow(res_new), 1)
})

test_that("enrich_kegg accepts legacy dot.case via shim (same result)", {
  dbs <- make_enrich_columns_dbs()
  kegg_db_legacy <- setNames(dbs$kegg_db, c("gene", "kegg.id", "kegg.term"))
  input_genes <- paste0("G", 1:15)
  res_new <- suppressWarnings(suppressMessages(
    enrich_kegg(input_genes, dbs$kegg_db, p_adjust = 1.0)))
  warn <- NULL
  res_old <- suppressMessages(withCallingHandlers(
    enrich_kegg(input_genes, kegg_db_legacy, p_adjust = 1.0),
    warning = function(w) {
      warn <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }))

  expect_match(warn, "deprecated")
  expect_identical(res_new$id, res_old$id)
  expect_equal(res_new$pvalue, res_old$pvalue, tolerance = 1e-8)
})

test_that("enrich_go accepts snake_case db and runs", {
  dbs <- make_enrich_columns_dbs()
  input_genes <- paste0("G", 1:15)
  res_go <- suppressWarnings(suppressMessages(
    enrich_go(input_genes, dbs$go_db, p_adjust = 1.0, max_gene_set = 30)))

  expect_s3_class(res_go, "data.frame")
  expect_gte(nrow(res_go), 1)
})

test_that("enrich_go accepts legacy dot.case via shim", {
  dbs <- make_enrich_columns_dbs()
  go_db_legacy <- setNames(dbs$go_db, c("gene", "go.id", "go.term"))
  input_genes <- paste0("G", 1:15)
  res_go <- suppressWarnings(suppressMessages(
    enrich_go(input_genes, dbs$go_db, p_adjust = 1.0, max_gene_set = 30)))
  warn2 <- NULL
  res_go_old <- suppressMessages(withCallingHandlers(
    enrich_go(input_genes, go_db_legacy, p_adjust = 1.0, max_gene_set = 30),
    warning = function(w) {
      warn2 <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }))

  expect_match(warn2, "deprecated")
  expect_identical(res_go$id, res_go_old$id)
})
