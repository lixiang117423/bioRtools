# Migrated from tests/test_enrich_batch.R (bare script).
# enrich_kegg/enrich_go batch mode: pass a find_degs_deseq2()-style data
# frame, auto-split by comparison, enrich per group, combine.

make_enrich_batch_dbs <- function() {
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
  list(kegg_db = kegg_db, go_db = go_db)
}

make_deg_df <- function() {
  rbind(
    data.frame(gene = paste0("G", 1:15), comparison = "AvsC",
               regulation = c(rep("Up", 7), rep("Down", 7), "NS"),
               stringsAsFactors = FALSE),
    data.frame(gene = paste0("G", 6:20), comparison = "MvsC",
               regulation = c(rep("Down", 7), rep("Up", 7), "NS"),
               stringsAsFactors = FALSE),
    data.frame(gene = paste0("G", 1:10), comparison = "AMvsC",
               regulation = rep(c("Up", "Down"), 5), stringsAsFactors = FALSE)
  )
}

quiet <- function(expr) suppressMessages(suppressWarnings(force(expr)))

test_that("enrich_kegg batch mode combines all comparisons", {
  dbs <- make_enrich_batch_dbs()
  res_k <- quiet(enrich_kegg(make_deg_df(), dbs$kegg_db, p_adjust = 1.0))

  expect_s3_class(res_k, "data.frame")
  expect_gt(nrow(res_k), 0)
  expect_true("comparison" %in% names(res_k))
  expect_setequal(unique(res_k$comparison), c("AvsC", "MvsC", "AMvsC"))
})

test_that("batch matches single-mode per comparison", {
  dbs <- make_enrich_batch_dbs()
  deg_df <- make_deg_df()
  res_k <- quiet(enrich_kegg(deg_df, dbs$kegg_db, p_adjust = 1.0))

  avsc_genes <- unique(as.character(deg_df$gene[deg_df$comparison == "AvsC" &
                                                  deg_df$regulation != "NS"]))
  single <- quiet(enrich_kegg(avsc_genes, dbs$kegg_db, p_adjust = 1.0))
  expect_identical(sort(single$id), sort(res_k$id[res_k$comparison == "AvsC"]))
})

test_that("keep_regulation = 'Up' restricts to Up only", {
  dbs <- make_enrich_batch_dbs()
  deg_df <- make_deg_df()
  res_up <- quiet(enrich_kegg(deg_df, dbs$kegg_db, p_adjust = 1.0,
                              keep_regulation = "Up"))

  expect_true("comparison" %in% names(res_up))
  expect_gt(nrow(res_up), 0)
  up_genes_avsc <- unique(as.character(deg_df$gene[deg_df$comparison == "AvsC" &
                                                     deg_df$regulation == "Up"]))
  single_up <- quiet(enrich_kegg(up_genes_avsc, dbs$kegg_db, p_adjust = 1.0))
  expect_identical(sort(single_up$id),
                   sort(res_up$id[res_up$comparison == "AvsC"]))
})

test_that("a comparison with no db overlap is skipped, not fatal", {
  dbs <- make_enrich_batch_dbs()
  deg_with_empty <- rbind(make_deg_df(),
    data.frame(gene = paste0("X", 1:5), comparison = "Empty",
               regulation = rep("Up", 5), stringsAsFactors = FALSE))
  res_e <- quiet(enrich_kegg(deg_with_empty, dbs$kegg_db, p_adjust = 1.0))

  expect_setequal(unique(res_e$comparison), c("AvsC", "MvsC", "AMvsC"))
})

test_that("enrich_go batch mode", {
  dbs <- make_enrich_batch_dbs()
  res_g <- quiet(enrich_go(make_deg_df(), dbs$go_db, p_adjust = 1.0,
                           max_gene_set = 30))

  expect_s3_class(res_g, "data.frame")
  expect_gt(nrow(res_g), 0)
  expect_true("comparison" %in% names(res_g))
  expect_setequal(unique(res_g$comparison), c("AvsC", "MvsC", "AMvsC"))
})

test_that("vector input still works (backward compat)", {
  dbs <- make_enrich_batch_dbs()
  res_vec <- quiet(enrich_kegg(paste0("G", 1:15), dbs$kegg_db, p_adjust = 1.0))

  expect_s3_class(res_vec, "data.frame")
  expect_gt(nrow(res_vec), 0)
  expect_false("comparison" %in% names(res_vec))
})
