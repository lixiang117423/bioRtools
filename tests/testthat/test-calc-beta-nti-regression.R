# Migrated from tests/test_calc_beta_nti_regression.R (bare script).
# Regression: calc_beta_nti output must be numerically equivalent after the
# calc_ses_beta_core refactor.

make_beta_nti_inputs <- function() {
  otu <- matrix(c(10, 20, 5, 0, 8, 15, 20, 5, 10, 12, 0, 3, 0, 15, 20, 8, 10, 5,
                  5, 0, 3, 20, 15, 10, 15, 10, 8, 5, 20, 0, 8, 12, 0, 15, 3, 20),
                nrow = 6, ncol = 6, byrow = TRUE)
  rownames(otu) <- paste0("ASV", 1:6)
  colnames(otu) <- paste0("S", 1:6)
  tree <- ape::read.tree(text =
    "(((ASV1:0.1,ASV2:0.2):0.1,ASV3:0.15):0.2,((ASV4:0.1,ASV5:0.2):0.15,ASV6:0.2):0.1);")
  meta <- data.frame(sample = colnames(otu),
                     group = c("A", "A", "A", "B", "B", "B"),
                     stringsAsFactors = FALSE)
  list(otu = otu, tree = tree, meta = meta)
}

baseline <- readRDS(system.file("extdata", "calc_beta_nti_fixture.rds",
                                package = "bioRtools"))

test_that("ungrouped output matches baseline", {
  inputs <- make_beta_nti_inputs()
  set.seed(42)
  got <- calc_beta_nti(inputs$otu, inputs$tree, beta_reps = 50, verbose = FALSE)

  expect_equal(got, baseline$all, tolerance = 1e-12)
})

test_that("grouped output matches baseline", {
  inputs <- make_beta_nti_inputs()
  set.seed(42)
  got_g <- calc_beta_nti(inputs$otu, inputs$tree, beta_reps = 50,
                         verbose = FALSE, sample = inputs$meta,
                         group_col = "group")

  expect_equal(got_g, baseline$grouped, tolerance = 1e-12)
})

test_that("tree/data mismatch errors", {
  inputs <- make_beta_nti_inputs()
  # Tree with zero tip overlap -> triggers "No matching taxa" stop().
  # (A partial mismatch does not error in either old or new code, since
  # intersect() keeps the shared taxa, so use fully disjoint labels here.)
  bad_tree <- ape::read.tree(text = "((X1:0.1,X2:0.2):0.1,X3:0.15);")
  expect_error(calc_beta_nti(inputs$otu, bad_tree, beta_reps = 5,
                             verbose = FALSE))
})
