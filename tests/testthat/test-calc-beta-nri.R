# Migrated from tests/test_calc_beta_nri.R (bare script).
# calc_beta_nri(): betaNRI = ses.betaMPD (MPD-based).

make_beta_nri_inputs <- function() {
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

test_that("shape + columns", {
  inputs <- make_beta_nri_inputs()
  set.seed(11)
  res <- calc_beta_nri(inputs$otu, inputs$tree, beta_reps = 50, verbose = FALSE)

  expect_true(all(c("from", "to", "beta_mpd", "beta_nri") %in% names(res)))
  expect_equal(nrow(res), 6 * 5 / 2)
  expect_false(any(is.na(res$beta_nri)))
})

test_that("parity with manual ses.betaMPD (one pair, via picante::comdist)", {
  inputs <- make_beta_nri_inputs()
  set.seed(11)
  res <- calc_beta_nri(inputs$otu, inputs$tree, beta_reps = 50, verbose = FALSE)

  comm <- t(inputs$otu)
  coph <- ape::cophenetic.phylo(inputs$tree)[colnames(comm), colnames(comm)]
  set.seed(11)
  obs_mpd <- as.matrix(picante::comdist(comm, coph, abundance.weighted = TRUE))
  nulls <- replicate(50, as.matrix(picante::comdist(
    comm, picante::taxaShuffle(coph), abundance.weighted = TRUE)))
  nri_manual <- (obs_mpd["S2", "S1"] - mean(nulls["S2", "S1", ])) /
    sd(nulls["S2", "S1", ])

  got <- res$beta_nri[res$from == "S2" & res$to == "S1"]
  expect_length(got, 1)
  expect_equal(as.numeric(got), nri_manual, tolerance = 1e-9)
})

test_that("grouped output", {
  inputs <- make_beta_nri_inputs()
  set.seed(11)
  res_g <- calc_beta_nri(inputs$otu, inputs$tree, beta_reps = 50,
                         verbose = FALSE, sample = inputs$meta,
                         group_col = "group")

  expect_true("group" %in% names(res_g))
  expect_equal(nrow(res_g), 2 * (3 * 2 / 2))
})

test_that("no tree produces an error", {
  inputs <- make_beta_nri_inputs()
  expect_error(calc_beta_nri(inputs$otu, NULL, beta_reps = 5, verbose = FALSE))
})
