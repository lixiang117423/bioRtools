# Migrated from tests/test_normalize_abund.R (bare script).
# normalize_abund(): 8 normalization methods.

otu <- matrix(c(10L, 20L, 5L, 0L, 8L,
                20L, 5L, 10L, 12L, 0L,
                0L, 15L, 20L, 8L, 10L,
                5L, 0L, 3L, 20L, 15L,
                15L, 10L, 8L, 5L, 20L), nrow = 5, ncol = 5)
rownames(otu) <- paste0("ASV", 1:5)
colnames(otu) <- paste0("S", 1:5)

test_that("clr — exp(clr) per-column geometric mean ~ 1", {
  clr_out <- normalize_abund(otu, method = "clr", pseudocount = 0.5)

  expect_true(all(dim(clr_out) == dim(otu)))
  gm_col <- exp(colMeans(clr_out))   # geometric mean of exp(clr) = exp(mean(clr))
  expect_equal(as.numeric(gm_col), rep(1, ncol(otu)), tolerance = 1e-8)
})

test_that("tss — colSums == 1", {
  tss_out <- normalize_abund(otu, method = "tss")

  expect_equal(as.numeric(colSums(tss_out)), rep(1, ncol(otu)),
               tolerance = 1e-10)
})

test_that("hellinger — values in [0,1]", {
  hel_out <- normalize_abund(otu, method = "hellinger")

  expect_true(min(hel_out) >= 0)
  expect_true(max(hel_out, na.rm = TRUE) <= 1.0000001)
})

test_that("tmm + rle — positive size_factor attr", {
  for (mm in c("tmm", "rle")) {
    out <- normalize_abund(otu, method = mm)
    expect_true(all(dim(out) == dim(otu)))
    sf <- attr(out, "size_factor")
    expect_false(is.null(sf))
    expect_length(sf, ncol(otu))
    expect_true(all(sf > 0))
  }
})

test_that("gmpr — dispatches + attaches size_factor attr", {
  # 5-taxon fixture < default intersect_no=10 -> gmpr warns + some NAs
  # (expected; real data has 1000s of taxa). gmpr itself is tested
  # rigorously in test-gmpr.R.
  out <- suppressWarnings(normalize_abund(otu, method = "gmpr"))

  expect_true(all(dim(out) == dim(otu)))
  sf <- attr(out, "size_factor")
  expect_false(is.null(sf))
  expect_length(sf, ncol(otu))
  expect_true(all(sf[!is.na(sf)] > 0))
})

test_that("rarefy — all colSums == depth; reproducible", {
  depth <- min(colSums(otu))
  # vegan::rrarefy prints a cautionary warning on this tiny fixture;
  # suppressed.
  rar_out <- suppressWarnings(normalize_abund(otu, method = "rarefy",
                                              depth = depth, seed = 7))
  expect_equal(as.numeric(colSums(rar_out)), rep(depth, ncol(otu)),
               tolerance = 1e-8)

  rar_out2 <- suppressWarnings(normalize_abund(otu, method = "rarefy",
                                               depth = depth, seed = 7))
  expect_identical(rar_out, rar_out2)
})

test_that("css — conditional on metagenomeSeq", {
  skip_if_not_installed("metagenomeSeq")

  out <- normalize_abund(otu, method = "css")
  expect_true(all(dim(out) == dim(otu)))
  expect_false(is.null(attr(out, "size_factor")))
})

test_that("errors", {
  expect_error(normalize_abund(otu, method = "nope"))
  expect_error(normalize_abund(otu + 0.5, method = "rarefy"))   # non-integer
})
