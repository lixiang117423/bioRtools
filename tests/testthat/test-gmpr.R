# Migrated from tests/test_gmpr.R (bare script).
# Internal gmpr(): GMPR size factors (Chen et al. 2018).

make_gmpr_otu <- function() {
  otu <- matrix(c(10L, 20L, 5L, 0L, 8L,
                  20L, 5L, 10L, 12L, 0L,
                  0L, 15L, 20L, 8L, 10L,
                  5L, 0L, 3L, 20L, 15L,
                  15L, 10L, 8L, 5L, 20L), nrow = 5, ncol = 5)
  rownames(otu) <- paste0("ASV", 1:5)
  colnames(otu) <- paste0("S", 1:5)
  otu
}

test_that("positive size factors, correct length + names", {
  otu <- make_gmpr_otu()
  sf <- bioRtools:::gmpr(otu, intersect_no = 2)

  expect_length(sf, 5)
  expect_true(all(sf > 0, na.rm = TRUE))
  expect_identical(names(sf), colnames(otu))
  expect_false(is.null(attr(sf, "nss")))
})

test_that("identical samples -> identical size factors", {
  otu_id <- matrix(rep(1:5, 5), nrow = 5, ncol = 5)
  rownames(otu_id) <- paste0("ASV", 1:5)
  colnames(otu_id) <- paste0("S", 1:5)
  sf_id <- bioRtools:::gmpr(otu_id, intersect_no = 2)

  expect_equal(length(unique(round(sf_id, 8))), 1)   # all equal
})

test_that("scale invariance — doubling one sample raises its factor", {
  otu <- make_gmpr_otu()
  sf <- bioRtools:::gmpr(otu, intersect_no = 2)

  otu2 <- otu
  otu2[, 1] <- otu2[, 1] * 2
  sf2 <- bioRtools:::gmpr(otu2, intersect_no = 2)

  expect_gt(sf2[1], sf[1])
})
