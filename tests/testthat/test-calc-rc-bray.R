# Migrated from tests/test_calc_rc_bray.R (bare script).
# calc_rc_bray(): Raup-Crick on Bray-Curtis.

make_rc_bray_inputs <- function() {
  otu <- matrix(c(10, 20, 5, 0, 8, 15, 20, 5, 10, 12, 0, 3, 0, 15, 20, 8, 10, 5,
                  5, 0, 3, 20, 15, 10, 15, 10, 8, 5, 20, 0, 8, 12, 0, 15, 3, 20),
                nrow = 6, ncol = 6, byrow = TRUE)
  rownames(otu) <- paste0("ASV", 1:6)
  colnames(otu) <- paste0("S", 1:6)
  meta <- data.frame(sample = colnames(otu),
                     group = c("A", "A", "A", "B", "B", "B"),
                     stringsAsFactors = FALSE)
  list(otu = otu, meta = meta)
}

test_that("shape + range", {
  inputs <- make_rc_bray_inputs()
  set.seed(7)
  res <- calc_rc_bray(inputs$otu, runs = 50, verbose = FALSE)

  expect_true(all(c("from", "to", "rc_bray") %in% names(res)))
  expect_equal(nrow(res), 6 * 5 / 2)
  expect_true(all(res$rc_bray >= -1.0001 & res$rc_bray <= 1.0001))
})

test_that("reproducibility (same seed -> identical result)", {
  inputs <- make_rc_bray_inputs()
  set.seed(7)
  a <- calc_rc_bray(inputs$otu, runs = 50, verbose = FALSE)
  set.seed(7)
  b <- calc_rc_bray(inputs$otu, runs = 50, verbose = FALSE)

  expect_identical(a, b)
})

test_that("grouped output has group col + per-group pair count", {
  inputs <- make_rc_bray_inputs()
  set.seed(7)
  res_g <- calc_rc_bray(inputs$otu, runs = 50, verbose = FALSE,
                        sample = inputs$meta, group_col = "group")

  expect_true("group" %in% names(res_g))
  expect_equal(nrow(res_g), 2 * (3 * 2 / 2))
})

test_that("parity with manual picante/vegan computation (one pair)", {
  inputs <- make_rc_bray_inputs()
  set.seed(7)
  res <- calc_rc_bray(inputs$otu, runs = 50, verbose = FALSE)

  comm <- t(inputs$otu)
  set.seed(7)
  obs <- as.matrix(vegan::vegdist(comm, "bray"))
  greater <- 0
  for (r in seq_len(50)) {
    nd <- as.matrix(vegan::vegdist(
      picante::randomizeMatrix(comm, null.model = "independentswap"), "bray"))
    greater <- greater + (nd > obs)
  }
  rc_manual <- (greater / 51 - 0.5) * 2
  # pick pair (S1,S2): rc_bray fills lower triangle (i>j), so from=S2,to=S1
  s1 <- "S1"; s2 <- "S2"
  got_val <- res$rc_bray[res$from == s2 & res$to == s1]
  expect_length(got_val, 1)
  expect_equal(as.numeric(got_val), rc_manual["S2", "S1"], tolerance = 1e-9)
})

test_that("validation rejects bad input", {
  expect_error(calc_rc_bray(123))
})
