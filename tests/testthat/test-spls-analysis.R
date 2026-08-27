# Migrated from tests/test_spls_analysis.R (bare script).
# spls_analysis() orientation handling: a feature x sample matrix (with or
# without feature_as_row = TRUE) and the auto-detect path must reproduce the
# sample x feature result exactly.

detect_feature_as_row <- bioRtools:::detect_feature_as_row

set.seed(1)
n_sample <- 12L
n_feat <- 30L
samples <- sprintf("S%02d", seq_len(n_sample))
features <- sprintf("F%02d", seq_len(n_feat))

# Sample x feature abundance matrix (the historical layout)
X <- matrix(stats::rnorm(n_sample * n_feat), nrow = n_sample,
            dimnames = list(samples, features))

sample_meta <- tibble::tibble(
  sample_id = samples,
  group = rep(c("A", "B", "C"), each = 4),
  batch = rep(c("b1", "b2"), times = 6)
)

run <- function(data, ...) {
  suppressMessages(suppressWarnings(
    spls_analysis(data, sample = sample_meta, sample_col = "sample_id",
                  group_col = "group", ncomp = 2, verbose = FALSE, ...)
  ))
}

# Reference: samples as rows
ref <- run(X)

score_mat <- function(res) {
  res$sample_scores %>%
    dplyr::arrange(sample) %>%
    dplyr::select(dplyr::starts_with("comp")) %>%
    as.matrix()
}
ref_mat <- score_mat(ref)

# Features-as-rows layout used by the explicit and auto-detect paths
Xf <- t(X)

test_that("detect_feature_as_row identifies the matrix orientation", {
  expect_false(detect_feature_as_row(X, sample_meta, "sample_id"))   # sample-row
  expect_true(detect_feature_as_row(t(X), sample_meta, "sample_id")) # feature-row
  expect_false(detect_feature_as_row(X, NULL, "sample_id"))          # no metadata
})

test_that("feature_as_row = TRUE reproduces the sample-row result", {
  res_force <- run(Xf, feature_as_row = TRUE)

  expect_identical(dim(score_mat(res_force)), dim(ref_mat))
  expect_equal(score_mat(res_force), ref_mat, tolerance = 1e-10)
})

test_that("auto-detect (feature_as_row = NA) reproduces the reference", {
  res_auto <- run(Xf)   # no feature_as_row -> auto-detect must catch feature-row

  expect_equal(score_mat(res_auto), ref_mat, tolerance = 1e-10)
})

test_that("backward compatibility: sample-row input + default NA", {
  res_back <- run(X, feature_as_row = NA)

  expect_equal(score_mat(res_back), ref_mat, tolerance = 1e-10)
})
