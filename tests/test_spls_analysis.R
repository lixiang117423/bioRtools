# Test script for spls_analysis() orientation handling.
# A feature x sample matrix (with or without feature_as_row = TRUE) and the
# auto-detect path must reproduce the sample x feature result exactly.

source("R/utils-orientation.R")
source("R/spls_analysis.R")

suppressPackageStartupMessages({
  library(mixOmics)
  library(tibble)
  library(dplyr)
})

set.seed(1)
n_sample <- 12L
n_feat <- 30L
samples <- sprintf("S%02d", seq_len(n_sample))
features <- sprintf("F%02d", seq_len(n_feat))

# Sample x feature abundance matrix (the historical layout)
X <- matrix(rnorm(n_sample * n_feat), nrow = n_sample,
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
ref_scores <- ref$sample_scores %>%
  dplyr::arrange(sample) %>%
  dplyr::select(sample, dplyr::starts_with("comp"))

score_mat <- function(res) {
  res$sample_scores %>%
    dplyr::arrange(sample) %>%
    dplyr::select(dplyr::starts_with("comp")) %>%
    as.matrix()
}
ref_mat <- score_mat(ref)

# --- Test 1: detect_feature_as_row direct checks ----------------------------
cat("Test 1: detect_feature_as_row direct checks\n")
stopifnot(detect_feature_as_row(X, sample_meta, "sample_id") == FALSE)   # sample-row
stopifnot(detect_feature_as_row(t(X), sample_meta, "sample_id") == TRUE) # feature-row
stopifnot(detect_feature_as_row(X, NULL, "sample_id") == FALSE)          # no metadata
cat("✓ Test 1 passed\n\n")

# --- Test 2: feature_as_row = TRUE reproduces sample-row result -------------
cat("Test 2: explicit feature_as_row = TRUE reproduces reference\n")
Xf <- t(X)   # features as rows, samples as columns
res_force <- run(Xf, feature_as_row = TRUE)
stopifnot(identical(dim(score_mat(res_force)), dim(ref_mat)))
stopifnot(all.equal(score_mat(res_force), ref_mat, tolerance = 1e-10))
cat("✓ Test 2 passed — scores identical to sample-row input\n\n")

# --- Test 3: auto-detect (default NA) reproduces reference ------------------
cat("Test 3: auto-detect (feature_as_row = NA) reproduces reference\n")
res_auto <- run(Xf)   # no feature_as_row -> auto-detect must catch feature-row
stopifnot(all.equal(score_mat(res_auto), ref_mat, tolerance = 1e-10))
cat("✓ Test 3 passed — auto-detect transposed correctly\n\n")

# --- Test 4: backward compatibility — sample-row still works with NA ---------
cat("Test 4: backward compatibility, sample-row input + default NA\n")
res_back <- run(X, feature_as_row = NA)
stopifnot(all.equal(score_mat(res_back), ref_mat, tolerance = 1e-10))
cat("✓ Test 4 passed — no spurious transpose on sample-row data\n\n")

cat("=====================================\n")
cat("All spls_analysis orientation tests passed! ✓\n")
cat("=====================================\n")
