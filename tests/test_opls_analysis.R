# Tests for opls_analysis() modes + internal pairwise helpers.
source("R/opls_pairwise_fit.R")
source("R/opls_analysis.R")

suppressPackageStartupMessages({
  library(ropls); library(tibble); library(dplyr)
})

set.seed(1)
n_feat <- 30L
features <- sprintf("F%02d", seq_len(n_feat))
# Three groups of 6 samples each, with a clear group signal.
grp <- rep(c("A", "B", "C"), each = 6)
Xm <- matrix(rnorm(length(grp) * n_feat), nrow = length(grp),
             dimnames = list(sprintf("S%02d", seq_along(grp)), features))
for (g in seq_along(unique(grp))) Xm[grp == unique(grp)[g], ] <-
  Xm[grp == unique(grp)[g], ] + g * 3

# --- Test A1: fit_pairwise_opls structure ---------------------------------
cat("Test A1: fit_pairwise_opls returns correct structure\n")
pairs <- utils::combn(unique(grp), 2, simplify = FALSE)  # A-B, A-C, B-C
fit <- fit_pairwise_opls(Xm, grp, pairs, features,
                         ortho_components = 1, scaling = "standard",
                         validation = "CV", cv_folds = 5, verbose = FALSE)

stopifnot(length(fit$models) == 3)                       # 3 comparisons
stopifnot(all(c("B vs A", "C vs A", "C vs B") %in% names(fit$models)))
# Each comparison keeps both groups' samples (6+6 = 12); 3 comparisons -> 36 rows.
stopifnot(nrow(fit$scores) == 36)
stopifnot("comparison" %in% names(fit$scores))
stopifnot("vip" %in% names(fit$vip) && "comparison" %in% names(fit$vip))
stopifnot("group" %in% names(fit$vip) && "ref_group" %in% names(fit$vip))
cat("✓ Test A1 passed\n\n")
