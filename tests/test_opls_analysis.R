# Tests for opls_analysis() modes + internal pairwise helpers.
source("R/opls_pairwise_fit.R")
source("R/utils-orientation.R")
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
# summaries must carry real model-quality values (ropls cols are R2Y(cum)/Q2(cum))
smry <- fit$summaries[[1]]
stopifnot(!is.na(smry$R2Y), !is.na(smry$Q2Y))
cat("✓ Test A1 passed\n\n")

# --- Test A2: compute_pairwise_diff structure + per-comparison p-adjust -----
cat("Test A2: compute_pairwise_diff computes log2FC/p-value/regulation\n")
diff <- compute_pairwise_diff(Xm, grp, pairs, features, fit$vip,
                              test_method = "t-test", p_adjust_method = "BH",
                              p_threshold = 0.05)

expected_cols <- c("feature", "group", "ref_group", "log2_fc", "p_value",
                   "p_adjust", "comparison", "vip", "regulation")
stopifnot(all(expected_cols %in% names(diff)))
stopifnot(nrow(diff) == length(pairs) * n_feat)          # one row per pair x variable
stopifnot(length(unique(diff$comparison)) == 3)
stopifnot(all(diff$regulation %in% c("Up", "Down", "NS")))
# Direction sanity: group B has higher mean than A (shift +2*3 in construction)
b_vs_a <- diff[diff$comparison == "B vs A", ]
stopifnot(all(b_vs_a$log2_fc > 0))                       # B up vs A
cat("✓ Test A2 passed\n\n")

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# --- Test A3a: pairwise = TRUE -> all-pairs -------------------------------
cat("Test A3a: pairwise = TRUE produces all-pairs comparisons\n")
sample_meta <- tibble::tibble(
  sample_id = rownames(Xm), group = grp
)
res_ap <- quiet(opls_analysis(Xm, sample = sample_meta, sample_col = "sample_id",
                              group_col = "group", pairwise = TRUE, verbose = FALSE))
stopifnot(length(unique(res_ap$vip_scores$comparison)) == 3)   # choose(3,2) = 3
stopifnot(!is.null(res_ap$models) && is.list(res_ap$models))
stopifnot(length(res_ap$models) == 3)
stopifnot(!is.null(res_ap$differential_analysis))
stopifnot(length(unique(res_ap$differential_analysis$comparison)) == 3)
cat("✓ Test A3a passed\n\n")

# --- Test A3b: defaults unchanged (single binary model) --------------------
# ropls OPLS-DA is binary-only, so single-model mode is exercised on 2 groups.
cat("Test A3b: default (2 groups, no pairwise) is a single binary model\n")
two_grp <- grp %in% c("A", "B")
Xm2 <- Xm[two_grp, ]
meta2 <- tibble::tibble(sample_id = rownames(Xm2), group = grp[two_grp])
res_single <- quiet(opls_analysis(Xm2, sample = meta2, sample_col = "sample_id",
                                  group_col = "group", verbose = FALSE))
stopifnot(is.null(res_single$models))                          # not pairwise
stopifnot(is.null(res_single$differential_analysis))           # no ref_group
stopifnot(!"comparison" %in% names(res_single$vip_scores))
cat("✓ Test A3b passed\n\n")

# --- Test A3c: ref_group ref-anchored unchanged ---------------------------
cat("Test A3c: ref_group produces ref-anchored pairwise + diff\n")
res_ref <- quiet(opls_analysis(Xm, sample = sample_meta, sample_col = "sample_id",
                               group_col = "group", ref_group = "A", verbose = FALSE))
stopifnot(length(unique(res_ref$vip_scores$comparison)) == 2)  # B vs A, C vs A
stopifnot(all(c("B vs A", "C vs A") %in% unique(res_ref$vip_scores$comparison)))
stopifnot(!is.null(res_ref$differential_analysis))
stopifnot(all(unique(res_ref$differential_analysis$ref_group) == "A"))
cat("✓ Test A3c passed\n\n")

cat("=====================================\n")
cat("All opls_analysis tests passed! ✓\n")
cat("=====================================\n")
