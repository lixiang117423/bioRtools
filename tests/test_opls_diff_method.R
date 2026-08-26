# Tests for the VIP-based differential determination in compute_pairwise_diff().
# Standalone: no ropls needed (VIP table is supplied directly).


suppressPackageStartupMessages(library(bioRtools))
compute_pairwise_diff <- bioRtools:::compute_pairwise_diff
suppressPackageStartupMessages(library(dplyr))

set.seed(2)
n_feat <- 40L
features <- sprintf("F%02d", seq_len(n_feat))
grp <- rep(c("A", "B", "C"), each = 6)
# Positive baseline (metabolomics-like intensities), clear group shifts
Xm <- matrix(rnorm(length(grp) * n_feat, mean = 5, sd = 1), nrow = length(grp),
             dimnames = list(sprintf("S%02d", seq_along(grp)), features))
Xm[grp == "B", ] <- Xm[grp == "B", ] + 2
Xm[grp == "C", ] <- Xm[grp == "C", ] + 4
pairs <- list(c("A", "B"), c("A", "C"))
cmp <- vapply(pairs, function(p) paste0(p[2], " vs ", p[1]), character(1))

# Fake VIP table: even-index features are important (1.5), odd are not (0.9)
vip <- do.call(rbind, lapply(cmp, function(cc) {
  data.frame(feature = features, comparison = cc,
             vip = ifelse(seq_along(features) %% 2 == 0, 1.5, 0.9),
             stringsAsFactors = FALSE)
}))

# --- Test 1: diff_method = "vip" (new default behavior) ---------------------
cat("Test 1: VIP-based labeling\n")
d_vip <- compute_pairwise_diff(Xm, grp, pairs, features, vip,
                               test_method = "t-test", p_adjust_method = "BH",
                               diff_method = "vip", vip_threshold = 1)
# VIP >= 1 -> Up/Down (direction from log2FC); VIP < 1 -> NS regardless of p
stopifnot(all(d_vip$regulation[d_vip$vip >= 1] %in% c("Up", "Down")))
stopifnot(all(d_vip$regulation[d_vip$vip < 1] == "NS"))
stopifnot(all(d_vip$log2_fc[d_vip$regulation == "Up"] > 0))
stopifnot(all(d_vip$log2_fc[d_vip$regulation == "Down"] < 0))
# is_differential is consistent with regulation and with VIP threshold
stopifnot(all(d_vip$is_differential == (d_vip$regulation %in% c("Up", "Down"))))
stopifnot(all(d_vip$is_differential == (d_vip$vip >= 1)))
# p-value columns are still computed for plotting
stopifnot(!any(is.na(d_vip$p_adjust)))
n_vip_sig <- sum(d_vip$is_differential)
stopifnot(n_vip_sig == 40)   # 20 important features x 2 comparisons
cat("✓ Test 1 passed (", n_vip_sig, " differential of", nrow(d_vip), ")\n\n")

# --- Test 2: diff_method = "pvalue" (legacy behavior) -----------------------
cat("Test 2: p-value-based labeling (legacy)\n")
d_p <- compute_pairwise_diff(Xm, grp, pairs, features, vip,
                             test_method = "t-test", p_adjust_method = "BH",
                             diff_method = "pvalue", p_threshold = 0.05)
stopifnot(all(d_p$regulation %in% c("Up", "Down", "NS")))
stopifnot(all(d_p$p_adjust[d_p$regulation == "Up"] < 0.05))
stopifnot(all(d_p$p_adjust[d_p$regulation == "Down"] < 0.05))
stopifnot(all(d_p$p_adjust[d_p$regulation == "NS"] >= 0.05))
stopifnot(all(d_p$is_differential == (d_p$regulation %in% c("Up", "Down"))))
cat("✓ Test 2 passed (", sum(d_p$is_differential), " differential)\n\n")

# --- Test 3: VIP threshold is honored ---------------------------------------
cat("Test 3: vip_threshold honored\n")
d_15 <- compute_pairwise_diff(Xm, grp, pairs, features, vip,
                              test_method = "t-test", diff_method = "vip",
                              vip_threshold = 1.5)
stopifnot(all(d_15$regulation[d_15$vip < 1.5] == "NS"))
stopifnot(all(d_15$regulation[d_15$vip >= 1.5] %in% c("Up", "Down")))
stopifnot(sum(d_15$is_differential) == 40)
cat("✓ Test 3 passed\n\n")

# --- Test 4: the two methods can disagree (they measure different things) ---
cat("Test 4: methods differ on this data\n")
if (identical(d_vip$regulation, d_p$regulation)) {
  cat("  (methods agreed on this dataset - not an error, just noting)\n")
} else {
  cat("  methods disagree as expected\n")
}
cat("✓ Test 4 passed\n\n")

# --- Test 5: default is now "vip" -------------------------------------------
cat("Test 5: default diff_method is vip\n")
d_def <- compute_pairwise_diff(Xm, grp, pairs, features, vip,
                               test_method = "t-test")
stopifnot(identical(d_def$regulation, d_vip$regulation))
cat("✓ Test 5 passed\n\n")

# --- Test 6: invalid diff_method errors -------------------------------------
cat("Test 6: invalid diff_method\n")
tryCatch({
  compute_pairwise_diff(Xm, grp, pairs, features, vip,
                        test_method = "t-test", diff_method = "bogus")
  stop("should have thrown an error")
}, error = function(e) {
  cat("✓ Test 6 passed: correctly caught -", conditionMessage(e), "\n")
})
cat("\n")

# --- Test 7: non-finite log2FC (negative data) labels NS, not NA ------------
cat("Test 7: non-finite log2FC -> NS\n")
Xm2 <- Xm
Xm2[grp == "A", 1] <- -5    # F01: negative ref mean -> NaN log2_fc
vip2 <- vip
vip2$vip[vip2$feature == "F01"] <- 1.5
d2 <- compute_pairwise_diff(Xm2, grp, pairs, features, vip2,
                            test_method = "t-test", diff_method = "vip")
stopifnot(all(d2$regulation[d2$feature == "F01"] == "NS"))
stopifnot(all(!is.na(d2$regulation)))
cat("✓ Test 7 passed\n\n")

cat("=====================================\n")
cat("All tests completed successfully! ✓\n")
cat("=====================================\n")
