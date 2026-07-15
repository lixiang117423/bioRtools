# Tests for normalize_abund(): 8 normalization methods.
source("R/gmpr.R")
source("R/normalize_abund.R")

otu <- matrix(c(10L,20L,5L,0L,8L,
                20L,5L,10L,12L,0L,
                0L,15L,20L,8L,10L,
                5L,0L,3L,20L,15L,
                15L,10L,8L,5L,20L), nrow = 5, ncol = 5)
rownames(otu) <- paste0("ASV", 1:5); colnames(otu) <- paste0("S", 1:5)

cat("Test 1: clr — exp(clr) per-column geometric mean ~ 1\n")
clr_out <- normalize_abund(otu, method = "clr", pseudocount = 0.5)
stopifnot(all(dim(clr_out) == dim(otu)))
gm_col <- exp(colMeans(clr_out))   # geometric mean of exp(clr) = exp(mean(clr))
stopifnot(isTRUE(all.equal(as.numeric(gm_col), rep(1, ncol(otu)), tolerance = 1e-8)))
cat("✓ Test 1 passed\n\n")

cat("Test 2: tss — colSums == 1\n")
tss_out <- normalize_abund(otu, method = "tss")
stopifnot(isTRUE(all.equal(as.numeric(colSums(tss_out)), rep(1, ncol(otu)), tolerance = 1e-10)))
cat("✓ Test 2 passed\n\n")

cat("Test 3: hellinger — values in [0,1]\n")
hel_out <- normalize_abund(otu, method = "hellinger")
stopifnot(min(hel_out) >= 0, max(hel_out, na.rm = TRUE) <= 1.0000001)
cat("✓ Test 3 passed\n\n")

cat("Test 4: tmm + rle — positive size_factor attr\n")
for (mm in c("tmm", "rle")) {
  out <- normalize_abund(otu, method = mm)
  stopifnot(all(dim(out) == dim(otu)))
  sf <- attr(out, "size_factor")
  stopifnot(!is.null(sf), length(sf) == ncol(otu), all(sf > 0))
}
cat("✓ Test 4 passed\n\n")

cat("Test 5: gmpr — dispatches + attaches size_factor attr\n")
# 5-taxon fixture < default intersect_no=10 -> gmpr warns + some NAs (expected;
# real data has 1000s of taxa). gmpr itself is tested rigorously in test_gmpr.R.
out <- suppressWarnings(normalize_abund(otu, method = "gmpr"))
stopifnot(all(dim(out) == dim(otu)))
sf <- attr(out, "size_factor")
stopifnot(!is.null(sf), length(sf) == ncol(otu))
stopifnot(all(sf[!is.na(sf)] > 0))
cat("✓ Test 5 passed\n\n")

cat("Test 6: rarefy — all colSums == depth; reproducible\n")
depth <- min(colSums(otu))
# vegan::rrarefy prints a cautionary warning on this tiny fixture; suppressed.
rar_out <- suppressWarnings(normalize_abund(otu, method = "rarefy", depth = depth, seed = 7))
stopifnot(isTRUE(all.equal(as.numeric(colSums(rar_out)), rep(depth, ncol(otu)))))
rar_out2 <- suppressWarnings(normalize_abund(otu, method = "rarefy", depth = depth, seed = 7))
stopifnot(identical(rar_out, rar_out2))
cat("✓ Test 6 passed\n\n")

cat("Test 7: css — conditional (skip if metagenomeSeq absent)\n")
if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
  cat("  (metagenomeSeq not installed — css skipped)\n")
} else {
  out <- normalize_abund(otu, method = "css")
  stopifnot(all(dim(out) == dim(otu)))
  stopifnot(!is.null(attr(out, "size_factor")))
  cat("✓ css passed\n")
}
cat("\n")

cat("Test 8: errors\n")
stopifnot(inherits(try(normalize_abund(otu, method = "nope"), silent = TRUE), "try-error"))
stopifnot(inherits(try(normalize_abund(otu + 0.5, method = "rarefy"), silent = TRUE), "try-error"))  # non-integer
cat("✓ Test 8 passed\n\n")

cat("=====================================\n")
cat("All normalize_abund tests passed! ✓\n")
cat("=====================================\n")
