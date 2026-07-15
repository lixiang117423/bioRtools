# Tests for calc_rc_bray(): Raup-Crick on Bray-Curtis.
source("R/calc_rc_bray.R")
suppressPackageStartupMessages({ library(picante); library(vegan); library(magrittr) })

otu <- matrix(c(10,20,5,0,8,15, 20,5,10,12,0,3, 0,15,20,8,10,5,
                5,0,3,20,15,10, 15,10,8,5,20,0, 8,12,0,15,3,20),
              nrow = 6, ncol = 6, byrow = TRUE)
rownames(otu) <- paste0("ASV", 1:6)
colnames(otu) <- paste0("S", 1:6)
meta <- data.frame(sample = colnames(otu),
                   group = c("A","A","A","B","B","B"), stringsAsFactors = FALSE)

cat("Test 1: shape + range\n")
set.seed(7)
res <- calc_rc_bray(otu, runs = 50, verbose = FALSE)
stopifnot(all(c("from","to","rc_bray") %in% names(res)))
stopifnot(nrow(res) == 6 * 5 / 2)              # 15 pairs
stopifnot(all(res$rc_bray >= -1.0001 & res$rc_bray <= 1.0001))
cat("✓ Test 1 passed\n\n")

cat("Test 2: reproducibility (same seed -> identical result)\n")
set.seed(7); a <- calc_rc_bray(otu, runs = 50, verbose = FALSE)
set.seed(7); b <- calc_rc_bray(otu, runs = 50, verbose = FALSE)
stopifnot(identical(a, b))
cat("✓ Test 2 passed\n\n")

cat("Test 3: grouped output has group col + per-group pair count\n")
set.seed(7)
res_g <- calc_rc_bray(otu, runs = 50, verbose = FALSE, sample = meta, group_col = "group")
stopifnot("group" %in% names(res_g))
stopifnot(nrow(res_g) == 2 * (3 * 2 / 2))     # 2 groups x 3 pairs
cat("✓ Test 3 passed\n\n")

cat("Test 4: parity with manual picante/vegan computation (one pair)\n")
comm <- t(otu)
set.seed(7)
obs <- as.matrix(vegan::vegdist(comm, "bray"))
greater <- 0
for (r in seq_len(50)) {
  nd <- as.matrix(vegan::vegdist(picante::randomizeMatrix(comm, null.model = "independentswap"), "bray"))
  greater <- greater + (nd > obs)
}
rc_manual <- (greater / 51 - 0.5) * 2
# pick pair (S1,S2): rc_bray fills lower triangle (i>j), so from=S2,to=S1
s1 <- "S1"; s2 <- "S2"
got_val <- res$rc_bray[res$from == s2 & res$to == s1]
stopifnot(length(got_val) == 1)
stopifnot(isTRUE(all.equal(as.numeric(got_val), rc_manual["S2","S1"], tolerance = 1e-9)))
cat("✓ Test 4 passed\n\n")

cat("Test 5: validation\n")
stopifnot(inherits(try(calc_rc_bray(123), silent = TRUE), "try-error"))
cat("✓ Test 5 passed\n\n")

cat("=====================================\n")
cat("All calc_rc_bray tests passed! ✓\n")
cat("=====================================\n")
