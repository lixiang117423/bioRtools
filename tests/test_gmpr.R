# Tests for internal gmpr(): GMPR size factors (Chen et al. 2018).
source("R/gmpr.R")

otu <- matrix(c(10L,20L,5L,0L,8L,
                20L,5L,10L,12L,0L,
                0L,15L,20L,8L,10L,
                5L,0L,3L,20L,15L,
                15L,10L,8L,5L,20L), nrow = 5, ncol = 5)
rownames(otu) <- paste0("ASV", 1:5)
colnames(otu) <- paste0("S", 1:5)

cat("Test 1: positive size factors, correct length + names\n")
sf <- gmpr(otu, intersect_no = 2)
stopifnot(length(sf) == 5)
stopifnot(all(sf > 0, na.rm = TRUE))
stopifnot(identical(names(sf), colnames(otu)))
stopifnot(!is.null(attr(sf, "nss")))
cat("✓ Test 1 passed\n\n")

cat("Test 2: identical samples -> identical size factors\n")
otu_id <- matrix(rep(1:5, 5), nrow = 5, ncol = 5)
rownames(otu_id) <- paste0("ASV", 1:5); colnames(otu_id) <- paste0("S", 1:5)
sf_id <- gmpr(otu_id, intersect_no = 2)
stopifnot(length(unique(round(sf_id, 8))) == 1)   # all equal
cat("✓ Test 2 passed\n\n")

cat("Test 3: scale invariance — doubling one sample raises its factor\n")
otu2 <- otu
otu2[, 1] <- otu2[, 1] * 2
sf2 <- gmpr(otu2, intersect_no = 2)
stopifnot(sf2[1] > sf[1])
cat("✓ Test 3 passed\n\n")

cat("=====================================\n")
cat("All gmpr tests passed! ✓\n")
cat("=====================================\n")
