# Tests for calc_beta_nri(): betaNRI = ses.betaMPD (MPD-based).
source("R/calc_ses_beta_core.R")
source("R/calc_beta_nri.R")
suppressPackageStartupMessages({ library(ape); library(picante); library(magrittr) })

otu <- matrix(c(10,20,5,0,8,15, 20,5,10,12,0,3, 0,15,20,8,10,5,
                5,0,3,20,15,10, 15,10,8,5,20,0, 8,12,0,15,3,20),
              nrow = 6, ncol = 6, byrow = TRUE)
rownames(otu) <- paste0("ASV", 1:6)
colnames(otu) <- paste0("S", 1:6)
tree <- ape::read.tree(text =
  "(((ASV1:0.1,ASV2:0.2):0.1,ASV3:0.15):0.2,((ASV4:0.1,ASV5:0.2):0.15,ASV6:0.2):0.1);")
meta <- data.frame(sample = colnames(otu),
                   group = c("A","A","A","B","B","B"), stringsAsFactors = FALSE)

cat("Test 1: shape + columns\n")
set.seed(11)
res <- calc_beta_nri(otu, tree, beta_reps = 50, verbose = FALSE)
stopifnot(all(c("from","to","beta_mpd","beta_nri") %in% names(res)))
stopifnot(nrow(res) == 6 * 5 / 2)
stopifnot(!any(is.na(res$beta_nri)))
cat("✓ Test 1 passed\n\n")

cat("Test 2: parity with manual ses.betaMPD (one pair, via picante::comdist)\n")
comm <- t(otu)
coph <- ape::cophenetic.phylo(tree)[colnames(comm), colnames(comm)]
set.seed(11)
obs_mpd <- as.matrix(picante::comdist(comm, coph, abundance.weighted = TRUE))
nulls <- replicate(50, as.matrix(picante::comdist(
  comm, picante::taxaShuffle(coph), abundance.weighted = TRUE)))
nri_manual <- (obs_mpd["S2","S1"] - mean(nulls["S2","S1",])) / sd(nulls["S2","S1",])
got <- res$beta_nri[res$from == "S2" & res$to == "S1"]
stopifnot(length(got) == 1)
stopifnot(isTRUE(all.equal(as.numeric(got), nri_manual, tolerance = 1e-9)))
cat("✓ Test 2 passed\n\n")

cat("Test 3: grouped output\n")
set.seed(11)
res_g <- calc_beta_nri(otu, tree, beta_reps = 50, verbose = FALSE,
                       sample = meta, group_col = "group")
stopifnot("group" %in% names(res_g))
stopifnot(nrow(res_g) == 2 * (3 * 2 / 2))
cat("✓ Test 3 passed\n\n")

cat("Test 4: no tree -> error\n")
stopifnot(inherits(try(calc_beta_nri(otu, NULL, beta_reps = 5,
                                     verbose = FALSE), silent = TRUE), "try-error"))
cat("✓ Test 4 passed\n\n")

cat("=====================================\n")
cat("All calc_beta_nri tests passed! ✓\n")
cat("=====================================\n")
