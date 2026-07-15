# Regression: calc_beta_nti output must be byte-identical after the
# calc_ses_beta_core refactor.
suppressPackageStartupMessages({ library(ape); library(dplyr); library(tidyr); library(picante) })
source("R/calc_ses_beta_core.R")
source("R/calc_beta_nti.R")

otu <- matrix(c(10,20,5,0,8,15, 20,5,10,12,0,3, 0,15,20,8,10,5,
                5,0,3,20,15,10, 15,10,8,5,20,0, 8,12,0,15,3,20),
              nrow = 6, ncol = 6, byrow = TRUE)
rownames(otu) <- paste0("ASV", 1:6)
colnames(otu) <- paste0("S", 1:6)
tree <- ape::read.tree(text =
  "(((ASV1:0.1,ASV2:0.2):0.1,ASV3:0.15):0.2,((ASV4:0.1,ASV5:0.2):0.15,ASV6:0.2):0.1);")
meta <- data.frame(sample = colnames(otu),
                   group = c("A","A","A","B","B","B"), stringsAsFactors = FALSE)

baseline <- readRDS("tests/baselines/calc_beta_nti_fixture.rds")

cat("Test 1: ungrouped output identical to baseline\n")
set.seed(42)
got <- calc_beta_nti(otu, tree, beta_reps = 50, verbose = FALSE)
stopifnot(identical(got, baseline$all))
cat("Test 1 passed\n\n")

cat("Test 2: grouped output identical to baseline\n")
set.seed(42)
got_g <- calc_beta_nti(otu, tree, beta_reps = 50, verbose = FALSE,
                       sample = meta, group_col = "group")
stopifnot(identical(got_g, baseline$grouped))
cat("Test 2 passed\n\n")

cat("Test 3: tree/data mismatch errors\n")
# Tree with zero tip overlap -> triggers "No matching taxa" stop().
# (A partial mismatch does not error in either old or new code, since
# intersect() keeps the shared taxa, so use fully disjoint labels here.)
bad_tree <- ape::read.tree(text = "((X1:0.1,X2:0.2):0.1,X3:0.15);")
stopifnot(inherits(try(calc_beta_nti(otu, bad_tree, beta_reps = 5,
                                     verbose = FALSE), silent = TRUE), "try-error"))
cat("Test 3 passed\n\n")

cat("=====================================\n")
cat("All calc_beta_nti regression tests passed!\n")
cat("=====================================\n")
