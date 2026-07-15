# Tests for calc_assembly_process(): 5-process Stegen classification.
source("R/calc_assembly_process.R")
suppressPackageStartupMessages({ library(magrittr) })

proc_levels <- c("Variable selection","Homogeneous selection",
                 "Dispersal limitation","Homogeneous dispersal","Drift")

# Hand-built inputs with known betaNTI/RCbray to exercise all 5 buckets.
beta_nti <- data.frame(
  from = c("S2","S3","S4","S5","S6","S3","S4","S5","S6","S4","S5","S6","S5","S6","S6"),
  to   = c("S1","S1","S1","S1","S1","S2","S2","S2","S2","S3","S3","S3","S4","S4","S5"),
  beta_nti = c( 3.0,  -3.0, 0.5,  0.0, 0.2,  2.5, -2.5, 0.9, -0.9, 1.0, -1.0, 0.3, 1.5, -1.5, 0.1),
  stringsAsFactors = FALSE)
rc_bray <- data.frame(
  from = beta_nti$from, to = beta_nti$to,
  # S2-S1(bNTI 3 -> var sel), S3-S1(bNTI -3 -> hom sel),
  # S4-S1(bNTI .5, rc .97 -> disp lim), S5-S1(bNTI 0, rc -.97 -> hom disp),
  # S6-S1(bNTI .2, rc .1 -> drift); rest arbitrary within |bNTI|<=2
  rc_bray = c(0.0, 0.0, 0.97, -0.97, 0.10,
              0.98, -0.98, 0.5, -0.5, 0.99, -0.99, 0.0, 0.97, -0.97, 0.0),
  stringsAsFactors = FALSE)

cat("Test 1: the 5 known pairs map to the 5 distinct processes\n")
res <- calc_assembly_process(beta_nti, rc_bray, verbose = FALSE)
stopifnot(inherits(res, "list"))
stopifnot(all(c("pair_process","summary") %in% names(res)))
pp <- res$pair_process
pick <- function(a, b) as.character(pp$process[pp$from == a & pp$to == b])
stopifnot(pick("S2","S1") == "Variable selection")
stopifnot(pick("S3","S1") == "Homogeneous selection")
stopifnot(pick("S4","S1") == "Dispersal limitation")
stopifnot(pick("S5","S1") == "Homogeneous dispersal")
stopifnot(pick("S6","S1") == "Drift")
cat("✓ Test 1 passed\n\n")

cat("Test 2: summary percentages sum to 100\n")
stopifnot(isTRUE(all.equal(sum(res$summary$percentage), 100)))
stopifnot(setequal(res$summary$process, proc_levels))   # all 5 levels present
cat("✓ Test 2 passed\n\n")

cat("Test 3: betaNRI auto-detection\n")
beta_nri <- beta_nti; names(beta_nri)[3] <- "beta_nri"
res2 <- calc_assembly_process(beta_nri, rc_bray, verbose = FALSE)
stopifnot(nrow(res2$pair_process) == nrow(beta_nri))
cat("✓ Test 3 passed\n\n")

cat("Test 4: custom thresholds\n")
res3 <- calc_assembly_process(beta_nti, rc_bray,
                              thresholds = c(phylo = 10, comm = 0.99), verbose = FALSE)
stopifnot(!"Variable selection" %in% as.character(res3$pair_process$process)) # bNTI 3 < 10 -> no var sel
cat("✓ Test 4 passed\n\n")

cat("Test 5: grouped summary\n")
beta_nti$group <- c(rep("A",5), rep("B",10))
rc_bray$group  <- c(rep("A",5), rep("B",10))
res_g <- calc_assembly_process(beta_nti, rc_bray, verbose = FALSE)
stopifnot("group" %in% names(res_g$summary))
stopifnot(length(unique(res_g$summary$group)) == 2)
# per-group percentages sum to 100 within each group
s <- split(res_g$summary$percentage, res_g$summary$group)
stopifnot(isTRUE(all.equal(sum(s$A), 100)))
stopifnot(isTRUE(all.equal(sum(s$B), 100)))
cat("✓ Test 5 passed\n\n")

cat("Test 6: errors\n")
no_col <- beta_nti; no_col$beta_nti <- NULL
stopifnot(inherits(try(calc_assembly_process(no_col, rc_bray), silent = TRUE), "try-error"))
rc_half <- rc_bray[1:3,]
stopifnot(inherits(try(calc_assembly_process(beta_nti, rc_half), silent = TRUE), "try-error"))
beta_nti$group <- NULL   # group mismatch
stopifnot(inherits(try(calc_assembly_process(beta_nti, rc_bray), silent = TRUE), "try-error"))
cat("✓ Test 6 passed\n\n")

cat("=====================================\n")
cat("All calc_assembly_process tests passed! ✓\n")
cat("=====================================\n")
