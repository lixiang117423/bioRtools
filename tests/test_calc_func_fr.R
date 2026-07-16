# Tests for calc_func_fr() + calc_func_fr_comm().
source("R/calc_func_fr.R")

# 3 ASVs x 3 functions (ASV x function):
asv_func <- matrix(c(1, 1, 0,
                     0, 1, 1,
                     1, 0, 1), nrow = 3, ncol = 3)
rownames(asv_func) <- c("ASV1", "ASV2", "ASV3")
colnames(asv_func) <- c("funcA", "funcB", "funcC")
# funcA = ASV1,ASV2 ; funcB = ASV2,ASV3 ; funcC = ASV1,ASV3
abundance <- matrix(c(10, 20, 30,
                      5,  5,  5), nrow = 3, ncol = 2,
                    dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2")))

cat("Test 1: unweighted FR in [0,1]; funcA/S1 = 2/3\n")
fr <- calc_func_fr(asv_func, abundance, abundance_weighted = FALSE)
stopifnot(all(fr$fr >= 0 & fr$fr <= 1 + 1e-9, na.rm = TRUE))
fa <- fr$fr[fr$sample == "S1" & fr[["function"]] == "funcA"]
stopifnot(isTRUE(all.equal(as.numeric(fa), 2 / 3)))
cat("✓ Test 1 passed\n\n")

cat("Test 2: weighted funcA/S1 = (10+20)/(10+20+30) = 0.5\n")
frw <- calc_func_fr(asv_func, abundance, abundance_weighted = TRUE)
faw <- frw$fr[frw$sample == "S1" & frw[["function"]] == "funcA"]
stopifnot(isTRUE(all.equal(as.numeric(faw), 30 / 60)))
cat("✓ Test 2 passed\n\n")

cat("Test 3: a function held by all present taxa -> FR = 1\n")
asv_func2 <- cbind(asv_func, funcD = c(1, 1, 1))
fr2 <- calc_func_fr(asv_func2, abundance)
fd <- fr2$fr[fr2$sample == "S1" & fr2[["function"]] == "funcD"]
stopifnot(isTRUE(all.equal(as.numeric(fd), 1)))
cat("✓ Test 3 passed\n\n")

cat("Test 4: calc_func_fr_comm — one row per sample; geometric mean\n")
comm <- calc_func_fr_comm(fr2)
stopifnot(nrow(comm) == 2)
stopifnot(all(c("sample", "fr_comm") %in% names(comm)))
s1_frs <- fr2$fr[fr2$sample == "S1"]
manual <- exp(mean(log(pmax(s1_frs, 1e-6))))
stopifnot(isTRUE(all.equal(comm$fr_comm[comm$sample == "S1"], manual)))
cat("✓ Test 4 passed\n\n")

cat("Test 5: adj_tax=TRUE without taxonomy -> error\n")
stopifnot(inherits(try(calc_func_fr(asv_func, abundance, adj_tax = TRUE),
                       silent = TRUE), "try-error"))
cat("✓ Test 5 passed\n\n")

cat("=====================================\n")
cat("All calc_func_fr tests passed! ✓\n")
cat("=====================================\n")
