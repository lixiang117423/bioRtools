# Migrated from tests/test_calc_func_fr.R (bare script).
# calc_func_fr() + calc_func_fr_comm().

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

test_that("unweighted FR in [0,1]; funcA/S1 = 2/3", {
  fr <- calc_func_fr(asv_func, abundance, abundance_weighted = FALSE)

  expect_true(all(fr$fr >= 0 & fr$fr <= 1 + 1e-9, na.rm = TRUE))
  fa <- fr$fr[fr$sample == "S1" & fr[["function"]] == "funcA"]
  expect_equal(as.numeric(fa), 2 / 3, tolerance = 1e-8)
})

test_that("weighted funcA/S1 = (10+20)/(10+20+30) = 0.5", {
  frw <- calc_func_fr(asv_func, abundance, abundance_weighted = TRUE)
  faw <- frw$fr[frw$sample == "S1" & frw[["function"]] == "funcA"]

  expect_equal(as.numeric(faw), 30 / 60, tolerance = 1e-8)
})

test_that("a function held by all present taxa has FR = 1", {
  asv_func2 <- cbind(asv_func, funcD = c(1, 1, 1))
  fr2 <- calc_func_fr(asv_func2, abundance)
  fd <- fr2$fr[fr2$sample == "S1" & fr2[["function"]] == "funcD"]

  expect_equal(as.numeric(fd), 1, tolerance = 1e-8)
})

test_that("calc_func_fr_comm returns one row per sample; geometric mean", {
  asv_func2 <- cbind(asv_func, funcD = c(1, 1, 1))
  fr2 <- calc_func_fr(asv_func2, abundance)
  comm <- calc_func_fr_comm(fr2)

  expect_equal(nrow(comm), 2)
  expect_true(all(c("sample", "fr_comm") %in% names(comm)))
  s1_frs <- fr2$fr[fr2$sample == "S1"]
  manual <- exp(mean(log(pmax(s1_frs, 1e-6))))
  expect_equal(comm$fr_comm[comm$sample == "S1"], manual, tolerance = 1e-8)
})

test_that("adj_tax=TRUE without taxonomy produces an error", {
  expect_error(calc_func_fr(asv_func, abundance, adj_tax = TRUE))
})
