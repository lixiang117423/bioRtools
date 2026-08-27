# Migrated from tests/test_calc_assembly_process.R (bare script).
# calc_assembly_process(): 5-process Stegen classification.

proc_levels <- c("Variable selection", "Homogeneous selection",
                 "Dispersal limitation", "Homogeneous dispersal", "Drift")

# Hand-built inputs with known betaNTI/RCbray to exercise all 5 buckets.
make_assembly_inputs <- function() {
  beta_nti <- data.frame(
    from = c("S2","S3","S4","S5","S6","S3","S4","S5","S6","S4","S5","S6","S5","S6","S6"),
    to   = c("S1","S1","S1","S1","S1","S2","S2","S2","S2","S3","S3","S3","S4","S4","S5"),
    beta_nti = c( 3.0,  -3.0, 0.5,  0.0, 0.2,  2.5, -2.5, 0.9, -0.9, 1.0, -1.0, 0.3, 1.5, -1.5, 0.1),
    stringsAsFactors = FALSE
  )
  rc_bray <- data.frame(
    from = beta_nti$from, to = beta_nti$to,
    # S2-S1(bNTI 3 -> var sel), S3-S1(bNTI -3 -> hom sel),
    # S4-S1(bNTI .5, rc .97 -> disp lim), S5-S1(bNTI 0, rc -.97 -> hom disp),
    # S6-S1(bNTI .2, rc .1 -> drift); rest arbitrary within |bNTI|<=2
    rc_bray = c(0.0, 0.0, 0.97, -0.97, 0.10,
                0.98, -0.98, 0.5, -0.5, 0.99, -0.99, 0.0, 0.97, -0.97, 0.0),
    stringsAsFactors = FALSE
  )
  list(beta_nti = beta_nti, rc_bray = rc_bray)
}

test_that("the 5 known pairs map to the 5 distinct processes", {
  inputs <- make_assembly_inputs()
  res <- calc_assembly_process(inputs$beta_nti, inputs$rc_bray, verbose = FALSE)

  expect_type(res, "list")
  expect_true(all(c("pair_process", "summary") %in% names(res)))

  pp <- res$pair_process
  pick <- function(a, b) as.character(pp$process[pp$from == a & pp$to == b])
  expect_equal(pick("S2", "S1"), "Variable selection")
  expect_equal(pick("S3", "S1"), "Homogeneous selection")
  expect_equal(pick("S4", "S1"), "Dispersal limitation")
  expect_equal(pick("S5", "S1"), "Homogeneous dispersal")
  expect_equal(pick("S6", "S1"), "Drift")
})

test_that("summary percentages sum to 100", {
  inputs <- make_assembly_inputs()
  res <- calc_assembly_process(inputs$beta_nti, inputs$rc_bray, verbose = FALSE)

  expect_equal(sum(res$summary$percentage), 100, tolerance = 1e-8)
  expect_setequal(res$summary$process, proc_levels)
})

test_that("betaNRI auto-detection", {
  inputs <- make_assembly_inputs()
  beta_nri <- inputs$beta_nti
  names(beta_nri)[3] <- "beta_nri"
  res2 <- calc_assembly_process(beta_nri, inputs$rc_bray, verbose = FALSE)

  expect_equal(nrow(res2$pair_process), nrow(beta_nri))
})

test_that("custom thresholds", {
  inputs <- make_assembly_inputs()
  res3 <- calc_assembly_process(inputs$beta_nti, inputs$rc_bray,
                                thresholds = c(phylo = 10, comm = 0.99),
                                verbose = FALSE)

  expect_false("Variable selection" %in% as.character(res3$pair_process$process))
})

test_that("grouped summary", {
  inputs <- make_assembly_inputs()
  beta_nti <- inputs$beta_nti
  rc_bray <- inputs$rc_bray
  beta_nti$group <- c(rep("A", 5), rep("B", 10))
  rc_bray$group  <- c(rep("A", 5), rep("B", 10))
  res_g <- calc_assembly_process(beta_nti, rc_bray, verbose = FALSE)

  expect_true("group" %in% names(res_g$summary))
  expect_equal(length(unique(res_g$summary$group)), 2)

  # per-group percentages sum to 100 within each group
  s <- split(res_g$summary$percentage, res_g$summary$group)
  expect_equal(sum(s$A), 100, tolerance = 1e-8)
  expect_equal(sum(s$B), 100, tolerance = 1e-8)
})

test_that("errors", {
  inputs <- make_assembly_inputs()
  beta_nti <- inputs$beta_nti
  rc_bray <- inputs$rc_bray

  no_col <- beta_nti
  no_col$beta_nti <- NULL
  expect_error(calc_assembly_process(no_col, rc_bray))

  rc_half <- rc_bray[1:3, ]
  expect_error(calc_assembly_process(beta_nti, rc_half))

  # group mismatch: rc_bray grouped but beta_nti not (state carried over from
  # the grouped test in the original bare script)
  rc_grp <- rc_bray
  rc_grp$group <- c(rep("A", 5), rep("B", 10))
  expect_error(calc_assembly_process(beta_nti, rc_grp))
})
