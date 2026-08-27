# Migrated from tests/test_opls_diff_method.R (bare script).
# VIP-based differential determination in compute_pairwise_diff().
# Standalone: no ropls needed (VIP table is supplied directly).

make_diff_method_fixture <- function() {
  set.seed(2)
  n_feat <- 40L
  features <- sprintf("F%02d", seq_len(n_feat))
  grp <- rep(c("A", "B", "C"), each = 6)
  # Positive baseline (metabolomics-like intensities), clear group shifts
  Xm <- matrix(rnorm(length(grp) * n_feat, mean = 5, sd = 1),
               nrow = length(grp),
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
  list(Xm = Xm, grp = grp, pairs = pairs, features = features, vip = vip)
}

test_that("VIP-based labeling (diff_method = 'vip')", {
  fx <- make_diff_method_fixture()
  d_vip <- bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                             fx$features, fx$vip,
                                             test_method = "t-test",
                                             p_adjust_method = "BH",
                                             diff_method = "vip",
                                             vip_threshold = 1)
  # VIP >= 1 -> Up/Down (direction from log2FC); VIP < 1 -> NS regardless of p
  expect_true(all(d_vip$regulation[d_vip$vip >= 1] %in% c("Up", "Down")))
  expect_true(all(d_vip$regulation[d_vip$vip < 1] == "NS"))
  expect_true(all(d_vip$log2_fc[d_vip$regulation == "Up"] > 0))
  expect_true(all(d_vip$log2_fc[d_vip$regulation == "Down"] < 0))
  # is_differential is consistent with regulation and with VIP threshold
  expect_true(all(d_vip$is_differential ==
                    (d_vip$regulation %in% c("Up", "Down"))))
  expect_true(all(d_vip$is_differential == (d_vip$vip >= 1)))
  # p-value columns are still computed for plotting
  expect_false(any(is.na(d_vip$p_adjust)))
  expect_equal(sum(d_vip$is_differential), 40)  # 20 important x 2 comparisons
})

test_that("p-value-based labeling (legacy)", {
  fx <- make_diff_method_fixture()
  d_p <- bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                           fx$features, fx$vip,
                                           test_method = "t-test",
                                           p_adjust_method = "BH",
                                           diff_method = "pvalue",
                                           p_threshold = 0.05)

  expect_true(all(d_p$regulation %in% c("Up", "Down", "NS")))
  expect_true(all(d_p$p_adjust[d_p$regulation == "Up"] < 0.05))
  expect_true(all(d_p$p_adjust[d_p$regulation == "Down"] < 0.05))
  expect_true(all(d_p$p_adjust[d_p$regulation == "NS"] >= 0.05))
  expect_true(all(d_p$is_differential ==
                    (d_p$regulation %in% c("Up", "Down"))))
})

test_that("vip_threshold honored", {
  fx <- make_diff_method_fixture()
  d_15 <- bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                            fx$features, fx$vip,
                                            test_method = "t-test",
                                            diff_method = "vip",
                                            vip_threshold = 1.5)

  expect_true(all(d_15$regulation[d_15$vip < 1.5] == "NS"))
  expect_true(all(d_15$regulation[d_15$vip >= 1.5] %in% c("Up", "Down")))
  expect_equal(sum(d_15$is_differential), 40)
})

# (Original "Test 4: methods can disagree" only printed an informational
# message and asserted nothing, so it is not carried over.)

test_that("default diff_method is vip", {
  fx <- make_diff_method_fixture()
  d_vip <- bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                             fx$features, fx$vip,
                                             test_method = "t-test",
                                             p_adjust_method = "BH",
                                             diff_method = "vip",
                                             vip_threshold = 1)
  d_def <- bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                             fx$features, fx$vip,
                                             test_method = "t-test")

  expect_identical(d_def$regulation, d_vip$regulation)
})

test_that("invalid diff_method errors", {
  fx <- make_diff_method_fixture()
  expect_error(bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                                 fx$features, fx$vip,
                                                 test_method = "t-test",
                                                 diff_method = "bogus"))
})

test_that("non-finite log2FC (negative data) labels NS, not NA", {
  fx <- make_diff_method_fixture()
  Xm2 <- fx$Xm
  Xm2[fx$grp == "A", 1] <- -5    # F01: negative ref mean -> NaN log2_fc
  vip2 <- fx$vip
  vip2$vip[vip2$feature == "F01"] <- 1.5
  # negative ref mean makes log() emit "NaNs produced" warnings by design
  d2 <- suppressWarnings(bioRtools:::compute_pairwise_diff(Xm2, fx$grp, fx$pairs,
                                          fx$features, vip2,
                                          test_method = "t-test",
                                          diff_method = "vip"))

  expect_true(all(d2$regulation[d2$feature == "F01"] == "NS"))
  expect_false(any(is.na(d2$regulation)))
})
