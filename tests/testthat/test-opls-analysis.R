# Migrated from tests/test_opls_analysis.R (bare script).
# opls_analysis() modes + internal pairwise helpers.

# Three groups of 6 samples each, with a clear group signal.
make_opls_fixture <- function() {
  set.seed(1)
  n_feat <- 30L
  features <- sprintf("F%02d", seq_len(n_feat))
  grp <- rep(c("A", "B", "C"), each = 6)
  Xm <- matrix(rnorm(length(grp) * n_feat), nrow = length(grp),
               dimnames = list(sprintf("S%02d", seq_along(grp)), features))
  for (g in seq_along(unique(grp))) {
    Xm[grp == unique(grp)[g], ] <- Xm[grp == unique(grp)[g], ] + g * 3
  }
  list(Xm = Xm, grp = grp, features = features,
       pairs = utils::combn(unique(grp), 2, simplify = FALSE))  # A-B, A-C, B-C
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("fit_pairwise_opls returns correct structure", {
  fx <- make_opls_fixture()
  fit <- bioRtools:::fit_pairwise_opls(fx$Xm, fx$grp, fx$pairs, fx$features,
                                       ortho_components = 1,
                                       scaling = "standard",
                                       validation = "CV", cv_folds = 5,
                                       verbose = FALSE)

  expect_length(fit$models, 3)   # 3 comparisons
  expect_true(all(c("B vs A", "C vs A", "C vs B") %in% names(fit$models)))
  # Each comparison keeps both groups' samples (6+6 = 12);
  # 3 comparisons -> 36 rows.
  expect_equal(nrow(fit$scores), 36)
  expect_true("comparison" %in% names(fit$scores))
  expect_true("vip" %in% names(fit$vip) && "comparison" %in% names(fit$vip))
  expect_true("group" %in% names(fit$vip) && "ref_group" %in% names(fit$vip))
  # summaries must carry real model-quality values (ropls cols are R2Y(cum)/Q2(cum))
  smry <- fit$summaries[[1]]
  expect_false(is.na(smry$R2Y))
  expect_false(is.na(smry$Q2Y))
})

test_that("zero-variance-in-subset variables are dropped, not errored", {
  fx <- make_opls_fixture()
  Xm_zv <- cbind(fx$Xm, ZV = rep(5, nrow(fx$Xm)))   # constant across all samples
  feat_zv <- c(fx$features, "ZV")
  fit_zv <- bioRtools:::fit_pairwise_opls(Xm_zv, fx$grp, fx$pairs, feat_zv,
                                          ortho_components = 1,
                                          scaling = "standard",
                                          validation = "CV", cv_folds = 5,
                                          verbose = FALSE)

  expect_length(fit_zv$models, 3)
  # ZV is zero-variance in every subset, so excluded from every comparison's VIP
  expect_false("ZV" %in% fit_zv$vip$feature)
})

test_that("compute_pairwise_diff computes log2FC/p-value/regulation", {
  fx <- make_opls_fixture()
  fit <- bioRtools:::fit_pairwise_opls(fx$Xm, fx$grp, fx$pairs, fx$features,
                                       ortho_components = 1,
                                       scaling = "standard",
                                       validation = "CV", cv_folds = 5,
                                       verbose = FALSE)
  diff <- bioRtools:::compute_pairwise_diff(fx$Xm, fx$grp, fx$pairs,
                                            fx$features, fit$vip,
                                            test_method = "t-test",
                                            p_adjust_method = "BH",
                                            p_threshold = 0.05)

  expected_cols <- c("feature", "group", "ref_group", "log2_fc", "p_value",
                     "p_adjust", "comparison", "vip", "regulation")
  expect_true(all(expected_cols %in% names(diff)))
  expect_equal(nrow(diff), length(fx$pairs) * 30)   # one row per pair x variable
  expect_equal(length(unique(diff$comparison)), 3)
  expect_true(all(diff$regulation %in% c("Up", "Down", "NS")))
  # Direction sanity: group B has higher mean than A (shift +2*3 in construction)
  b_vs_a <- diff[diff$comparison == "B vs A", ]
  expect_true(all(b_vs_a$log2_fc > 0))              # B up vs A
})

test_that("pairwise = TRUE produces all-pairs comparisons", {
  fx <- make_opls_fixture()
  sample_meta <- tibble::tibble(sample_id = rownames(fx$Xm), group = fx$grp)
  res_ap <- quiet(opls_analysis(fx$Xm, sample = sample_meta,
                                sample_col = "sample_id", group_col = "group",
                                pairwise = TRUE, verbose = FALSE))

  expect_equal(length(unique(res_ap$vip_scores$comparison)), 3)  # choose(3,2) = 3
  expect_false(is.null(res_ap$models))
  expect_type(res_ap$models, "list")
  expect_length(res_ap$models, 3)
  expect_false(is.null(res_ap$differential_analysis))
  expect_equal(length(unique(res_ap$differential_analysis$comparison)), 3)
})

test_that("default (2 groups, no pairwise) is a single binary model", {
  # ropls OPLS-DA is binary-only, so single-model mode is exercised on 2 groups.
  fx <- make_opls_fixture()
  two_grp <- fx$grp %in% c("A", "B")
  Xm2 <- fx$Xm[two_grp, ]
  meta2 <- tibble::tibble(sample_id = rownames(Xm2), group = fx$grp[two_grp])
  res_single <- quiet(opls_analysis(Xm2, sample = meta2,
                                    sample_col = "sample_id",
                                    group_col = "group", verbose = FALSE))

  expect_null(res_single$models)                     # not pairwise
  expect_null(res_single$differential_analysis)      # no ref_group
  expect_false("comparison" %in% names(res_single$vip_scores))
})

test_that("ref_group produces ref-anchored pairwise + diff", {
  fx <- make_opls_fixture()
  sample_meta <- tibble::tibble(sample_id = rownames(fx$Xm), group = fx$grp)
  res_ref <- quiet(opls_analysis(fx$Xm, sample = sample_meta,
                                 sample_col = "sample_id", group_col = "group",
                                 ref_group = "A", verbose = FALSE))

  expect_equal(length(unique(res_ref$vip_scores$comparison)), 2)  # B vs A, C vs A
  expect_true(all(c("B vs A", "C vs A") %in%
                    unique(res_ref$vip_scores$comparison)))
  expect_false(is.null(res_ref$differential_analysis))
  expect_true(all(unique(res_ref$differential_analysis$ref_group) == "A"))
})

test_that(">2 groups, no ref_group/pairwise -> auto all-pairs", {
  fx <- make_opls_fixture()
  sample_meta <- tibble::tibble(sample_id = rownames(fx$Xm), group = fx$grp)
  res_auto <- quiet(opls_analysis(fx$Xm, sample = sample_meta,
                                  sample_col = "sample_id",
                                  group_col = "group", verbose = FALSE))

  expect_false(is.null(res_auto$models))
  expect_equal(length(unique(res_auto$vip_scores$comparison)), 3)
})

test_that("pairwise_oplsda delegates to opls_analysis(pairwise=TRUE)", {
  fx <- make_opls_fixture()
  sample_meta <- tibble::tibble(sample_id = rownames(fx$Xm), group = fx$grp)
  res_ap <- quiet(opls_analysis(fx$Xm, sample = sample_meta,
                                sample_col = "sample_id", group_col = "group",
                                pairwise = TRUE, verbose = FALSE))
  dep <- quiet(pairwise_oplsda(fx$Xm, sample_meta, group_col = "group"))

  expect_true(all(c("vip_scores", "model_summary", "scores", "models") %in%
                    names(dep)))
  # Same comparisons as the direct all-pairs call
  expect_identical(sort(unique(dep$vip_scores$comparison)),
                   sort(unique(res_ap$vip_scores$comparison)))
})
