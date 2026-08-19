cv_test <- data.frame(
  K = 1:10,
  CV_error = c(0.58, 0.52, 0.50, 0.49, 0.485, 0.482, 0.481, 0.4805,
               0.4802, 0.4801)
)

test_that("find_admixture_k returns the expected K picks", {
  res <- find_admixture_k(cv_test)

  expect_s3_class(res, "admixture_k_results")
  expect_equal(res$best_k, 3)
  expect_equal(res$max_curvature_k, 2)
  expect_equal(res$threshold_k, c(rel_pct_lt_3 = 4, rel_pct_lt_2 = 5))
})

test_that("find_admixture_k tables have the right shape and columns", {
  res <- find_admixture_k(cv_test)

  expect_equal(names(res$kneedle_table), c("K", "CV_error", "distance"))
  expect_equal(nrow(res$kneedle_table), 10)
  expect_equal(names(res$diff_table), c("K_from", "K_to", "d1"))
  expect_equal(nrow(res$diff_table), 9)
  expect_equal(names(res$curvature_table), c("K", "d2"))
  expect_equal(nrow(res$curvature_table), 8)
  expect_equal(names(res$rel_change_table), c("K_from", "K_to", "rel_pct"))
  expect_equal(nrow(res$rel_change_table), 9)
})

test_that("find_admixture_k builds the diagnostic plots", {
  res <- find_admixture_k(cv_test)
  expect_named(res$plots, c("cv_plot", "distance_plot", "combined_plot"))
  expect_s3_class(res$plots$cv_plot, "ggplot")
  expect_s3_class(res$plots$distance_plot, "ggplot")
  expect_s3_class(res$plots$combined_plot, "ggplot")

  # combined_plot carries a secondary (dual) y-axis
  p <- res$plots$combined_plot
  has_sec <- any(vapply(p$scales$scales, function(sc) {
    !is.null(sc$secondary.axis) && !inherits(sc$secondary.axis, "waiver")
  }, logical(1)))
  expect_true(has_sec)

  res_np <- find_admixture_k(cv_test, plot = FALSE)
  expect_null(res_np$plots)
})

test_that("find_admixture_k reads a file path via read_data", {
  dir <- tempfile("cv")
  dir.create(dir)
  f <- file.path(dir, "cv_results.csv")
  readr::write_csv(cv_test, f)

  res <- find_admixture_k(f)
  expect_equal(res$best_k, 3)
})

test_that("find_admixture_k sorts unsorted K values", {
  shuffled <- cv_test[c(5, 1, 9, 3, 7, 2, 4, 6, 8, 10), ]
  res <- find_admixture_k(shuffled)
  expect_equal(res$best_k, 3)
})

test_that("find_admixture_k honors custom column names", {
  df <- data.frame(ancestry = cv_test$K, cv = cv_test$CV_error)
  res <- find_admixture_k(df, k_col = "ancestry", cv_col = "cv")
  expect_equal(res$best_k, 3)
  expect_equal(names(res$kneedle_table), c("K", "CV_error", "distance"))
})

test_that("find_admixture_k reports custom relative thresholds", {
  res <- find_admixture_k(cv_test, rel_thresholds = c(2, 1))
  expect_equal(names(res$threshold_k), c("rel_pct_lt_2", "rel_pct_lt_1"))
  expect_equal(res$threshold_k, c(rel_pct_lt_2 = 5, rel_pct_lt_1 = 6))

  # a threshold no K ever reaches yields NA plus a warning
  expect_warning(
    res_na <- find_admixture_k(cv_test, rel_thresholds = 0.01),
    "never fell below"
  )
  expect_true(is.na(res_na$threshold_k[["rel_pct_lt_0.01"]]))
})

test_that("find_admixture_k print method shows the key picks", {
  res <- find_admixture_k(cv_test)
  expect_output(print(res), "Best K")
  expect_output(print(res), "Max curvature K")
})

test_that("find_admixture_k validates inputs", {
  expect_error(find_admixture_k(42), "must be a data frame")
  expect_error(find_admixture_k(cv_test[, 1, drop = FALSE]),
               "must name a column")
  expect_error(find_admixture_k(data.frame(K = 1:2, CV_error = c(1, 2))),
               "At least 3 K values")
  expect_error(find_admixture_k(data.frame(K = c(1, 1, 2), CV_error = c(1, 2, 3))),
               "Duplicate K values")
  expect_error(find_admixture_k(data.frame(K = 1:3, CV_error = c(1, 1, 1))),
               "constant across K")
  expect_error(find_admixture_k(cv_test, rel_thresholds = -1),
               "positive numbers")
  expect_error(find_admixture_k(cv_test, plot = "yes"), "TRUE or FALSE")
})
