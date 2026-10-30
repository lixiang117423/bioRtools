test_that("base correlation fallback agrees with WGCNA", {
  testthat::skip_if_not_installed("WGCNA")

  x <- iris[, c("Sepal.Length", "Sepal.Width")]
  y <- iris[, c("Petal.Length", "Petal.Width")]

  for (method in c("pearson", "spearman", "kendall")) {
    fallback <- bioRtools:::base_cor_and_pvalue(x, y, method = method)
    native <- WGCNA::corAndPvalue(x, y, method = method)

    expect_equal(
      unname(fallback$cor),
      unname(native$cor),
      tolerance = 1e-12
    )
    expect_equal(
      unname(fallback$p),
      unname(native$p),
      tolerance = 1e-10
    )
  }
})
