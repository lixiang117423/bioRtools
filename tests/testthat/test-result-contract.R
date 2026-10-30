make_contract_fixture <- function() {
  set.seed(42)
  n_samples <- 12L
  n_features <- 10L
  group <- rep(c("A", "B"), each = 6L)
  X <- matrix(
    rnorm(n_samples * n_features),
    nrow = n_samples,
    dimnames = list(paste0("S", seq_len(n_samples)), paste0("F", seq_len(n_features)))
  )
  X[group == "B", ] <- X[group == "B", ] + 1.5
  meta <- data.frame(
    sample = rownames(X),
    group = group,
    row.names = rownames(X),
    stringsAsFactors = FALSE
  )
  physicochemical <- data.frame(
    batch = as.numeric(group == "B"),
    row.names = rownames(X)
  )
  list(X = X, meta = meta, physicochemical = physicochemical)
}

expect_common_contract <- function(result, engine) {
  expect_s3_class(result, "bioRtools_fit")
  expect_false(is.null(result$model))
  expect_true(is.list(result$tables))
  expect_true(length(result$tables) > 0)
  expect_true(is.list(result$plots))
  expect_true(is.list(result$provenance))
  expect_identical(result$provenance$engine, engine)
  expect_true(nzchar(result$provenance$package_version))
  expect_true(nzchar(result$provenance$r_version))
}

test_that("core multivariate workflows return the shared result contract", {
  fx <- make_contract_fixture()

  pca <- suppressMessages(suppressWarnings(
    pca_analysis(fx$X, fx$meta, color_by = "group", plot_type = "none")
  ))
  expect_common_contract(pca, "FactoMineR::PCA")
  expect_equal(pca$tables$sample_scores$sample, rownames(fx$X))

  pcoa <- suppressMessages(suppressWarnings(pcoa_analysis(fx$X, fx$meta)))
  expect_common_contract(pcoa, "vegan::vegdist + ape::pcoa")
  expect_equal(nrow(pcoa$tables$scores), nrow(fx$X))

  permanova <- suppressMessages(suppressWarnings(
    permanova_test(fx$X ~ group, fx$meta, permutations = 99)
  ))
  expect_common_contract(permanova, "vegan::adonis2")
  expect_true("permanova" %in% names(permanova$tables))

  spls <- suppressMessages(suppressWarnings(
    spls_analysis(
      fx$X,
      sample = fx$meta,
      sample_col = "sample",
      group_col = "group",
      ncomp = 2,
      keep_x = c(4, 3),
      verbose = FALSE
    )
  ))
  expect_common_contract(spls, "mixOmics::splsda")
  expect_equal(nrow(spls$tables$sample_scores), nrow(fx$X))

  opls <- suppressMessages(suppressWarnings(
    opls_analysis(
      fx$X,
      sample = fx$meta,
      sample_col = "sample",
      group_col = "group",
      validation = "CV",
      cv_folds = 3,
      verbose = FALSE
    )
  ))
  expect_common_contract(opls, "ropls::opls")
  expect_equal(nrow(opls$tables$scores), nrow(fx$X))

  expect_output(print(pcoa), "<bioRtools_fit>")
  expect_s3_class(summary(pcoa), "summary.bioRtools_fit")
  expect_s3_class(tibble::as_tibble(pcoa), "tbl_df")
})

test_that("as_tibble.bioRtools_fit validates the requested table", {
  fx <- make_contract_fixture()
  pcoa <- suppressMessages(suppressWarnings(pcoa_analysis(fx$X, fx$meta)))

  scores <- tibble::as_tibble(pcoa, which = "scores")
  expect_equal(nrow(scores), nrow(fx$X))
  expect_error(tibble::as_tibble(pcoa, which = "missing"), "must name one of")
  pcoa$tables$bad <- 1:3
  expect_error(tibble::as_tibble(pcoa, which = "bad"), "not a data frame")
})
