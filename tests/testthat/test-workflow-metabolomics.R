test_that("public nutrimouse workflow runs end to end", {
  data(nutrimouse, package = "mixOmics")
  X <- as.data.frame(nutrimouse$lipid)
  meta <- data.frame(
    sample = rownames(X),
    diet = nutrimouse$diet,
    genotype = nutrimouse$genotype,
    row.names = rownames(X),
    stringsAsFactors = FALSE
  )

  set.seed(4)
  opls <- suppressMessages(suppressWarnings(
    opls_analysis(
      X,
      sample = meta,
      sample_col = "sample",
      group_col = "diet",
      pairwise = TRUE,
      validation = "CV",
      cv_folds = 4,
      verbose = FALSE
    )
  ))
  expect_s3_class(opls, "bioRtools_fit")
  expect_length(opls$model, choose(5, 2))
  expect_true(nrow(opls$tables$differential) > 0)

  volcano <- suppressWarnings(plot_metab_volcano(opls$tables$differential))
  expect_s3_class(volcano, "ggplot")

  set.seed(9)
  fit <- suppressMessages(suppressWarnings(
    spls_analysis(
      X,
      sample = meta,
      sample_col = "sample",
      group_col = "diet",
      ncomp = 2,
      keep_x = c(10, 5),
      verbose = FALSE
    )
  ))
  expect_s3_class(fit, "bioRtools_fit")
  expect_equal(nrow(fit$tables$sample_scores), nrow(X))

  scaled <- scale(as.matrix(X), center = TRUE, scale = TRUE)
  set.seed(9)
  native <- mixOmics::splsda(
    X = scaled,
    Y = as.factor(nutrimouse$diet),
    ncomp = 2,
    keepX = c(10, 5),
    scale = FALSE,
    max.iter = 500,
    tol = 1e-6
  )
  expect_equal(
    unname(as.matrix(fit$tables$sample_scores[, c("comp1", "comp2")])),
    unname(native$variates$X),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(fit$tables$variable_loadings[, c("comp1", "comp2")])),
    unname(native$loadings$X),
    tolerance = 1e-12
  )
})
