test_that("public GlobalPatterns workflow runs end to end", {
  testthat::skip_if_not_installed("mia")
  data(GlobalPatterns, package = "mia")

  X <- t(as.matrix(SummarizedExperiment::assay(GlobalPatterns)))
  meta <- as.data.frame(SummarizedExperiment::colData(GlobalPatterns))
  rownames(X) <- meta$X.SampleID
  meta$sample <- meta$X.SampleID
  rownames(meta) <- meta$X.SampleID
  taxo <- as.data.frame(SummarizedExperiment::rowData(GlobalPatterns))

  alpha <- calc_alpha_diversity(
    X,
    meta,
    index = c("shannon", "observed"),
    verbose = FALSE
  )
  expect_equal(nrow(alpha), nrow(X) * 2L)
  expect_true(all(c("SampleType", "sample", "value") %in% names(alpha)))

  pcoa <- suppressMessages(pcoa_analysis(
    X,
    meta,
    color = "SampleType",
    size = 2,
    alpha = 0.9
  ))
  expect_s3_class(pcoa, "bioRtools_fit")
  expect_equal(nrow(pcoa$tables$scores), nrow(X))

  set.seed(1)
  permanova <- permanova_test(
    X ~ SampleType,
    meta,
    permutations = 99
  )
  expect_s3_class(permanova, "bioRtools_fit")
  expect_gt(permanova$tables$permanova$R2[1], 0.5)

  composition <- suppressMessages(top_tax_n(
    X,
    meta,
    taxo,
    which = "Phylum",
    group_col = "SampleType",
    plot_by = "group",
    plot_type = "stacked_bar",
    ref = "Feces",
    verbose = FALSE
  ))
  expect_s3_class(composition$plot_abundance, "ggplot")
  expect_true(nrow(composition$data_processed) > 0)
})
