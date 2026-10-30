make_alignment_fixture <- function() {
  set.seed(7)
  samples <- paste0("S", 1:10)
  groups <- rep(c("A", "B"), each = 5)
  X <- matrix(
    rpois(length(samples) * 12, lambda = 20),
    nrow = length(samples),
    dimnames = list(samples, paste0("F", 1:12))
  )
  X[groups == "B", 1:4] <- X[groups == "B", 1:4] + 10
  meta <- data.frame(
    SampleID = samples,
    group = groups,
    batch = rep(c("b1", "b2"), 5),
    row.names = samples,
    stringsAsFactors = FALSE
  )
  list(X = X, meta = meta)
}

test_that("calc_alpha_diversity accepts a non-'sample' first ID column", {
  fx <- make_alignment_fixture()
  res <- calc_alpha_diversity(
    fx$X,
    fx$meta,
    index = c("shannon", "observed"),
    verbose = FALSE
  )

  expect_equal(nrow(res), nrow(fx$X) * 2)
  expect_true(all(c("SampleID", "group", "batch") %in% names(res)))
  expect_false(any(grepl("\\.x$|\\.y$", names(res))))
  expect_setequal(res$sample, rownames(fx$X))
})

test_that("permanova_test aligns metadata by sample ID instead of row order", {
  fx <- make_alignment_fixture()
  set.seed(11)
  aligned <- permanova_test(fx$X ~ group, fx$meta, permutations = 199)
  set.seed(11)
  shuffled <- permanova_test(
    fx$X ~ group,
    fx$meta[sample(nrow(fx$meta)), , drop = FALSE],
    permutations = 199
  )

  expect_equal(
    aligned$result_permanova$R2,
    shuffled$result_permanova$R2,
    tolerance = 1e-12
  )
  expect_equal(
    aligned$result_permanova[["Pr(>F)"]],
    shuffled$result_permanova[["Pr(>F)"]],
    tolerance = 1e-12
  )
})

test_that("metadata alignment rejects incomplete and duplicate IDs", {
  fx <- make_alignment_fixture()
  missing <- fx$meta[-1, , drop = FALSE]
  expect_error(
    permanova_test(fx$X ~ group, missing, permutations = 99),
    "Number of response rows|do not match"
  )

  duplicate <- fx$meta
  duplicate$SampleID[2] <- duplicate$SampleID[1]
  expect_error(
    calc_alpha_diversity(fx$X, duplicate, index = "shannon", verbose = FALSE),
    "Duplicate sample IDs"
  )
})

test_that("permanova_test requires response row names", {
  fx <- make_alignment_fixture()
  no_names <- unname(fx$X)
  expect_error(
    permanova_test(no_names ~ group, fx$meta, permutations = 99),
    "non-default sample IDs"
  )
})
