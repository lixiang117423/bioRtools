# Migrated from tests/test_orientation.R (bare script).
# Verifies the shared feature_as_row orientation helpers directly, then checks
# each ordination function: a feature x sample matrix (auto-detected AND forced)
# must reproduce the sample x feature result exactly.

detect_feature_as_row <- bioRtools:::detect_feature_as_row
orient_to_sample_row <- bioRtools:::orient_to_sample_row

set.seed(1)
n_sample <- 12L
n_feat <- 30L
samples <- sprintf("S%02d", seq_len(n_sample))
features <- sprintf("F%02d", seq_len(n_feat))
# Two groups: OPLS-DA is binary-only, so the same dataset must work for every fn
grp <- rep(c("A", "B"), each = 6)

# Sample x feature matrix with a clear group signal (stable for supervised models)
Xm <- matrix(rnorm(n_sample * n_feat), nrow = n_sample,
             dimnames = list(samples, features))
for (g in seq_along(unique(grp))) Xm[grp == unique(grp)[g], ] <-
  Xm[grp == unique(grp)[g], ] + g * 3
X <- as.data.frame(Xm)              # rda_analysis requires a data frame
Xf <- as.data.frame(t(Xm))          # features as rows, samples as columns

sample_meta <- tibble::tibble(
  sample_id = samples, group = grp, batch = rep(c("b1", "b2"), 6)
)
physico <- data.frame(
  pH = rnorm(n_sample), temp = rnorm(n_sample), row.names = samples
)
sample_rda <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)

# Numeric coordinate matrix from a scores data frame, rows sorted by sample id
coord_mat <- function(df) {
  id_cols <- intersect(names(df), c("sample", "sample_id", "sample_name", "ID", "id"))
  nums <- df[, sapply(df, is.numeric), drop = FALSE]
  if (length(id_cols)) nums <- nums[order(df[[id_cols[1]]]), , drop = FALSE]
  as.matrix(nums)
}
same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10))
quiet <- function(expr) suppressMessages(suppressWarnings(force(expr)))

test_that("shared orientation helpers", {
  expect_false(detect_feature_as_row(X, sample_meta, "sample_id"))   # sample-row
  expect_true(detect_feature_as_row(Xf, sample_meta, "sample_id"))   # feature-row
  expect_false(detect_feature_as_row(X, NULL, "sample_id"))          # no metadata
  expect_true(detect_feature_as_row(Xf, physico, NULL))              # rda-style (rownames)
  expect_identical(orient_to_sample_row(X, sample_meta, "sample_id", NA, FALSE), X)  # no transpose
  expect_identical(
    orient_to_sample_row(Xf, sample_meta, "sample_id", FALSE, FALSE), Xf)             # forced off
  rn <- rownames(orient_to_sample_row(Xf, sample_meta, "sample_id", NA, FALSE))
  expect_identical(rn, samples)                                        # transposed -> samples
  expect_error(orient_to_sample_row(X, sample_meta, "sample_id", "maybe", FALSE))
})

test_that("pca_analysis: feature-row reproduces sample-row", {
  ref <- quiet(pca_analysis(X, sample_meta, color_by = "group", plot_type = "none"))
  auto <- quiet(pca_analysis(Xf, sample_meta, color_by = "group", plot_type = "none"))
  force <- quiet(pca_analysis(Xf, sample_meta, color_by = "group", plot_type = "none",
                              feature_as_row = TRUE))
  expect_true(same(coord_mat(ref$sample_scores), coord_mat(auto$sample_scores)))
  expect_true(same(coord_mat(ref$sample_scores), coord_mat(force$sample_scores)))
})

test_that("pcoa_analysis: feature-row reproduces sample-row", {
  ref <- quiet(pcoa_analysis(X, sample_meta))
  auto <- quiet(pcoa_analysis(Xf, sample_meta))
  expect_true(same(coord_mat(ref$point_data), coord_mat(auto$point_data)))
})

test_that("opls_analysis: feature-row reproduces sample-row", {
  ref <- quiet(opls_analysis(X, sample_meta, sample_col = "sample_id",
                             group_col = "group", verbose = FALSE))
  auto <- quiet(opls_analysis(Xf, sample_meta, sample_col = "sample_id",
                              group_col = "group", verbose = FALSE))
  expect_true(same(coord_mat(ref$scores), coord_mat(auto$scores)))
})

test_that("rda_analysis: feature-row reproduces sample-row", {
  ref <- quiet(rda_analysis(X, physico, sample_rda))
  auto <- quiet(rda_analysis(Xf, physico, sample_rda))
  expect_true(same(coord_mat(ref$sample_scores), coord_mat(auto$sample_scores)))
})

test_that("tsne_analysis: feature-row reproduces sample-row", {
  ref <- quiet(tsne_analysis(X, sample_meta, perplexity = 3, pca_dims = 10,
                             plot_type = "none"))
  auto <- quiet(tsne_analysis(Xf, sample_meta, perplexity = 3, pca_dims = 10,
                              plot_type = "none"))
  expect_true(same(coord_mat(ref$sample_scores), coord_mat(auto$sample_scores)))
})

test_that("pairwise_oplsda: feature-row reproduces sample-row (VIP + scores)", {
  ref <- quiet(pairwise_oplsda(X, sample_meta, group_col = "group"))
  auto <- quiet(pairwise_oplsda(Xf, sample_meta, group_col = "group"))
  vip_ref <- ref$vip_scores[order(ref$vip_scores$feature), ]
  vip_auto <- auto$vip_scores[order(auto$vip_scores$feature), ]
  expect_true(same(vip_ref, vip_auto))
  expect_equal(nrow(ref$scores), nrow(auto$scores))
})
