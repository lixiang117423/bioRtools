# Migrated from tests/test_plot_multi_volcano.R (bare script).
# plot_multi_volcano() with synthetic data, so it runs in R CMD check without
# external files.

# Synthetic data with the auto-detected column names (gene / avg_log2FC /
# p_val_adj / cluster), mirroring the Seurat-style input the function targets
set.seed(1)
test_data <- data.frame(
  gene       = paste0("Gene", 1:600),
  avg_log2FC = stats::rnorm(600, 0, 1.5),
  p_val_adj  = stats::runif(600, 0, 1),
  cluster    = sample(paste0("Cluster", 1:4), 600, replace = TRUE),
  stringsAsFactors = FALSE
)

test_that("basic usage auto-detects the columns and returns a ggplot", {
  p1 <- plot_multi_volcano(test_data)
  expect_true(ggplot2::is_ggplot(p1))
})

test_that("custom thresholds and top_n", {
  p2 <- plot_multi_volcano(
    test_data,
    log2fc_thr = 1.5,
    padj_thr = 0.01,
    top_n = 3
  )
  expect_true(ggplot2::is_ggplot(p2))
})

test_that("custom colors", {
  p3 <- plot_multi_volcano(
    test_data,
    up_color = "#FF6B6B",
    down_color = "#6B6BFF",
    group_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  )
  expect_true(ggplot2::is_ggplot(p3))
})

test_that("custom axis limits", {
  p4 <- plot_multi_volcano(
    test_data,
    y_limits = c(-5, 10),
    y_breaks = c(-5, -2, 0, 2, 5, 10)
  )
  expect_true(ggplot2::is_ggplot(p4))
})

test_that("undetectable columns raise an error", {
  expect_error(
    plot_multi_volcano(data.frame(a = 1:3, b = 1:3)),
    "Cannot auto-detect group column"
  )
})

test_that("plot can be saved to a file", {
  p1 <- plot_multi_volcano(test_data)
  out_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(out_file, p1, width = 10, height = 6, dpi = 100)
  expect_true(file.exists(out_file))
})
