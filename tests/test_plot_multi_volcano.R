# Test script for plot_multi_volcano function
# This script tests plot_multi_volcano with synthetic data, so it runs in
# R CMD check without external files.

suppressPackageStartupMessages(library(bioRtools))
library(dplyr)

# Synthetic data with the auto-detected column names (gene / avg_log2FC /
# p_val_adj / cluster), mirroring the Seurat-style input the function targets
set.seed(1)
test_data <- data.frame(
  gene       = paste0("Gene", 1:600),
  avg_log2FC = rnorm(600, 0, 1.5),
  p_val_adj  = runif(600, 0, 1),
  cluster    = sample(paste0("Cluster", 1:4), 600, replace = TRUE),
  stringsAsFactors = FALSE
)

is_ggplot <- function(x) inherits(x, "ggplot")

cat("Test 1: Basic usage (auto-detected columns)\n")
p1 <- plot_multi_volcano(test_data)
stopifnot(is_ggplot(p1))
cat("✓ Test 1 passed\n\n")

cat("Test 2: Custom thresholds and top_n\n")
p2 <- plot_multi_volcano(
  test_data,
  log2fc_thr = 1.5,
  padj_thr = 0.01,
  top_n = 3
)
stopifnot(is_ggplot(p2))
cat("✓ Test 2 passed\n\n")

cat("Test 3: Custom colors\n")
p3 <- plot_multi_volcano(
  test_data,
  up_color = "#FF6B6B",
  down_color = "#6B6BFF",
  group_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
)
stopifnot(is_ggplot(p3))
cat("✓ Test 3 passed\n\n")

cat("Test 4: Custom axis limits\n")
p4 <- plot_multi_volcano(
  test_data,
  y_limits = c(-5, 10),
  y_breaks = c(-5, -2, 0, 2, 5, 10)
)
stopifnot(is_ggplot(p4))
cat("✓ Test 4 passed\n\n")

cat("Test 5: Error handling - undetectable columns\n")
tryCatch({
  plot_multi_volcano(data.frame(a = 1:3, b = 1:3))
  cat("✗ Test 5 failed: should have thrown an error\n")
  stop("Test 5 failed: no error thrown")
}, error = function(e) {
  cat("✓ Test 5 passed: correctly caught error -", conditionMessage(e), "\n\n")
})

cat("Test 6: Save plot\n")
out_file <- tempfile(fileext = ".png")
ggplot2::ggsave(out_file, p1, width = 10, height = 6, dpi = 100)
stopifnot(file.exists(out_file))
cat("✓ Test 6 passed: plot saved to", out_file, "\n\n")

cat("=====================================\n")
cat("All plot_multi_volcano tests passed! ✓\n")
cat("=====================================\n")
