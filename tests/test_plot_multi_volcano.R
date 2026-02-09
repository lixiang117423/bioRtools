# Test script for plot_multi_volcano function
# This script tests the new plot_multi_volcano function

# Source the function directly (for testing before package installation)
source("R/plot_multi_volcano.R")

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggprism)

# Read test data
test_data <- read_tsv("/Users/lixiang/Downloads/PlotCase/多组火山图/data.tsv")

cat("Data loaded successfully!\n")
cat("Dimensions:", nrow(test_data), "rows,", ncol(test_data), "columns\n")
cat("Columns:", paste(names(test_data), collapse = ", "), "\n\n")

# Test 1: Basic usage
cat("Test 1: Basic usage\n")
p1 <- plot_multi_volcano(
  test_data,
  fc_column = "avg_log2FC",
  cluster_column = "cluster",
  gene_column = "gene"
)
print(p1)
cat("✓ Test 1 passed\n\n")

# Test 2: Custom thresholds
cat("Test 2: Custom thresholds\n")
p2 <- plot_multi_volcano(
  test_data,
  fc_threshold = 1.5,
  pval_threshold = 0.01,
  label_n = 3
)
print(p2)
cat("✓ Test 2 passed\n\n")

# Test 3: Custom colors (need enough colors for all clusters)
cat("Test 3: Custom colors\n")
p3 <- plot_multi_volcano(
  test_data,
  cluster_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                     "#FFFF33", "#A65628", "#F781BF", "#999999", "#000000"),
  up_color = "#FF6B6B",
  down_color = "#6B6BFF"
)
print(p3)
cat("✓ Test 3 passed\n\n")

# Test 4: Custom axis limits
cat("Test 4: Custom axis limits\n")
p4 <- plot_multi_volcano(
  test_data,
  y_limits = c(-5, 10),
  y_breaks = c(-5, -2, 0, 2, 5, 10)
)
print(p4)
cat("✓ Test 4 passed\n\n")

# Test 5: Error handling - missing column
cat("Test 5: Error handling\n")
tryCatch({
  plot_multi_volcano(
    test_data,
    fc_column = "nonexistent_column"
  )
  cat("✗ Test 5 failed: should have thrown an error\n")
}, error = function(e) {
  cat("✓ Test 5 passed: correctly caught error -", conditionMessage(e), "\n\n")
})

# Test 6: Error handling - invalid threshold
cat("Test 6: Invalid pval_threshold\n")
tryCatch({
  plot_multi_volcano(
    test_data,
    pval_threshold = 2  # Invalid: must be between 0 and 1
  )
  cat("✗ Test 6 failed: should have thrown an error\n")
}, error = function(e) {
  cat("✓ Test 6 passed: correctly caught error -", conditionMessage(e), "\n\n")
})

# Test 7: Save plot
cat("Test 7: Save plot\n")
ggplot2::ggsave(
  "/Users/lixiang/NutstoreFiles/03.编程相关/bioRtools/test_multi_volcano.png",
  p1,
  width = 10,
  height = 6,
  dpi = 300
)
cat("✓ Test 7 passed: plot saved to test_multi_volcano.png\n\n")

cat("=====================================\n")
cat("All tests completed successfully! ✓\n")
cat("=====================================\n")
