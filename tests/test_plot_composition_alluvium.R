# Test script for plot_composition_alluvium function
# This script tests the new alluvium plot function

# Source the function directly (for testing before package installation)
source("R/plot_composition_alluvium.R")

library(dplyr)
library(readxl)

# Read test data
test_data <- read_excel(
  "~/NutstoreFiles/03.编程相关/bioRtools/PlotCase/物种组成冲积图/NCOMMS-24-29043C_Source Data to Main Figures.xlsx",
  sheet = "Fig.1a",
  skip = 1
)

cat("Data loaded successfully!\n")
cat("Dimensions:", nrow(test_data), "rows,", ncol(test_data), "columns\n")
cat("Columns:", paste(names(test_data), collapse = ", "), "\n\n")

# Test 1: Basic usage
cat("Test 1: Basic usage\n")
p1 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)"
)
print(p1)
cat("✓ Test 1 passed\n\n")

# Test 2: Custom color palette
cat("Test 2: Custom color palette\n")
p2 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  palette = "Set2"
)
print(p2)
cat("✓ Test 2 passed\n\n")

# Test 3: Custom alpha
cat("Test 3: Custom alpha\n")
p3 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  alpha = 0.9
)
print(p3)
cat("✓ Test 3 passed\n\n")

# Test 4: Custom y-axis formatting
cat("Test 4: Custom y-axis formatting\n")
p4 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  y_accuracy = 0.1
)
print(p4)
cat("✓ Test 4 passed\n\n")

# Test 5: Custom legend position
cat("Test 5: Custom legend position\n")
p5 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  legend_position = "bottom"
)
print(p5)
cat("✓ Test 5 passed\n\n")

# Test 6: Using id_column parameter
cat("Test 6: Using id_column parameter\n")
p6 <- plot_composition_alluvium(
  test_data,
  taxon_column = "family",
  id_column = "Viral Taxonomy (family)"
)
print(p6)
cat("✓ Test 6 passed\n\n")

# Test 7: Save plot
cat("Test 7: Save plot\n")
ggplot2::ggsave(
  "/Users/lixiang/NutstoreFiles/03.编程相关/bioRtools/test_composition_alluvium.png",
  p1,
  width = 8,
  height = 6,
  dpi = 300
)
cat("✓ Test 7 passed: plot saved to test_composition_alluvium.png\n\n")

cat("=====================================\n")
cat("All tests completed successfully! ✓\n")
cat("=====================================\n")
