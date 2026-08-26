# Test script for plot_composition_alluvium function
# This script tests plot_composition_alluvium with synthetic data, so it runs
# in R CMD check without external files.

suppressPackageStartupMessages(library(bioRtools))

# Synthetic wide-format composition data (taxon column + sample columns)
set.seed(1)
test_data <- data.frame(
  "Viral Taxonomy (family)" = paste0("Family_", 1:18),
  Sample1 = sample(0:100, 18, replace = TRUE),
  Sample2 = sample(0:100, 18, replace = TRUE),
  Sample3 = sample(0:100, 18, replace = TRUE),
  Sample4 = sample(0:100, 18, replace = TRUE),
  check.names = FALSE
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
stopifnot(inherits(p1, "ggplot"))
cat("✓ Test 1 passed\n\n")

# Test 2: Custom color palette
cat("Test 2: Custom color palette\n")
p2 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  palette = "Set2"
)
stopifnot(inherits(p2, "ggplot"))
cat("✓ Test 2 passed\n\n")

# Test 3: Custom alpha
cat("Test 3: Custom alpha\n")
p3 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  alpha = 0.9
)
stopifnot(inherits(p3, "ggplot"))
cat("✓ Test 3 passed\n\n")

# Test 4: Custom y-axis formatting
cat("Test 4: Custom y-axis formatting\n")
p4 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  y_accuracy = 0.1
)
stopifnot(inherits(p4, "ggplot"))
cat("✓ Test 4 passed\n\n")

# Test 5: Custom legend position
cat("Test 5: Custom legend position\n")
p5 <- plot_composition_alluvium(
  test_data,
  taxon_column = "Viral Taxonomy (family)",
  legend_position = "bottom"
)
stopifnot(inherits(p5, "ggplot"))
cat("✓ Test 5 passed\n\n")

# Test 6: Using id_column parameter
cat("Test 6: Using id_column parameter\n")
p6 <- plot_composition_alluvium(
  test_data,
  taxon_column = "family",
  id_column = "Viral Taxonomy (family)"
)
stopifnot(inherits(p6, "ggplot"))
cat("✓ Test 6 passed\n\n")

# Test 7: Save plot
cat("Test 7: Save plot\n")
out_file <- tempfile(fileext = ".png")
ggplot2::ggsave(out_file, p1, width = 8, height = 6, dpi = 100)
stopifnot(file.exists(out_file))
cat("✓ Test 7 passed: plot saved to", out_file, "\n\n")

cat("=====================================\n")
cat("All tests completed successfully! ✓\n")
cat("=====================================\n")
