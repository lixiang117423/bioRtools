# Tests for bioRtools core functions
# Uses testthat framework for unit testing

# Load required libraries
library(testthat)
library(bioRtools)

# ============================================================================
# Tests for error_handlers.R
# ============================================================================

test_that("error handlers generate appropriate messages", {
  # Test invalid input error
  msg <- err_invalid_input("param1", "a matrix")
  expect_true(grepl("Invalid.*param1", msg))
  expect_true(grepl("matrix", msg))

  # Test dimension mismatch
  msg <- err_dimension_mismatch("data", "sample", 100, 50)
  expect_true(grepl("Dimension mismatch", msg))
  expect_true(grepl("100", msg))
  expect_true(grepl("50", msg))

  # Test missing required
  msg <- err_missing_required("data")
  expect_true(grepl("Required.*data", msg))

  # Test non-numeric
  msg <- err_non_numeric("column1")
  expect_true(grepl("numeric", msg))
  expect_true(grepl("column1", msg))
})

test_that("warning functions work correctly", {
  # Test that warnings are issued (doesn't stop execution)
  expect_warning(
    warn_data_processing("Test warning message"),
    "bioRtools.*Test warning"
  )

  expect_warning(
    warn_genes_removed(5, "test reason"),
    "5.*genes.*removed"
  )
})

# ============================================================================
# Tests for utils_validation.R
# ============================================================================

test_that("validate_count_matrix works correctly", {
  # Valid matrix
  valid_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  expect_invisible(validate_count_matrix(valid_matrix))

  # Valid data frame
  valid_df <- data.frame(sample1 = c(1, 2), sample2 = c(3, 4))
  expect_invisible(validate_count_matrix(valid_df))

  # Invalid: not matrix or data frame
  expect_error(validate_count_matrix(c(1, 2, 3)))

  # Invalid: empty matrix
  expect_error(validate_count_matrix(matrix(nrow = 0, ncol = 0)))
})

test_that("validate_integer_counts detects non-integer data", {
  # Valid integer matrix
  integer_matrix <- matrix(c(1L, 2L, 3L, 4L), nrow = 2, ncol = 2)
  expect_invisible(validate_integer_counts(integer_matrix))

  # Invalid: non-integer values
  float_matrix <- matrix(c(1.5, 2.3, 3.1, 4.2), nrow = 2, ncol = 2)
  expect_error(validate_integer_counts(float_matrix))

  # Invalid: negative values
  negative_matrix <- matrix(c(-1L, 2L, 3L, 4L), nrow = 2, ncol = 2)
  expect_error(validate_integer_counts(negative_matrix))
})

test_that("validate_formula_vars checks formula variables", {
  # Valid formula
  sample_meta <- data.frame(group = c("A", "B"), batch = c(1, 2))
  formula <- ~ group + batch
  expect_invisible(validate_formula_vars(formula, sample_meta))

  # Invalid: missing variable
  bad_formula <- ~ group + condition
  expect_error(validate_formula_vars(bad_formula, sample_meta))
})

test_that("validate_threshold checks parameter ranges", {
  # Valid threshold
  expect_invisible(validate_threshold(0.5, min = 0, max = 1))

  # Invalid: too low
  expect_error(validate_threshold(-0.1, min = 0, max = 1))

  # Invalid: too high
  expect_error(validate_threshold(1.5, min = 0, max = 1))
})

test_that("validate_logical checks for single boolean", {
  # Valid
  expect_invisible(validate_logical(TRUE))
  expect_invisible(validate_logical(FALSE))

  # Invalid: not logical
  expect_error(validate_logical("TRUE"))

  # Invalid: vector
  expect_error(validate_logical(c(TRUE, FALSE)))
})

# ============================================================================
# Tests for color_palettes.R
# ============================================================================

test_that("pal_regulation returns correct colors", {
  colors <- pal_regulation()

  # Check all expected categories are present
  expect_true("Up-regulated" %in% names(colors))
  expect_true("Down-regulated" %in% names(colors))
  expect_true("Not significant" %in% names(colors))

  # Check colors are valid hex codes
  expect_match(colors["Up-regulated"], "^#[0-9A-F]{6}$")
  expect_match(colors["Down-regulated"], "^#[0-9A-F]{6}$")
  expect_match(colors["Not significant"], "^#[0-9A-F]{6}$")

  # Different up/down colors
  expect_false(colors["Up-regulated"] == colors["Down-regulated"])
})

test_that("pal_groups returns correct number of colors", {
  # Request 5 colors
  colors <- pal_groups(5)
  expect_length(colors, 5)

  # All colors should be valid hex codes
  expect_true(all(grepl("^#[0-9A-F]{6}$", colors)))

  # Different palettes available
  colors_pastel <- pal_groups(5, palette = "pastel")
  expect_length(colors_pastel, 5)
})

test_that("pal_significance returns appropriate categories", {
  colors <- pal_significance()

  expect_true("****" %in% names(colors))  # p < 0.0001
  expect_true("***" %in% names(colors))   # p < 0.001
  expect_true("**" %in% names(colors))    # p < 0.01
  expect_true("*" %in% names(colors))     # p < 0.05
  expect_true("NS" %in% names(colors))    # Not significant
})

# ============================================================================
# Tests for utils_plotting.R
# ============================================================================

test_that("apply_bioRtools_theme applies theme to plot", {
  library(ggplot2)

  p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point()

  p_themed <- apply_bioRtools_theme(p)

  # Check that plot still has all essential layers
  expect_true(inherits(p_themed, "ggplot"))

  # Check that theme was applied (theme object should be different)
  expect_false(is.null(p_themed$theme))
})

test_that("prepare_volcano_data creates correct structure", {
  # Create sample DESeq2-like results
  results <- data.frame(
    log2FoldChange = c(2, -1.5, 0.5, -0.3, 1.2),
    padj = c(0.001, 0.01, 0.5, 0.2, 0.05),
    baseMean = c(100, 200, 50, 150, 120)
  )

  volcano_data <- prepare_volcano_data(results, fc_threshold = 1, p_threshold = 0.05)

  # Check required columns exist
  expect_true("neg_log10_p" %in% colnames(volcano_data))
  expect_true("regulation" %in% colnames(volcano_data))

  # Check regulation classification
  expect_equal(
    volcano_data[volcano_data$log2FoldChange > 1 & volcano_data$padj < 0.05, "regulation"][1],
    "Significant"
  )
})

# ============================================================================
# Tests for helpers_deseq2.R
# ============================================================================

test_that("check_library_size_ratio handles large ratios", {
  # Create count matrix with uneven library sizes
  counts <- matrix(c(10000, 100, 5000, 200), nrow = 2, ncol = 2)

  # This should issue a warning
  expect_warning(
    check_library_size_ratio(counts, warn_ratio = 50),
    NA  # No warning if ratio is acceptable
  )

  # This should warn about large ratio
  expect_warning(
    check_library_size_ratio(counts, warn_ratio = 10)
  )
})

test_that("identify_low_abundance_features finds low features", {
  # Create count matrix
  counts <- matrix(c(1, 2, 500, 600), nrow = 2, ncol = 2)
  rownames(counts) <- c("gene1", "gene2")

  low_features <- identify_low_abundance_features(counts, min_count_threshold = 10)

  # gene1 has total count = 3 (below threshold of 10)
  expect_true("gene1" %in% low_features)

  # gene2 has total count = 1100 (above threshold)
  expect_false("gene2" %in% low_features)
})

# ============================================================================
# Tests for helpers_wgcna.R
# ============================================================================

test_that("validate_wgcna_inputs checks parameters", {
  expr_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  sample_info <- data.frame(group = rep(c("A", "B"), 5))
  traits <- c("group")

  # Valid inputs
  expect_invisible(
    validate_wgcna_inputs(expr_data, sample_info, traits)
  )

  # Invalid: missing required parameter
  expect_error(
    validate_wgcna_inputs(NULL, sample_info, traits)
  )

  # Invalid: invalid tom_type
  expect_error(
    validate_wgcna_inputs(expr_data, sample_info, traits, tom_type = "invalid")
  )

  # Invalid: invalid correlation method
  expect_error(
    validate_wgcna_inputs(expr_data, sample_info, traits,
      correlation_method = "invalid")
  )
})

test_that("filter_wgcna_genes removes low variance genes", {
  # Create expression matrix with varying gene variance
  set.seed(42)
  expr <- rbind(
    low_var = c(5, 5.1, 5.2, 5.0, 5.1),  # Low variance
    high_var = c(1, 10, 2, 9, 3)           # High variance
  )

  filtered <- filter_wgcna_genes(expr, min_sd = 2, verbose = FALSE)

  # Low variance gene should be removed
  expect_false("low_var" %in% rownames(filtered))

  # High variance gene should be kept
  expect_true("high_var" %in% rownames(filtered))
})

# ============================================================================
# Summary
# ============================================================================

# Run all tests
test_results <- testthat::test_dir("tests")
print(test_results)
