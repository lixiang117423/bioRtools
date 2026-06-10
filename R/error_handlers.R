#' Error and Warning Handler Functions
#'
#' Centralized error and warning messages for consistent user feedback across bioRtools.
#' These functions ensure consistent error messaging and improve code maintainability.
#'
#' @keywords internal
#' @noRd

# ============================================================================
# Validation Errors
# ============================================================================

#' @keywords internal
#' @noRd
err_invalid_input <- function(param_name, expected) {
  sprintf(
    "Invalid '%s' argument. Expected %s.",
    param_name, expected
  )
}

#' @keywords internal
#' @noRd
err_missing_required <- function(param_name) {
  sprintf(
    "Required parameter '%s' is missing.",
    param_name
  )
}

#' @keywords internal
#' @noRd
err_dimension_mismatch <- function(dim1_name, dim2_name, dim1, dim2) {
  sprintf(
    "Dimension mismatch: %s has %d elements, but %s has %d elements.",
    dim1_name, dim1, dim2_name, dim2
  )
}

#' @keywords internal
#' @noRd
err_row_col_mismatch <- function(matrix_name, data_name) {
  sprintf(
    "Row names of %s must match row names (sample IDs) in %s.",
    matrix_name, data_name
  )
}

# ============================================================================
# Data Quality Errors
# ============================================================================

#' @keywords internal
#' @noRd
err_non_numeric <- function(column_name) {
  sprintf(
    "Column '%s' must contain numeric values only.",
    column_name
  )
}

#' @keywords internal
#' @noRd
err_non_integer_counts <- function(matrix_name) {
  sprintf(
    "%s must contain non-negative integer counts. Raw counts are required, not normalized values (FPKM, TPM, etc.).",
    matrix_name
  )
}

#' @keywords internal
#' @noRd
err_negative_values <- function(variable_name, min_value) {
  sprintf(
    "%s contains negative values (minimum: %.2f). All values must be non-negative.",
    variable_name, min_value
  )
}

#' @keywords internal
#' @noRd
err_missing_values <- function(location, count) {
  sprintf(
    "Found %d missing/NaN values in %s. Please handle missing values before analysis.",
    count, location
  )
}

#' @keywords internal
#' @noRd
err_infinite_values <- function(genes, samples) {
  sprintf(
    "Found infinite values (Inf/-Inf) in %d genes and %d samples:\n%s\nConsider removing these genes or imputing values.",
    length(unique(genes)), length(unique(samples)),
    paste(sprintf("  Gene %s in sample(s): %s", unique(genes), 
                  sapply(unique(genes), function(g) paste(samples[genes == g], collapse = ", "))), 
          collapse = "\n")
  )
}

# ============================================================================
# Package/Dependency Errors
# ============================================================================

#' @keywords internal
#' @noRd
err_missing_packages <- function(missing_pkgs, optional = FALSE) {
  type <- if (optional) "Optional" else "Required"
  sprintf(
    "%s packages not installed: %s\nInstall with: install.packages(c(%s))",
    type,
    paste(missing_pkgs, collapse = ", "),
    paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
  )
}

#' @keywords internal
#' @noRd
err_file_not_found <- function(file_path, file_type = "file") {
  sprintf(
    "%s not found at: %s\nPlease verify the path is correct.",
    file_type, file_path
  )
}

# ============================================================================
# Design Formula Errors
# ============================================================================

#' @keywords internal
#' @noRd
err_invalid_formula <- function(formula_str, reason) {
  sprintf(
    "Invalid design formula: %s\nReason: %s\nExample: ~group or ~condition + batch",
    formula_str, reason
  )
}

#' @keywords internal
#' @noRd
err_missing_design_variable <- function(var_name, available_vars) {
  sprintf(
    "Design formula references variable '%s' which is not in sample metadata.\nAvailable variables: %s",
    var_name, paste(available_vars, collapse = ", ")
  )
}

# ============================================================================
# Threshold/Parameter Errors
# ============================================================================

#' @keywords internal
#' @noRd
err_invalid_threshold <- function(param_name, value, valid_range) {
  sprintf(
    "'%s' = %.4f is invalid. Must be in range: %s",
    param_name, value, valid_range
  )
}

#' @keywords internal
#' @noRd
err_insufficient_replicates <- function(group_name, n_replicates, min_required) {
  sprintf(
    "Group '%s' has only %d replicate(s), but at least %d replicates are recommended for reliable statistics.",
    group_name, n_replicates, min_required
  )
}

# ============================================================================
# Analysis-Specific Errors
# ============================================================================

#' @keywords internal
#' @noRd
err_deseq2_workflow <- function(step, error_msg) {
  sprintf(
    "DESeq2 analysis failed at step '%s':\n%s\nCommon causes: problematic genes, insufficient replicates, or rank-deficient design.",
    step, error_msg
  )
}

#' @keywords internal
#' @noRd
err_no_significant_results <- function(analysis_type, threshold_name, threshold_value) {
  sprintf(
    "No significant %s found using %s threshold of %.4f.\nConsider relaxing thresholds or checking data quality.",
    analysis_type, threshold_name, threshold_value
  )
}

# ============================================================================
# Warning Functions
# ============================================================================

#' @keywords internal
#' @noRd
warn_data_processing <- function(message) {
  warning("bioRtools: ", message, call. = FALSE)
}

#' @keywords internal
#' @noRd
warn_genes_removed <- function(n_genes, reason) {
  warn_data_processing(
    sprintf(
      "%d genes were removed (%s).",
      n_genes, reason
    )
  )
}

#' @keywords internal
#' @noRd
warn_small_effect_sizes <- function(n_genes, threshold) {
  warn_data_processing(
    sprintf(
      "%d genes have |log2FC| < %.2f. These may represent biologically small differences despite statistical significance.",
      n_genes, threshold
    )
  )
}

#' @keywords internal
#' @noRd
warn_multiple_testing <- function(n_tests) {
  warn_data_processing(
    sprintf(
      "Multiple testing correction applied across %d tests. FDR thresholds are more stringent with larger numbers of tests.",
      n_tests
    )
  )
}

#' @keywords internal
#' @noRd
warn_imputed_values <- function(n_missing, method) {
  warn_data_processing(
    sprintf(
      "%d missing values were imputed using '%s' method.",
      n_missing, method
    )
  )
}
