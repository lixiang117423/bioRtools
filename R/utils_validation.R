#' Data Validation Utility Functions
#'
#' Helper functions for common validation checks used across bioRtools functions.
#' These utilities improve code consistency and reduce duplication.
#'
#' @keywords internal
#' @noRd

# ============================================================================
# Matrix/Data Frame Validation
# ============================================================================

#' Validate count matrix format and content
#'
#' @param data Matrix or data frame of counts
#' @param param_name Parameter name for error messages (default: "data")
#' @keywords internal
#' @noRd
validate_count_matrix <- function(data, param_name = "data") {
  # Check type
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(err_invalid_input(param_name, "a matrix or data frame"))
  }

  # Check dimensions
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop(err_missing_required(param_name))
  }

  invisible(data)
}

#' Validate sample metadata format
#'
#' @param sample Data frame of sample metadata
#' @param param_name Parameter name for error messages (default: "sample")
#' @keywords internal
#' @noRd
validate_sample_metadata <- function(sample, param_name = "sample") {
  if (!is.data.frame(sample)) {
    stop(err_invalid_input(param_name, "a data frame"))
  }

  if (nrow(sample) == 0) {
    stop(err_missing_required(param_name))
  }

  invisible(sample)
}

#' Check if data and sample metadata dimensions match
#'
#' @param data Count matrix
#' @param sample Sample metadata data frame
#' @param match_by "columns" (samples in columns of data) or "rows" (samples as rows in data)
#'
#' @return Invisibly returns TRUE if valid, stops with error if not
#' @keywords internal
#' @noRd
validate_data_sample_match <- function(data, sample, match_by = "columns") {
  data_samples <- if (match_by == "columns") ncol(data) else nrow(data)
  sample_samples <- nrow(sample)

  if (data_samples != sample_samples) {
    stop(err_dimension_mismatch(
      "data", "sample",
      data_samples, sample_samples
    ))
  }

  # Check row/column name alignment if available
  data_names <- if (match_by == "columns") colnames(data) else rownames(data)
  sample_names <- rownames(sample)

  if (!is.null(data_names) && !is.null(sample_names)) {
    if (!all(data_names == sample_names)) {
      stop(err_row_col_mismatch("data", "sample metadata"))
    }
  }

  invisible(TRUE)
}

# ============================================================================
# Numeric Data Validation
# ============================================================================

#' Check for valid numeric values in matrix
#'
#' @param data Matrix/data frame
#' @param allow_negative Allow negative values (default: FALSE)
#' @param allow_zero Allow zero values (default: TRUE)
#' @param allow_na Allow NA values (default: FALSE)
#' @param param_name Parameter name for error messages
#'
#' @return Invisibly returns TRUE if valid, stops with error if not
#' @keywords internal
#' @noRd
validate_numeric_data <- function(data, allow_negative = FALSE, allow_zero = TRUE,
                                   allow_na = FALSE, param_name = "data") {
  # Check if numeric
  if (!all(sapply(data, is.numeric))) {
    stop(err_non_numeric(param_name))
  }

  # Check for NA values
  if (!allow_na && anyNA(data)) {
    n_missing <- sum(is.na(data))
    stop(err_missing_values(param_name, n_missing))
  }

  # Check for negative values
  if (!allow_negative && any(data < 0, na.rm = TRUE)) {
    min_val <- min(data, na.rm = TRUE)
    stop(err_negative_values(param_name, min_val))
  }

  # Check for infinite values
  if (any(!is.finite(data), na.rm = TRUE)) {
    n_inf <- sum(!is.finite(data), na.rm = TRUE)
    stop(sprintf("Found %d infinite values in %s", n_inf, param_name))
  }

  invisible(TRUE)
}

#' Check for integer count data
#'
#' @param data Matrix/data frame of counts
#' @param param_name Parameter name for error messages
#'
#' @return Invisibly returns TRUE if valid, stops with error if not
#' @keywords internal
#' @noRd
validate_integer_counts <- function(data, param_name = "data") {
  validate_numeric_data(data, allow_negative = FALSE, allow_zero = TRUE, param_name = param_name)

  # Check if values are integers (or very close to integers)
  data_matrix <- as.matrix(data)
  is_integer_like <- all(abs(data_matrix - round(data_matrix)) < 1e-10, na.rm = TRUE)

  if (!is_integer_like) {
    stop(err_non_integer_counts(param_name))
  }

  invisible(TRUE)
}

# ============================================================================
# Design Formula Validation
# ============================================================================

#' Validate design formula references
#'
#' @param formula Design formula
#' @param sample Sample metadata data frame
#'
#' @return Invisibly returns TRUE if valid, stops with error if not
#' @keywords internal
#' @noRd
validate_formula_vars <- function(formula, sample) {
  # Extract variable names from formula (excluding ~)
  formula_vars <- all.vars(formula)

  # Get available variables
  available_vars <- colnames(sample)

  # Check if all formula variables are in sample
  missing_vars <- setdiff(formula_vars, available_vars)

  if (length(missing_vars) > 0) {
    stop(err_missing_design_variable(missing_vars[1], available_vars))
  }

  invisible(TRUE)
}

# ============================================================================
# Threshold/Parameter Validation
# ============================================================================

#' Validate numeric threshold parameter
#'
#' @param value Numeric value to check
#' @param min Minimum allowed value
#' @param max Maximum allowed value
#' @param param_name Parameter name for error messages
#'
#' @return Invisibly returns TRUE if valid, stops with error if not
#' @keywords internal
#' @noRd
validate_threshold <- function(value, min = -Inf, max = Inf, param_name = "threshold") {
  if (value < min || value > max) {
    range_str <- sprintf("[%.4f, %.4f]", min, max)
    stop(err_invalid_threshold(param_name, value, range_str))
  }

  invisible(TRUE)
}

#' Validate logical parameter
#'
#' @param value Logical value to check
#' @param param_name Parameter name for error messages
#'
#' @return Invisibly returns TRUE if valid, stops with error if not
#' @keywords internal
#' @noRd
validate_logical <- function(value, param_name = "parameter") {
  if (!is.logical(value) || length(value) != 1) {
    stop(err_invalid_input(param_name, "a single logical value (TRUE or FALSE)"))
  }

  invisible(TRUE)
}

# ============================================================================
# Sample Size Validation
# ============================================================================

#' Check minimum replicates per group
#'
#' @param sample Sample metadata
#' @param group_col Column name containing group assignments
#' @param min_replicates Minimum required replicates per group
#' @param warn_only If TRUE, issue warning instead of error (default: TRUE)
#'
#' @return Invisibly returns TRUE, issues warning/error if insufficient replicates
#' @keywords internal
#' @noRd
validate_replicates <- function(sample, group_col = "group", min_replicates = 3, warn_only = TRUE) {
  if (!group_col %in% colnames(sample)) {
    stop(err_missing_design_variable(group_col, colnames(sample)))
  }

  group_counts <- table(sample[[group_col]])
  insufficient <- group_counts < min_replicates

  if (any(insufficient)) {
    problem_groups <- names(group_counts)[insufficient]
    for (group in problem_groups) {
      n_rep <- as.numeric(group_counts[group])
      msg <- err_insufficient_replicates(group, n_rep, min_replicates)

      if (warn_only) {
        warn_data_processing(msg)
      } else {
        stop(msg)
      }
    }
  }

  invisible(TRUE)
}

# ============================================================================
# Package Dependency Validation
# ============================================================================

#' Check if required packages are installed
#'
#' @param packages Character vector of package names to check
#' @param stop_on_missing If TRUE (default), stop with error if packages missing
#'
#' @return Invisibly returns TRUE if all present, stops/warns if missing
#' @keywords internal
#' @noRd
check_required_packages <- function(packages, stop_on_missing = TRUE) {
  installed <- sapply(packages, function(pkg) requireNamespace(pkg, quietly = TRUE))
  missing_pkgs <- packages[!installed]

  if (length(missing_pkgs) > 0) {
    if (stop_on_missing) {
      stop(err_missing_packages(missing_pkgs, optional = FALSE))
    } else {
      warn_data_processing(err_missing_packages(missing_pkgs, optional = TRUE))
    }
  }

  invisible(TRUE)
}

# ============================================================================
# File Validation
# ============================================================================

#' Check if file exists
#'
#' @param file_path Path to file
#' @param file_type Type of file for error message (default: "File")
#'
#' @return Invisibly returns TRUE if exists, stops with error if not
#' @keywords internal
#' @noRd
validate_file_exists <- function(file_path, file_type = "File") {
  if (!file.exists(file_path)) {
    stop(err_file_not_found(file_path, file_type))
  }

  invisible(TRUE)
}
