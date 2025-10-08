#' Convert P-values to Significance Labels
#'
#' @description
#' Convert numeric p-values to significance labels using standard statistical
#' thresholds. Commonly used for adding significance annotations to plots.
#'
#' @param p_values Numeric vector of p-values. Values should be between 0 and 1.
#' @param ns Character string. Label for non-significant results (p >= 0.05).
#'   Default is "NS".
#' @param thresholds Named numeric vector of significance thresholds in ascending
#'   order. Names are the labels, values are the upper bounds (exclusive).
#'   Default is \code{c("***" = 0.001, "**" = 0.01, "*" = 0.05)}.
#'
#' @details
#' The function assigns labels based on the following logic:
#' \itemize{
#'   \item p < 0.001: "***"
#'   \item 0.001 <= p < 0.01: "**"
#'   \item 0.01 <= p < 0.05: "*"
#'   \item p >= 0.05: "NS" (or custom label)
#' }
#'
#' Custom thresholds can be provided to change these boundaries. The function
#' handles NA values by preserving them in the output.
#'
#' @return A character vector of the same length as \code{p_values} containing
#'   significance labels. NA values in input are preserved as NA in output.
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @examples
#' # Basic usage
#' p_vals <- c(0.0001, 0.005, 0.03, 0.07, 0.5, NA)
#' label_signif(p_vals)
#' # [1] "***" "**"  "*"   "NS"  "NS"  NA
#'
#' # Custom non-significant label
#' label_signif(p_vals, ns = "ns")
#' # [1] "***" "**"  "*"   "ns"  "ns"  NA
#'
#' # Change non-significant to empty string for cleaner plots
#' label_signif(p_vals, ns = "")
#' # [1] "***" "**"  "*"   ""    ""    NA
#'
#' # Custom thresholds (4-level system)
#' custom_thresh <- c(
#'   "****" = 0.0001,
#'   "***" = 0.001,
#'   "**" = 0.01,
#'   "*" = 0.05
#' )
#' label_signif(c(0.00005, 0.0005, 0.005, 0.03, 0.1), thresholds = custom_thresh)
#' # [1] "****" "***"  "**"   "*"    "NS"
#'
#' # Use in a data analysis pipeline
#' library(dplyr)
#' data.frame(
#'   gene = c("Gene1", "Gene2", "Gene3"),
#'   p_value = c(0.0001, 0.03, 0.5)
#' ) %>%
#'   mutate(significance = label_signif(p_value))
#'
#' @export
label_signif <- function(p_values,
                         ns = "NS",
                         thresholds = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) {
  # Input validation
  validate_p_values(p_values)
  validate_ns_label(ns)
  validate_thresholds(thresholds)
  
  # Sort thresholds in ascending order for proper application
  thresholds <- sort(thresholds)
  
  # Apply significance labels
  apply_labels_to_p_values(p_values, ns, thresholds)
}


#' Validate P-values
#'
#' @param p_values Numeric vector to validate
#' @keywords internal
#' @noRd
validate_p_values <- function(p_values) {
  if (!is.numeric(p_values)) {
    stop(
      "p_values must be numeric, not ",
      class(p_values)[1],
      call. = FALSE
    )
  }
  
  # Check for values outside [0, 1] range (excluding NA)
  out_of_range <- p_values[!is.na(p_values) & (p_values < 0 | p_values > 1)]
  
  if (length(out_of_range) > 0) {
    warning(
      "Some p-values are outside the valid range [0, 1]. ",
      "Found values: ",
      paste(head(out_of_range, 3), collapse = ", "),
      if (length(out_of_range) > 3) "..." else "",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Validate Non-significant Label
#'
#' @param ns Character string to validate
#' @keywords internal
#' @noRd
validate_ns_label <- function(ns) {
  if (!is.character(ns) || length(ns) != 1) {
    stop(
      "ns must be a single character string",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Validate Thresholds
#'
#' @param thresholds Named numeric vector to validate
#' @keywords internal
#' @noRd
validate_thresholds <- function(thresholds) {
  if (!is.numeric(thresholds)) {
    stop(
      "thresholds must be numeric",
      call. = FALSE
    )
  }
  
  if (is.null(names(thresholds)) || any(names(thresholds) == "")) {
    stop(
      "thresholds must be a fully named numeric vector",
      call. = FALSE
    )
  }
  
  if (any(thresholds < 0 | thresholds > 1)) {
    stop(
      "All threshold values must be between 0 and 1",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Apply Significance Labels to P-values
#'
#' @param p_values Numeric vector of p-values
#' @param ns Non-significant label
#' @param thresholds Sorted named numeric vector of thresholds
#' @return Character vector of labels
#' @keywords internal
#' @noRd
apply_labels_to_p_values <- function(p_values, ns, thresholds) {
  # Initialize all values as non-significant
  labels <- rep(ns, length(p_values))
  
  # Apply labels from most stringent to least stringent threshold
  # This ensures that the most significant label is applied
  for (i in length(thresholds):1) {
    threshold_value <- thresholds[i]
    threshold_label <- names(thresholds)[i]
    
    # Update labels for values below this threshold
    labels[!is.na(p_values) & p_values < threshold_value] <- threshold_label
  }
  
  # Preserve NA values in output
  labels[is.na(p_values)] <- NA_character_
  
  return(labels)
}


#' @rdname label_signif
#' @export
label_significance <- label_signif