#' Grouped Summary Statistics
#'
#' Computes descriptive statistics for a numeric variable, optionally grouped
#' by one or more categorical columns. Designed to follow \code{dplyr::group_by()}
#' in a pipe: \code{df \%>\% group_by(group) \%>\% summarise_stats(value)}.
#'
#' @param data A grouped or ungrouped data frame.
#' @param value Column name of the numeric variable to summarise. Accepts a
#'   bare name (e.g. \code{value}) or a string (e.g. \code{"value"}).
#'   Default: \code{"value"} (the conventional output column from
#'   \code{tidyr::pivot_longer()}).
#' @param na.rm Logical; remove NA values before computing (default: TRUE).
#'
#' @return A data frame with columns: \code{max}, \code{min}, \code{mean},
#'   \code{median}, \code{sd}, \code{se}, \code{cv}, \code{iqr}, \code{range},
#'   \code{sum}, \code{n}, \code{n_na}.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' df <- data.frame(group = rep(c("A", "B"), each = 5), value = c(1:10))
#'
#' # Ungrouped (default column name "value")
#' summarise_stats(df)
#'
#' # Grouped
#' df %>% dplyr::group_by(group) %>% summarise_stats()
#'
#' # Bare name
#' df %>% dplyr::group_by(group) %>% summarise_stats(value)
#'
#' # String
#' df %>% dplyr::group_by(group) %>% summarise_stats("value")
summarise_stats <- function(data, value = "value", na.rm = TRUE) {
  value_col <- rlang::as_name(rlang::enquo(value))

  if (!value_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data frame", value_col))
  }

  vals <- data[[value_col]]
  n_total <- length(vals)
  n_miss  <- sum(is.na(vals))

  data %>%
    dplyr::summarise(
      max    = max(.data[[value_col]], na.rm = na.rm),
      min    = min(.data[[value_col]], na.rm = na.rm),
      mean   = mean(.data[[value_col]], na.rm = na.rm),
      median = stats::median(.data[[value_col]], na.rm = na.rm),
      sd     = stats::sd(.data[[value_col]], na.rm = na.rm),
      se     = stats::sd(.data[[value_col]], na.rm = na.rm) / sqrt(sum(!is.na(.data[[value_col]]))),
      cv     = stats::sd(.data[[value_col]], na.rm = na.rm) / mean(.data[[value_col]], na.rm = na.rm) * 100,
      iqr    = stats::IQR(.data[[value_col]], na.rm = na.rm),
      range  = max(.data[[value_col]], na.rm = na.rm) - min(.data[[value_col]], na.rm = na.rm),
      sum    = sum(.data[[value_col]], na.rm = na.rm),
      n      = sum(!is.na(.data[[value_col]])),
      n_na   = sum(is.na(.data[[value_col]])),
      .groups = "drop"
    )
}
