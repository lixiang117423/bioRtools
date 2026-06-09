#' Grouped Summary Statistics
#'
#' Computes descriptive statistics for a numeric variable, optionally grouped
#' by one or more categorical columns. Designed to follow \code{dplyr::group_by()}
#' in a pipe: \code{df \%>\% group_by(group) \%>\% summarise_stats(value)}.
#'
#' @param data A grouped or ungrouped data frame.
#' @param value Column name (unquoted) of the numeric variable to summarise.
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
#' # Ungrouped
#' summarise_stats(df, value)
#'
#' # Grouped
#' df %>% dplyr::group_by(group) %>% summarise_stats(value)
summarise_stats <- function(data, value, na.rm = TRUE) {
  value_col <- dplyr::pull(data, {{ value }})
  n_total <- length(value_col)
  n_miss  <- sum(is.na(value_col))

  data %>%
    dplyr::summarise(
      max    = max({{ value }}, na.rm = na.rm),
      min    = min({{ value }}, na.rm = na.rm),
      mean   = mean({{ value }}, na.rm = na.rm),
      median = stats::median({{ value }}, na.rm = na.rm),
      sd     = stats::sd({{ value }}, na.rm = na.rm),
      se     = stats::sd({{ value }}, na.rm = na.rm) / sqrt(sum(!is.na({{ value }}))),
      cv     = stats::sd({{ value }}, na.rm = na.rm) / mean({{ value }}, na.rm = na.rm) * 100,
      iqr    = stats::IQR({{ value }}, na.rm = na.rm),
      range  = max({{ value }}, na.rm = na.rm) - min({{ value }}, na.rm = na.rm),
      sum    = sum({{ value }}, na.rm = na.rm),
      n      = sum(!is.na({{ value }})),
      n_na   = sum(is.na({{ value }})),
      .groups = "drop"
    )
}
