#' Remove duplicate rows based on a column
#'
#' Convenience wrapper around \code{dplyr::distinct(col, .keep_all = TRUE)}.
#' Keeps the first row for each unique value of \code{col}, preserving all other
#' columns.
#'
#' @param data A data frame.
#' @param col <\code{\link[dplyr:dplyr_tidy_select]{tidy-select}}> Bare column
#'   name whose values define uniqueness.
#'
#' @return \code{data} with duplicate rows removed: only the first row for each
#'   unique value of \code{col} is kept, with all columns preserved.
#'
#' @export
#' @author Xiang LI <lixiang117423@gmail.com>
#'
#' @examples
#' \dontrun{
#' df <- data.frame(sample = c("A", "A", "B", "B", "C"), value = 1:5)
#'
#' # keep one row per sample (first occurrence)
#' dedup_by_col(df, sample)
#' }
dedup_by_col <- function(data, col) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  dplyr::distinct(data, {{ col }}, .keep_all = TRUE)
}
