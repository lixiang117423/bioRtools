#' Filter rows by matching a string/regex pattern in a column
#'
#' Convenience wrapper around \code{dplyr::filter(stringr::str_detect(...))}.
#' Keeps rows where \code{col} matches \code{pattern}; with \code{negate = TRUE}
#' it keeps the rows that do NOT match (the \code{!str_detect(...)} case).
#'
#' @param data A data frame.
#' @param col <\code{\link[dplyr:dplyr_tidy_select]{tidy-select}}> Bare column
#'   name to search, e.g. \code{kegg_term}.
#' @param pattern Regex pattern to match against the column.
#' @param negate Logical; if \code{TRUE}, keep rows that do NOT match.
#'
#' @return \code{data} filtered to matching rows. Rows where \code{col} is
#'   \code{NA} are dropped, matching \code{dplyr::filter()} semantics.
#'
#' @export
#' @author Xiang LI <lixiang117423@gmail.com>
#'
#' @examples
#' \dontrun{
#' df <- data.frame(name = c("MAPK1", "MAPK3", "AKT1", "mTOR"))
#'
#' # keep rows whose name contains "MAPK"
#' filter_str(df, name, "MAPK")
#'
#' # keep rows whose name does NOT contain "MAPK"
#' filter_str(df, name, "MAPK", negate = TRUE)
#' }
filter_str <- function(data, col, pattern, negate = FALSE) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.character(pattern) || length(pattern) != 1 || is.na(pattern)) {
    stop("'pattern' must be a single non-NA character string")
  }
  if (!is.logical(negate) || length(negate) != 1 || is.na(negate)) {
    stop("'negate' must be a single TRUE or FALSE")
  }

  if (negate) {
    return(dplyr::filter(data, !stringr::str_detect({{ col }}, pattern)))
  }

  dplyr::filter(data, stringr::str_detect({{ col }}, pattern))
}
