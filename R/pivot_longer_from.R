#' Pivot Longer From a Starting Column to the End
#'
#' Thin pipe-friendly wrapper around \code{tidyr::pivot_longer()} for the
#' recurring pattern \code{tidyr::pivot_longer(cols = start:ncol(.))}. Pass
#' the starting column index; everything from there to the last column is
#' pivoted.
#'
#' @param data A data frame.
#' @param start Integer; starting column index. Pivots columns
#'   \code{start:ncol(data)}.
#' @param ... Additional arguments passed to \code{tidyr::pivot_longer()}
#'   (e.g., \code{names_to}, \code{values_to}, \code{names_prefix}).
#'
#' @return A tibble in long format.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(id = 1, type = "A", x = 1, y = 2, z = 3)
#'
#' # Instead of: df %>% tidyr::pivot_longer(cols = 3:ncol(.))
#' df %>% pivot_longer_from(3)
#'
#' # With extra args forwarded to tidyr::pivot_longer
#' df %>%
#'   pivot_longer_from(3, names_to = "variable", values_to = "value")
#' }
#'
pivot_longer_from <- function(data, start, ...) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.numeric(start) || length(start) != 1) {
    stop("'start' must be a single integer")
  }
  if (start != round(start)) {
    stop("'start' must be an integer column index")
  }
  if (start < 1 || start > ncol(data)) {
    stop(sprintf("'start' must be between 1 and %d (ncol(data))", ncol(data)))
  }

  tidyr::pivot_longer(data, cols = start:ncol(data), ...)
}
