#' Convert Column Names to a Single-Column Data Frame
#'
#' @description
#' Extract column names of a data frame (or matrix) and return them as a
#' single-column data frame. Useful for building sample metadata tables from
#' the column names of a count matrix in a pipe chain:
#'
#' \preformatted{
#' df \%>\%
#'   colnames_to_df("sample") \%>\%
#'   dplyr::filter(sample != "description") \%>\%
#'   dplyr::mutate(group = stringr::str_split(sample, "R") \%>\% sapply("[", 1))
#' }
#'
#' @param df A data frame or matrix.
#' @param name Character string for the output column name (default: "sample").
#'
#' @return A single-column data frame.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(T1R1 = 1:3, T1R2 = 4:6, T2R1 = 7:9)
#' colnames_to_df(df, "sample")
#'
#' # Pipe into further wrangling
#' df %>%
#'   colnames_to_df("sample") %>%
#'   dplyr::mutate(group = stringr::str_split(sample, "R") %>% sapply("[", 1))
#' }
#'
colnames_to_df <- function(df, name = "sample") {

  if (!is.data.frame(df) && !is.matrix(df)) {
    stop("'df' must be a data frame or matrix")
  }

  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("'name' must be a single non-empty character string")
  }

  stats::setNames(
    data.frame(colnames(df), stringsAsFactors = FALSE),
    name
  )
}
