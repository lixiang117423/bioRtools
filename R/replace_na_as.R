#' Replace NA Values in a Data Frame
#'
#' Replaces all \code{NA} values in a data frame with a specified value.
#' A pipe-friendly shorthand for \code{replace(is.na(.), value)}.
#'
#' @param data A data frame.
#' @param value Replacement value for NA (default: 0).
#'
#' @return The data frame with NA values replaced.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' df <- data.frame(a = c(1, NA, 3), b = c(NA, 2, NA))
#' replace_na_as(df)
#' replace_na_as(df, "Absent")
replace_na_as <- function(data, value = 0) {
  data[is.na(data)] <- value
  data
}
