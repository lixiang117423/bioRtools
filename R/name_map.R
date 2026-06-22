#' Show Column Name Mapping for a Data Frame
#'
#' Returns a two-column data frame mapping sanitized column names
#' (\code{names(df)}) to the original column names from the file
#' (\code{attr(df, "raw_names")}, attached by \code{\link{read_data}}).
#' Useful for inspecting how a tabular file's headers were sanitized
#' (e.g., for Excel column headers like \code{"Plant height (PH)"}).
#'
#' @param df A data frame, ideally returned by \code{read_data()}.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{sanitized}: Current column name (\code{names(df)}).
#'     \item \code{raw}: Original column name from the file. If no
#'       \code{raw_names} attribute exists, falls back to \code{sanitized}.
#'   }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' df <- read_data("data.xlsx")
#' name_map(df)
#' #        sanitized             raw
#' # 1              id              id
#' # 2       treatment       treatment
#' # 3 Plant.height..PH. Plant height (PH)
#' # 4  Root.length..RL.  Root length (RL)
#' }
#'
name_map <- function(df) {

  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }

  sanitized <- names(df)
  raw <- attr(df, "raw_names")

  if (is.null(raw)) {
    warning("No 'raw_names' attribute. Did you read with bioRtools::read_data()? Falling back to sanitized as raw.")
    raw <- sanitized
  } else if (length(raw) != length(sanitized)) {
    warning(sprintf("Length mismatch: raw_names has %d entries but data frame has %d columns. Output rows may be misaligned.", length(raw), length(sanitized)))
  }

  data.frame(sanitized = sanitized, raw = raw, stringsAsFactors = FALSE)
}
