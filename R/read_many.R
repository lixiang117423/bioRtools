#' Read Many Delimited Files into One Data Frame
#'
#' Reads a vector of delimited text files with \code{\link[readr]{read_delim}}
#' and row-binds them into a single data frame, adding an identifier column
#' that records which file each row came from. This is the standard
#' \code{map + read_delim + unnest} pattern collapsed into one call.
#'
#' All files must have compatible column layouts (same column types in the
#' same order); a mismatch raises the usual \code{dplyr::bind_rows} error.
#' Column-type messages from \code{readr} are suppressed.
#'
#' @param paths Character vector of file paths. Every path must exist.
#' @param id_col Name of the identifier column added to the output
#'   (default: \code{"file"}).
#' @param id_from What the identifier column stores: the full path
#'   (\code{"path"}) or the file name only (\code{"basename"}, default).
#' @param delim Single-character field delimiter (default: \code{"\t"}).
#' @param ... Additional arguments passed to \code{\link[readr]{read_delim}}
#'   (e.g. \code{col_names}, \code{col_types}, \code{skip}). A
#'   \code{show_col_types} value passed here is honored.
#'
#' @return A data frame with the \code{id_col} identifier column followed by
#'   the columns of each file.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{read_data}} for reading a single file by extension
#'
#' @examples
#' \dontrun{
#' # Count matrices split by chromosome share the same layout
#' counts <- read_many(
#'   list.files("data/counts", "\\.txt$", full.names = TRUE),
#'   id_col = "sample",
#'   col_names = c("gene", "count")
#' )
#'
#' # Combine files and derive a cleaner sample label from the path
#' counts %>%
#'   dplyr::mutate(sample = split_pick(file, index = 2))
#' }
read_many <- function(paths, id_col = "file",
                      id_from = c("basename", "path"), delim = "\t", ...) {
  if (!is.character(paths) || length(paths) == 0L || any(!nzchar(paths))) {
    stop("'paths' must be a non-empty character vector of file paths")
  }
  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("'id_col' must be a single non-empty string")
  }
  id_from <- match.arg(id_from)
  if (!is.character(delim) || length(delim) != 1L || !nzchar(delim)) {
    stop("'delim' must be a single non-empty string")
  }
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0) {
    stop("File(s) not found: ", paste(head(missing_files, 5), collapse = ", "))
  }

  dots <- list(...)
  if (is.null(dots$show_col_types)) {
    dots$show_col_types <- FALSE
  }

  ids <- if (id_from == "basename") basename(paths) else paths

  out <- purrr::map_dfr(
    stats::setNames(paths, ids),
    function(p) do.call(readr::read_delim, c(list(p, delim = delim), dots)),
    .id = id_col
  )

  message("Read ", length(paths), " file(s): ",
          nrow(out), " rows x ", ncol(out), " columns")
  out
}
