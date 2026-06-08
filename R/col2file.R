#' Write a Data Frame Column to a Plain Text File
#'
#' Extracts a single column from a data frame and writes it to a text file
#' with no header and no quoting. Useful for exporting shell scripts,
#' sample lists, or any single-column data.
#'
#' @param data A data frame.
#' @param col Column name (string) or index to export.
#' @param file Output file path.
#' @param append Logical; if TRUE, append to existing file (default: TRUE).
#' @param eol End-of-line character (default: "\\n").
#'
#' @return The input \code{data}, invisibly. This makes the function
#'   pipe-friendly.
#'
#' @export
#'
#' @examples
#' df <- data.frame(cmd = c("bwa mem ref.fa R1.fq R2.fq > sample1.sam",
#'                          "bwa mem ref.fa R1.fq R2.fq > sample2.sam"))
#' tf <- tempfile(fileext = ".sh")
#' col2file(df, "cmd", tf)
#' readLines(tf)
#' unlink(tf)
col2file <- function(data, col, file, append = TRUE, eol = "\n") {
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (is.character(col)) {
    if (!col %in% names(data)) stop(sprintf("Column '%s' not found", col))
  } else if (is.numeric(col)) {
    if (col < 1 || col > ncol(data)) stop(sprintf("Column index %d out of range", col))
  } else {
    stop("'col' must be a column name (string) or index (integer)")
  }

  readr::write_delim(
    data[, col, drop = FALSE],
    file       = file,
    col_names  = FALSE,
    quote      = "none",
    append     = append,
    eol        = eol
  )
  invisible(data)
}
