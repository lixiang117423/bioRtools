#' Split a String and Pick One Segment
#'
#' Splits each string by a separator and picks the segment at a fixed
#' position — the generic version of the
#' \code{stringr::str_split(x, sep)} followed by picking the \code{i}-th
#' element pattern. Parsing file paths is just one application; the function
#' works on any string vector (sample IDs, taxonomy strings, column names).
#'
#' The split is literal — the separator is never treated as a regular
#' expression. A string with fewer than \code{index} segments yields
#' \code{NA} for that entry.
#'
#' Optionally, the picked segment can be split again by a second separator
#' (\code{split2}) to keep a sub-segment, e.g.
#' \code{split2 = ".", index2 = 1} turns \code{"clean.fq"} into
#' \code{"clean"}.
#'
#' @param x Character vector of strings.
#' @param sep Single non-empty separator for the split (default: \code{"/"}).
#' @param index 1-based position of the segment to pick (default: 3).
#' @param split2 Optional separator for a second split of the picked segment,
#'   or \code{NULL} (default) to skip the second step.
#' @param index2 1-based position to pick after the second split
#'   (default: 1).
#'
#' @return A character vector the same length as \code{x}; \code{NA} where a
#'   split yields too few pieces.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link[stringr]{str_split_i}} for the regex-based
#'   \pkg{stringr} equivalent
#'
#' @examples
#' # Path parsing (the typical use)
#' paths <- c("rawdata/S1/clean.fq", "rawdata/S2/clean.fq")
#' split_pick(paths, index = 2)                # "S1", "S2"
#'
#' # Any string, any separator, any position
#' split_pick(c("a_b_c", "d_e_f"), sep = "_", index = 3)  # "c", "f"
#' split_pick(c("chr1:100-200", "chr2:300-400"),
#'            sep = ":", index = 1)            # "chr1", "chr2"
#'
#' # Too-short strings become NA
#' split_pick("a/b", index = 5)                # NA
#'
#' # Optional second split: keep the part before the first dot
#' split_pick(paths, split2 = ".", index2 = 1) # "clean", "clean"
#'
#' # Typical use: add a sample column derived from the file path
#' \dontrun{
#' df %>% dplyr::mutate(sample = split_pick(path, index = 2))
#' }
split_pick <- function(x, sep = "/", index = 3, split2 = NULL, index2 = 1) {
  if (!is.character(x)) {
    stop("'x' must be a character vector")
  }
  if (!is.character(sep) || length(sep) != 1L || !nzchar(sep)) {
    stop("'sep' must be a single non-empty string")
  }
  if (!is.numeric(index) || length(index) != 1L || is.na(index) ||
      index < 1 || index != floor(index)) {
    stop("'index' must be a single positive integer")
  }
  if (!is.null(split2) &&
      (!is.character(split2) || length(split2) != 1L || !nzchar(split2))) {
    stop("'split2' must be NULL or a single non-empty string")
  }
  if (!is.numeric(index2) || length(index2) != 1L || is.na(index2) ||
      index2 < 1 || index2 != floor(index2)) {
    stop("'index2' must be a single positive integer")
  }

  parts <- stringr::str_split(x, stringr::fixed(sep))

  segment <- vapply(parts, function(pieces, idx) {
    if (length(pieces) >= idx) pieces[[idx]] else NA_character_
  }, character(1), idx = index)

  if (!is.null(split2)) {
    sub_parts <- stringr::str_split(segment, stringr::fixed(split2))
    segment <- vapply(sub_parts, function(pieces, idx) {
      if (length(pieces) >= idx) pieces[[idx]] else NA_character_
    }, character(1), idx = index2)
  }

  segment
}
