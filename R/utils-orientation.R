#' Detect whether abundance data has features as rows
#'
#' Compares sample IDs from metadata against \code{data}'s row and column names.
#' Features-as-rows is signalled when sample IDs align with the columns rather
#' than the rows. Sample IDs are pooled from the explicit \code{sample_col},
#' any character column of \code{sample}, and the metadata's own row names (the
#' last covers analyses such as RDA where sample IDs live in row names).
#' Returns FALSE with no metadata or an ambiguous signal, so the historical
#' samples-as-rows layout stays the safe default.
#'
#' @param data Matrix or data frame of abundance values.
#' @param sample Sample metadata data frame, or NULL.
#' @param sample_col Optional column in \code{sample} holding sample IDs.
#' @return TRUE if \code{data} looks like features-as-rows, else FALSE.
#' @keywords internal
detect_feature_as_row <- function(data, sample, sample_col) {
  if (is.null(sample)) {
    return(FALSE)
  }

  # Pool every place sample IDs could live: an explicit column, any character
  # column, and the metadata's own row names.
  pools <- list()
  if (!is.null(sample_col) && nzchar(sample_col) &&
      sample_col %in% names(sample)) {
    pools[[sample_col]] <- sample[[sample_col]]
  }
  char_cols <- names(sample)[vapply(sample, is.character, logical(1))]
  for (col in char_cols) {
    pools[[col]] <- sample[[col]]
  }
  sample_rn <- rownames(sample)
  if (!is.null(sample_rn) && length(sample_rn) > 0L) {
    pools[["__rownames__"]] <- sample_rn
  }
  if (length(pools) == 0L) {
    return(FALSE)
  }

  count_hits <- function(ids) {
    if (is.null(ids) || length(ids) == 0L) {
      return(0L)
    }
    max(vapply(pools, function(p) length(intersect(ids, p)), integer(1)))
  }

  col_hits <- count_hits(colnames(data))
  row_hits <- count_hits(rownames(data))
  data_cn <- colnames(data)

  # feature-row when sample IDs land in data's columns, not its rows
  col_hits > row_hits &&
    !is.null(data_cn) && length(data_cn) > 0L &&
    col_hits >= ceiling(0.5 * length(data_cn))
}


#' Orient abundance data to samples-as-rows
#'
#' Resolves the requested orientation (\code{NA} auto-detects via
#' \code{\link{detect_feature_as_row}}; \code{TRUE}/\code{FALSE} force) and
#' transposes feature-row input so downstream code always sees samples x
#' features — the historical layout every analysis expects.
#'
#' @param data Matrix or data frame of abundance values.
#' @param sample Sample metadata data frame, or NULL.
#' @param sample_col Optional column in \code{sample} holding sample IDs.
#' @param feature_as_row TRUE, FALSE, or NA.
#' @param verbose Logical; print a message when transposing.
#' @return \code{data} in samples-as-rows orientation.
#' @keywords internal
orient_to_sample_row <- function(data, sample, sample_col, feature_as_row, verbose) {
  if (!is.logical(feature_as_row) || length(feature_as_row) != 1L) {
    stop("'feature_as_row' must be TRUE, FALSE, or NA")
  }
  if (is.na(feature_as_row)) {
    feature_as_row <- detect_feature_as_row(data, sample, sample_col)
  }
  if (feature_as_row) {
    if (is.data.frame(data) && !all(sapply(data, is.numeric))) {
      stop("When 'data' has features as rows, every column must be numeric")
    }
    data <- t(as.matrix(data))
    if (verbose) {
      message("Features detected as rows; transposed to samples x features.")
    }
  }
  data
}


#' Rename legacy dot-case database columns to snake_case
#'
#' For each \code{new = old} pair in \code{map}: when \code{new} is absent but
#' \code{old} is present, rename \code{old} to \code{new} and warn once (the
#' dot-case names are deprecated). Columns already in snake_case are untouched.
#'
#' @param db Data frame.
#' @param map Named character vector: names are snake_case targets, values are
#'   legacy dot-case sources.
#' @param arg_name Argument name used in the deprecation warning.
#' @return \code{db} with legacy columns renamed to snake_case.
#' @keywords internal
reconcile_db_columns <- function(db, map, arg_name) {
  renamed <- character(0)
  for (new in names(map)) {
    old <- map[[new]]
    if (!new %in% names(db) && old %in% names(db)) {
      names(db)[names(db) == old] <- new
      renamed <- c(renamed, paste0(old, " -> ", new))
    }
  }
  if (length(renamed)) {
    warning(arg_name, ": dot-case column names are deprecated; rename ",
            paste(renamed, collapse = ", "), " to snake_case.")
  }
  db
}
