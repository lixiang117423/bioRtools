#' Construct a bioRtools analysis result
#'
#' Internal constructor used by the core multivariate workflows. It preserves
#' the historical list fields while adding a common model/table/plot contract
#' and provenance metadata.
#'
#' @param x Existing result list.
#' @param model Fitted model object, or a named list of fitted models.
#' @param tables Named list of analysis-ready data frames.
#' @param plots Named list of ggplot objects.
#' @param analysis_type Human-readable analysis name.
#' @param engine Engine used to fit the model.
#' @param engine_versions Named character vector of engine package versions.
#' @param data_info List describing the validated input dimensions and layout.
#' @param call The originating call.
#' @param seed Optional random seed used by the analysis.
#'
#' @return A list with class `bioRtools_fit`.
#' @keywords internal
#' @noRd
new_bioRtools_fit <- function(x,
                              model,
                              tables,
                              plots = list(),
                              analysis_type,
                              engine,
                              engine_versions = character(),
                              data_info = list(),
                              call = NULL,
                              seed = NULL) {
  if (!is.list(x)) {
    stop("'x' must be a list", call. = FALSE)
  }
  if (!is.list(tables) || is.null(names(tables)) || any(!nzchar(names(tables)))) {
    stop("'tables' must be a named list", call. = FALSE)
  }
  if (!is.list(plots)) {
    stop("'plots' must be a list", call. = FALSE)
  }

  x$model <- model
  x$tables <- tables
  x$plots <- plots
  x$provenance <- list(
    analysis_type = analysis_type,
    engine = engine,
    engine_versions = engine_versions,
    package_version = tryCatch(
      as.character(utils::packageVersion("bioRtools")),
      error = function(e) NA_character_
    ),
    r_version = R.version.string,
    data = data_info,
    call = call,
    seed = seed
  )

  class(x) <- unique(c("bioRtools_fit", class(x)))
  x
}

#' @export
print.bioRtools_fit <- function(x, ...) {
  provenance <- x$provenance
  data_info <- provenance$data

  cat("<bioRtools_fit>\n")
  cat("Analysis:", provenance$analysis_type, "\n")
  cat("Engine:", provenance$engine, "\n")
  if (length(data_info) > 0) {
    cat("Samples:", data_info$n_samples %||% NA_integer_,
        " Features:", data_info$n_features %||% NA_integer_, "\n")
  }
  cat("Tables:", paste(names(x$tables), collapse = ", "), "\n")
  cat("Plots:", if (length(x$plots)) paste(names(x$plots), collapse = ", ") else "none", "\n")
  invisible(x)
}

#' @export
summary.bioRtools_fit <- function(object, ...) {
  table_info <- data.frame(
    name = names(object$tables),
    class = vapply(object$tables, function(x) class(x)[1], character(1)),
    rows = vapply(object$tables, function(x) if (is.data.frame(x)) nrow(x) else NA_integer_, integer(1)),
    columns = vapply(object$tables, function(x) if (is.data.frame(x)) ncol(x) else NA_integer_, integer(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      provenance = object$provenance,
      tables = table_info,
      plots = names(object$plots)
    ),
    class = "summary.bioRtools_fit"
  )
}

#' @export
print.summary.bioRtools_fit <- function(x, ...) {
  print(x$provenance)
  cat("\nTables:\n")
  print(x$tables, row.names = FALSE)
  cat("\nPlots:", if (length(x$plots)) paste(x$plots, collapse = ", ") else "none", "\n")
  invisible(x)
}

#' Convert a bioRtools result table to a tibble
#'
#' @param x A `bioRtools_fit` object.
#' @param which Name of the table to return. By default, returns the first
#'   data-frame table in the result.
#' @param ... Passed to [tibble::as_tibble()].
#'
#' @return A tibble.
#' @export
as_tibble.bioRtools_fit <- function(x, ..., which = NULL) {
  if (is.null(which)) {
    candidates <- names(x$tables)[vapply(x$tables, is.data.frame, logical(1))]
    if (length(candidates) == 0) {
      stop("This result contains no data-frame tables", call. = FALSE)
    }
    which <- candidates[1]
  }
  if (!is.character(which) || length(which) != 1 || !which %in% names(x$tables)) {
    stop("'which' must name one of: ", paste(names(x$tables), collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(x$tables[[which]])) {
    stop("The selected table is not a data frame: ", which, call. = FALSE)
  }
  tibble::as_tibble(x$tables[[which]], ...)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x
}
