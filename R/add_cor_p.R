#' Add Per-Group Correlation and P-value Columns
#'
#' For each combination of grouping columns, compute the correlation
#' coefficient and p-value between two numeric columns, and add them as
#' new columns (\code{cor} and \code{p}) broadcast to all rows in the group.
#' Optionally also add linear regression \code{lm_r2} and \code{lm_pvalue}
#' from \code{lm(y ~ x)} on raw values.
#'
#' This is a pipe-friendly replacement for the verbose
#' \code{df \%>\% group_by(...) \%>\% mutate(cor = cor(x, y), p = cor.test(x, y)$p.value) \%>\% ungroup()}
#' pattern common in grouped correlation workflows.
#'
#' @param df A data frame.
#' @param x,y Value columns to correlate. Accept bare names (e.g. \code{otu.value})
#'   or strings (e.g. \code{"otu.value"}).
#' @param from,to Grouping columns identifying the two entities being
#'   correlated (e.g., \code{from = "otu"}, \code{to = "phe"}). Must be
#'   character strings (column names). When both are provided, the function
#'   groups by these two columns. When both are NULL (default), the function
#'   respects any existing grouping from a prior \code{dplyr::group_by()}.
#'   Providing only one of \code{from}/\code{to} is an error.
#' @param method Correlation method: \code{"pearson"} (default), \code{"kendall"},
#'   or \code{"spearman"}.
#' @param add_regression Logical; if TRUE, also add \code{lm_r2} and
#'   \code{lm_pvalue} from \code{lm(y ~ x)} on raw values. For
#'   \code{method = "pearson"}, \code{lm_r2} equals \code{cor^2} and
#'   \code{lm_pvalue} equals \code{p} by construction. Default FALSE.
#'
#' @return A data frame with \code{cor}, \code{p}, and \code{method} columns
#'   added. When \code{add_regression = TRUE}, also adds \code{lm_r2} and
#'   \code{lm_pvalue}. The \code{method} column records the analysis method
#'   used (e.g., \code{"pearson"} or \code{"pearson + lm"}), convenient for
#'   populating the Methods section in papers. Groups with fewer than 3
#'   non-NA observations get \code{NA} for the failed statistics.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # Explicit from/to grouping
#' df %>%
#'   add_cor_p(otu.value, phe.value, from = "otu", to = "phe")
#'
#' # Or rely on existing dplyr::group_by
#' df %>%
#'   dplyr::group_by(otu, phe) %>%
#'   add_cor_p(otu.value, phe.value)
#'
#' # Also add linear regression R^2 and p-value
#' df %>%
#'   add_cor_p(otu.value, phe.value,
#'             from = "otu", to = "phe",
#'             add_regression = TRUE)
#' }
#'
add_cor_p <- function(df, x, y, from = NULL, to = NULL,
                      method = "pearson", add_regression = FALSE) {

  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }

  x_name <- rlang::as_name(rlang::enquo(x))
  y_name <- rlang::as_name(rlang::enquo(y))

  if (!x_name %in% names(df)) {
    stop(sprintf("Column '%s' not found in data frame", x_name))
  }
  if (!y_name %in% names(df)) {
    stop(sprintf("Column '%s' not found in data frame", y_name))
  }
  if (!is.numeric(df[[x_name]])) {
    stop(sprintf("Column '%s' must be numeric", x_name))
  }
  if (!is.numeric(df[[y_name]])) {
    stop(sprintf("Column '%s' must be numeric", y_name))
  }
  if (!method %in% c("pearson", "kendall", "spearman")) {
    stop("'method' must be one of: 'pearson', 'kendall', 'spearman'")
  }
  if (!is.logical(add_regression) || length(add_regression) != 1L) {
    stop("'add_regression' must be a single TRUE or FALSE")
  }

  # from/to: both or neither
  if (xor(is.null(from), is.null(to))) {
    stop("'from' and 'to' must be provided together (both or neither)")
  }

  if (!is.null(from)) {
    if (!is.character(from) || length(from) != 1L || is.na(from)) {
      stop("'from' must be a single character string")
    }
    if (!is.character(to) || length(to) != 1L || is.na(to)) {
      stop("'to' must be a single character string")
    }
    if (from == to) {
      stop("'from' and 'to' must be different columns")
    }
    missing_cols <- setdiff(c(from, to), names(df))
    if (length(missing_cols) > 0L) {
      stop(sprintf("Group column(s) not found: %s",
                   paste(missing_cols, collapse = ", ")))
    }
    df <- df %>% dplyr::group_by(dplyr::across(dplyr::all_of(c(from, to))))
  }

  result <- df %>%
    dplyr::mutate(
      cor = stats::cor(.data[[x_name]], .data[[y_name]],
                       method = method, use = "pairwise.complete.obs"),
      p = safe_cor_test_p(.data[[x_name]], .data[[y_name]], method = method)
    )

  if (add_regression) {
    result <- result %>%
      dplyr::mutate(
        lm_r2 = safe_lm_r2(.data[[x_name]], .data[[y_name]]),
        lm_pvalue = safe_lm_pvalue(.data[[x_name]], .data[[y_name]])
      )
  }

  method_value <- if (add_regression) paste(method, "+ lm") else method

  result %>%
    dplyr::ungroup() %>%
    dplyr::mutate(method = method_value)
}

#' Compute correlation p-value, returning NA on failure
#' @keywords internal
safe_cor_test_p <- function(x, y, method) {
  tryCatch({
    stats::cor.test(x, y, method = method)$p.value
  }, error = function(e) NA_real_)
}

#' Compute lm R-squared, returning NA on failure
#' @keywords internal
safe_lm_r2 <- function(x, y) {
  tryCatch({
    summary(stats::lm(y ~ x))$r.squared
  }, error = function(e) NA_real_)
}

#' Compute lm F-test p-value, returning NA on failure
#' @keywords internal
safe_lm_pvalue <- function(x, y) {
  tryCatch({
    s <- summary(stats::lm(y ~ x))
    f_stat <- s$fstatistic
    if (is.na(f_stat[1])) NA_real_ else {
      stats::pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    }
  }, error = function(e) NA_real_)
}
