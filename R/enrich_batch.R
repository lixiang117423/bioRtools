#' Run enrichment per comparison on a DEG result table
#'
#' Used by \code{enrich_kegg()} / \code{enrich_go()} when \code{gene} is a data
#' frame (e.g. the output of \code{find_degs_deseq2()}). Splits the table by its
#' comparison column, pulls the gene list for each group (dropping
#' non-significant rows by default), runs the single-comparison enrichment
#' function, and row-binds the results with a comparison column. A group that
#' errors or yields no enriched terms is skipped with a warning.
#'
#' @param degs Data frame with gene / comparison / regulation columns.
#' @param db Annotation database (\code{kegg_db} / \code{go_db}).
#' @param enrich_fun Single-comparison enrichment function
#'   (\code{enrich_kegg} or \code{enrich_go}).
#' @param comparison_col Column naming each contrast.
#' @param gene_col Column holding gene identifiers.
#' @param regulation_col Column holding regulation labels.
#' @param keep_regulation \code{NULL} to drop rows whose regulation is
#'   \code{"NS"} (default); or a character vector of regulation levels to keep.
#' @param ... Passed to \code{enrich_fun} for every comparison
#'   (e.g. \code{p_adjust}, \code{min_pathway_size}).
#' @return Combined enrichment data frame with a comparison column (0 rows if
#'   nothing was enriched).
#' @keywords internal
enrich_batch_by_comparison <- function(degs, db, enrich_fun,
                                       comparison_col = "comparison",
                                       gene_col = "gene",
                                       regulation_col = "regulation",
                                       keep_regulation = NULL,
                                       ...) {
  if (!comparison_col %in% names(degs)) {
    stop("Batch mode needs a '", comparison_col,
         "' column in the 'gene' data frame (e.g. find_degs_deseq2() output)")
  }
  if (!gene_col %in% names(degs)) {
    stop("Batch mode needs a '", gene_col, "' column in the 'gene' data frame")
  }

  # Filter by regulation when the column is present
  if (regulation_col %in% names(degs)) {
    if (is.null(keep_regulation)) {
      degs <- degs[degs[[regulation_col]] != "NS", , drop = FALSE]
    } else {
      degs <- degs[degs[[regulation_col]] %in% keep_regulation, , drop = FALSE]
    }
  } else if (!is.null(keep_regulation)) {
    stop("keep_regulation is set but column '", regulation_col,
         "' is not present in the 'gene' data frame")
  }

  if (nrow(degs) == 0) {
    warning("No rows remain after regulation filtering; nothing to enrich.")
    return(dplyr::tibble(!!comparison_col := character(0)))
  }

  comparisons <- unique(as.character(degs[[comparison_col]]))
  results <- vector("list", length(comparisons))

  for (i in seq_along(comparisons)) {
    comp <- comparisons[i]
    genes <- degs[[gene_col]][as.character(degs[[comparison_col]]) == comp]
    genes <- unique(genes[!is.na(genes) & genes != ""])
    if (length(genes) == 0) {
      warning("Comparison '", comp, "': no genes after filtering; skipped.")
      next
    }

    res <- tryCatch(
      enrich_fun(genes, db, ...),
      error = function(e) {
        warning("Comparison '", comp, "': enrichment failed - ",
                conditionMessage(e))
        NULL
      }
    )

    if (is.null(res)) next
    if (nrow(res) == 0) {
      message("Comparison '", comp, "': no enriched terms found.")
      next
    }
    message("Comparison '", comp, "': ", nrow(res), " enriched terms.")
    res[[comparison_col]] <- comp
    results[[i]] <- res
  }

  combined <- dplyr::bind_rows(results)
  if (nrow(combined) == 0) {
    return(dplyr::tibble(!!comparison_col := character(0)))
  }
  combined
}
