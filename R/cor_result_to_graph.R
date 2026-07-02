#' Build an undirected igraph from a cor_analysis() edge list
#'
#' Internal helper shared by [cor2gephi()] and [find_hub_nodes()]. Validates
#' the edge list, deduplicates reciprocal edges (within-dataset cor_analysis
#' emits both A-B and B-A), builds an undirected igraph, and sets edge weights.
#'
#' @param cor_result Data frame returned by [cor_analysis()]; must contain
#'   columns \code{from}, \code{to}, \code{cor}.
#' @param weight \code{"absolute"} (edge weight = |cor|) or \code{"signed"}
#'   (edge weight = signed cor).
#'
#' @return Named list: \code{g} (undirected igraph with \code{E(g)$weight}
#'   set per \code{weight} and the always-signed correlation in \code{E(g)$cor};
#'   \code{E(g)$pvalue} added when present) and \code{edges} (the deduped edge
#'   data frame, for callers needing the table form).
#' @keywords internal
cor_result_to_graph <- function(cor_result, weight = c("absolute", "signed")) {
  if (!is.data.frame(cor_result)) {
    stop("'cor_result' must be a data frame (the output of cor_analysis())")
  }
  need <- c("from", "to", "cor")
  missing_cols <- setdiff(need, names(cor_result))
  if (length(missing_cols)) {
    stop("'cor_result' is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  weight <- match.arg(weight)
  if (nrow(cor_result) == 0) {
    stop("'cor_result' has no significant correlations to export")
  }

  edges <- cor_result
  # Dedup reciprocal edges (within-dataset mode emits A-B and B-A).
  edges$.key <- paste(pmin(edges$from, edges$to),
                      pmax(edges$from, edges$to), sep = "|")
  edges <- edges[!duplicated(edges$.key), ]
  edges$.key <- NULL

  g <- igraph::graph_from_data_frame(edges[, c("from", "to")], directed = FALSE)
  igraph::E(g)$weight <- if (weight == "absolute") abs(edges$cor) else edges$cor
  igraph::E(g)$cor <- edges$cor
  if ("pvalue" %in% names(edges)) igraph::E(g)$pvalue <- edges$pvalue
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                        edge.attr.comb = "first")

  list(g = g, edges = edges)
}
