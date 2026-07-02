#' Identify hub nodes in a correlation or association network
#'
#' Computes per-node centrality (degree, strength, betweenness, closeness,
#' eigenvector) and flags hub nodes by a strength-based criterion. Accepts
#' either a \code{\link{cor_analysis}} edge list or a \code{\link{microbiome_net}}
#' result.
#'
#' @param network Either a \code{cor_analysis()} edge-list data frame (columns
#'   \code{from}, \code{to}, \code{cor}) or a \code{microbiome_net()} result
#'   (a list with a \code{$networks} component).
#' @param group Character; required when \code{network} is a
#'   \code{microbiome_net()} result — selects \code{network$networks[[group]]}.
#'   Ignored for the edge-list path.
#' @param method Hub-selection criterion: \code{"top_percent"} (default),
#'   \code{"top_n"}, or \code{"threshold"}.
#' @param top_percent Fraction of nodes by strength (default 0.1 = 10%); used
#'   when \code{method = "top_percent"}. Must be in (0, 1].
#' @param top_n Integer; number of top-strength nodes; used when
#'   \code{method = "top_n"}.
#' @param strength_threshold Numeric; nodes with strength >= this are hubs; used
#'   when \code{method = "threshold"} (must be supplied).
#' @param weight \code{"absolute"} (default, strength = sum of |cor|) or
#'   \code{"signed"} (strength = sum of signed cor).
#'
#' @return A data frame with columns \code{node, degree, strength, betweenness,
#'   closeness, eigenvector, is_hub}, sorted by strength descending.
#' @details \code{betweenness}, \code{closeness}, and \code{eigenvector} use
#'   \code{|cor|} as edge weights (igraph requires non-negative weights for
#'   these). \code{strength} uses the signed or absolute correlations per
#'   \code{weight}. Transform edge weights beforehand if you need
#'   dissimilarity-based shortest paths. For \code{microbiome_net} input,
#'   centrality is recomputed on the full \code{$networks[[group]]} graph, so
#'   values may differ from \code{$node_props}: \code{microbiome_net} normalizes
#'   degree, uses a dissimilarity for betweenness, and restricts
#'   closeness/eigenvector to the largest connected component for disconnected
#'   graphs.
#' @export
#' @author Xiang LI <lixiang117423@gmail.com>
#' @seealso \code{\link{cor_analysis}}, \code{\link{microbiome_net}},
#'   \code{\link{cor2gephi}}
#'
#' @examples
#' \dontrun{
#' edges <- cor_analysis(iris[, 1:4], cor = 0.6, pvalue = 0.05)
#' hubs <- find_hub_nodes(edges, method = "top_percent", top_percent = 0.1)
#'
#' net <- microbiome_net(df.asv, df.sample, group_col = "group")
#' hubs <- find_hub_nodes(net, group = "T1", method = "top_n", top_n = 10)
#' }
find_hub_nodes <- function(network, group = NULL,
                           method = c("top_percent", "top_n", "threshold"),
                           top_percent = 0.1, top_n = 10,
                           strength_threshold = NULL,
                           weight = c("absolute", "signed")) {
  method <- match.arg(method)
  weight <- match.arg(weight)

  # --- Resolve input to an igraph object --------------------------------
  if (is.data.frame(network)) {
    g <- cor_result_to_graph(network, weight)$g
  } else if (is.list(network) && "networks" %in% names(network)) {
    if (is.null(group)) {
      stop("'group' is required when 'network' is a microbiome_net() result")
    }
    if (!group %in% names(network$networks)) {
      stop(sprintf("Group '%s' not found. Available: %s", group,
                   paste(names(network$networks), collapse = ", ")))
    }
    g <- network$networks[[group]]
    # microbiome_net stores absolute weights in E(g)$weight and sign in E(g)$sign
    if (weight == "signed" && !is.null(igraph::E(g)$sign)) {
      igraph::E(g)$weight <- igraph::E(g)$weight * igraph::E(g)$sign
    }
  } else {
    stop("'network' must be a cor_analysis() edge list or a microbiome_net() result")
  }

  if (igraph::ecount(g) == 0) {
    stop("network has no edges; cannot compute hub centrality")
  }

  # --- Centrality -------------------------------------------------------
  nodes <- igraph::V(g)$name
  cent <- data.frame(
    node = nodes,
    degree = igraph::degree(g),
    strength = round(as.numeric(igraph::strength(g)), 6),
    betweenness = round(as.numeric(igraph::betweenness(g, normalized = TRUE,
      weights = abs(igraph::E(g)$weight))), 6),
    closeness = round(as.numeric(igraph::closeness(g, normalized = TRUE,
      weights = abs(igraph::E(g)$weight))), 6),
    eigenvector = round(
      tryCatch(as.numeric(igraph::eigen_centrality(g,
        weights = abs(igraph::E(g)$weight))$vector),
        error = function(e) rep(NA_real_, igraph::vcount(g))), 6),
    is_hub = FALSE,
    stringsAsFactors = FALSE
  )

  # --- Hub selection by strength ----------------------------------------
  n <- nrow(cent)
  if (method == "top_percent") {
    if (!is.numeric(top_percent) || length(top_percent) != 1 ||
        top_percent <= 0 || top_percent > 1) {
      stop("'top_percent' must be a single number in (0, 1]")
    }
    n_hub <- ceiling(n * top_percent)
    cent$is_hub <- rank(-cent$strength, ties.method = "min") <= n_hub
  } else if (method == "top_n") {
    if (!is.numeric(top_n) || length(top_n) != 1 ||
        top_n < 1 || top_n != round(top_n)) {
      stop("'top_n' must be a single positive integer")
    }
    if (top_n > n) {
      message(sprintf("top_n (%d) > number of nodes (%d); selecting all", top_n, n))
      top_n <- n
    }
    cent$is_hub <- rank(-cent$strength, ties.method = "min") <= top_n
  } else {
    if (is.null(strength_threshold) || !is.numeric(strength_threshold) ||
        length(strength_threshold) != 1 || is.na(strength_threshold)) {
      stop('\'strength_threshold\' must be a single non-NA number when method = "threshold"')
    }
    cent$is_hub <- cent$strength >= strength_threshold
  }

  cent <- cent[order(-cent$strength), ]
  rownames(cent) <- NULL
  cent
}
