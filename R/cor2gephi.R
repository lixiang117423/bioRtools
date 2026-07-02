#' Export a cor_analysis Result to Gephi Format
#'
#' Converts the edge-list data frame returned by \code{\link{cor_analysis}} into
#' Gephi nodes and edges CSV files. Within-dataset correlations produce
#' reciprocal edges (A-B and B-A); these are de-duplicated so each unordered
#' pair appears once. Optionally enriches the nodes table with eigenvector
#' centrality, community (modularity) membership, and a per-community hub flag.
#'
#' @param cor_result Data frame returned by \code{\link{cor_analysis}}; must
#'   contain columns \code{from}, \code{to}, \code{cor}.
#' @param prefix File path prefix. \code{"-nodes.csv"} and \code{"-edges.csv"}
#'   are appended (default: "cor_network").
#' @param weight \code{"absolute"} (default) stores \code{|cor|} in the Gephi
#'   \code{Weight} column; \code{"signed"} stores the signed correlation. The
#'   signed value is always preserved in the \code{Correlation} column.
#' @param enrich Logical (default TRUE). When TRUE, adds eigenvector centrality,
#'   modularity, and IsHub columns to the nodes table.
#' @param community Community detection method when \code{enrich = TRUE}:
#'   \code{"auto"} (default, uses Louvain), \code{"louvain"},
#'   \code{"fast_greedy"}, or \code{"none"} (omit modularity).
#'
#' @return Invisibly, a named list with \code{nodes} and \code{edges} data
#'   frames.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{cor_analysis}}, \code{\link{net2gephi}}
#'
#' @examples
#' \dontrun{
#' res <- cor_analysis(iris[, 1:4], cor = 0.6, pvalue = 0.05)
#' cor2gephi(res, prefix = "result/cor_network")
#' }
cor2gephi <- function(cor_result, prefix = "cor_network",
                      weight = c("absolute", "signed"),
                      enrich = TRUE,
                      community = c("auto", "louvain", "fast_greedy", "none")) {
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
  community <- match.arg(community)
  if (!is.logical(enrich) || length(enrich) != 1) {
    stop("'enrich' must be a single TRUE or FALSE")
  }

  edges_raw <- cor_result
  if (nrow(edges_raw) == 0) {
    stop("'cor_result' has no significant correlations to export")
  }

  # --- Deduplicate reciprocal edges (within-dataset mode) ----------------
  # Symmetric correlations yield both (A,B) and (B,A); keep one per unordered pair.
  edges_raw$.key <- paste(pmin(edges_raw$from, edges_raw$to),
                          pmax(edges_raw$from, edges_raw$to), sep = "|")
  edges_raw <- edges_raw[!duplicated(edges_raw$.key), ]
  edges_raw$.key <- NULL

  # --- Build igraph ------------------------------------------------------
  g <- igraph::graph_from_data_frame(
    edges_raw[, c("from", "to", "cor")],
    directed = FALSE
  )
  igraph::E(g)$weight <- edges_raw$cor
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                        edge.attr.comb = "first")

  # --- Edges table -------------------------------------------------------
  el <- igraph::as_data_frame(g, what = "edges")
  signed_cor <- el$weight
  weight_col <- if (weight == "absolute") abs(signed_cor) else signed_cor

  edges <- data.frame(
    Source = el$from,
    Target = el$to,
    Weight = round(weight_col, 6),
    Type = "undirected",
    Correlation = round(signed_cor, 6),
    Direction = ifelse(signed_cor > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  if ("pvalue" %in% names(edges_raw)) {
    # Map pvalue direction-agnostically (sorted pair key) so it survives any
    # edge reordering introduced by simplify().
    pv_lookup <- stats::setNames(
      edges_raw$pvalue,
      paste(pmin(edges_raw$from, edges_raw$to),
            pmax(edges_raw$from, edges_raw$to), sep = "|")
    )
    el_key <- paste(pmin(el$from, el$to), pmax(el$from, el$to), sep = "|")
    edges$Pvalue <- round(as.numeric(pv_lookup[el_key]), 6)
  }

  # --- Nodes table -------------------------------------------------------
  node_names <- if (igraph::vcount(g)) igraph::V(g)$name else character(0)
  nodes <- data.frame(
    Id = seq_along(node_names),
    Label = node_names,
    Degree = igraph::degree(g),
    Strength = igraph::strength(g),
    stringsAsFactors = FALSE
  )

  if (enrich) {
    eig <- tryCatch(
      as.numeric(igraph::eigen_centrality(g, weights = abs(igraph::E(g)$weight))$vector),
      error = function(e) rep(NA_real_, igraph::vcount(g))
    )
    nodes$Eigenvector <- round(eig, 6)

    if (community != "none") {
      cl <- tryCatch(
        switch(community,
          louvain = igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight)),
          fast_greedy = igraph::cluster_fast_greedy(g, weights = abs(igraph::E(g)$weight)),
          auto = igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight))
        ),
        error = function(e) NULL
      )
      # Fall back to a single community if detection fails (e.g. tiny graphs).
      membership <- if (is.null(cl)) rep(1L, igraph::vcount(g))
                    else as.integer(igraph::membership(cl))
      nodes$Modularity <- membership

      # IsHub: node with the max Strength within its community (one per community).
      nodes$IsHub <- FALSE
      for (m in unique(membership)) {
        in_mod <- which(membership == m)
        hub <- in_mod[which.max(nodes$Strength[in_mod])]
        nodes$IsHub[hub] <- TRUE
      }
    }
  }

  # --- Write files -------------------------------------------------------
  node_file <- paste0(prefix, "-nodes.csv")
  edge_file <- paste0(prefix, "-edges.csv")
  readr::write_csv(nodes, node_file)
  readr::write_csv(edges, edge_file)

  message(sprintf("Gephi files exported:\n  Nodes: %s (%d nodes)\n  Edges: %s (%d edges)",
                  node_file, nrow(nodes), edge_file, nrow(edges)))

  invisible(list(nodes = nodes, edges = edges))
}
