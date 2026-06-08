#' Export microbiome_net Result to Gephi Format
#'
#' Exports nodes and edges tables as CSV files that can be directly imported
#' into Gephi for interactive network visualization.
#'
#' @param net_result Output from \code{\link{microbiome_net}}.
#' @param group Character string specifying which group to export.
#' @param prefix File path prefix. "-nodes.csv" and "-edges.csv" will be appended.
#'   For example, "output/network_T1" produces "output/network_T1-nodes.csv" and
#'   "output/network_T1-edges.csv".
#'
#' @return Invisibly returns a named list with \code{nodes} and \code{edges}
#'   data frames.
#' @export
#' @seealso \code{\link{microbiome_net}}
#'
#' @examples
#' \dontrun{
#' res <- microbiome_net(df.asv, df.sample, groupCol = "group")
#' net2gephi(res, group = "T1", prefix = "result/network_T1")
#' }
net2gephi <- function(net_result, group, prefix = "network") {
  if (!group %in% names(net_result$networks)) {
    stop(sprintf("Group '%s' not found. Available: %s", group,
                 paste(names(net_result$networks), collapse = ", ")))
  }

  g <- net_result$networks[[group]]
  props <- net_result$nodeProps[net_result$nodeProps$group == group, ]

  # --- Nodes table ---
  nodes <- data.frame(
    Id       = seq_along(igraph::V(g)$name),
    Label    = igraph::V(g)$name,
    Degree   = igraph::degree(g),
    Strength = igraph::strength(g),
    stringsAsFactors = FALSE
  )

  # Add cluster from nodeProps
  cluster <- props$cluster[match(nodes$Label, props$feature)]
  nodes$Modularity <- cluster

  # Add hub status
  nodes$IsHub <- props$is_hub[match(nodes$Label, props$feature)]

  # Add eigenvector centrality
  nodes$Eigenvector <- round(props$eigenvector[match(nodes$Label, props$feature)], 6)

  # Add taxonomy if available
  tax_cols <- setdiff(names(props), c("feature", "group", "degree", "betweenness",
                     "closeness", "eigenvector", "is_hub", "cluster"))
  if (length(tax_cols) > 0) {
    tax_data <- props[match(nodes$Label, props$feature), tax_cols, drop = FALSE]
    nodes <- cbind(nodes, tax_data)
  }

  # --- Edges table ---
  edges_raw <- igraph::as_data_frame(g, what = "edges")

  edges <- data.frame(
    Source = edges_raw$from,
    Target = edges_raw$to,
    Weight = round(edges_raw$weight, 6),
    Type   = "undirected",
    stringsAsFactors = FALSE
  )

  # Add correlation sign if available
  if ("sign" %in% names(edges_raw)) {
    edges$Correlation = round(edges_raw$weight * edges_raw$sign, 6)
    edges$Direction = ifelse(edges_raw$sign > 0, "Positive", "Negative")
  }

  # --- Write files ---
  node_file <- paste0(prefix, "-nodes.csv")
  edge_file <- paste0(prefix, "-edges.csv")

  readr::write_csv(nodes, node_file)
  readr::write_csv(edges, edge_file)

  message(sprintf("Gephi files exported:\n  Nodes: %s (%d nodes)\n  Edges: %s (%d edges)",
                  node_file, nrow(nodes), edge_file, nrow(edges)))

  invisible(list(nodes = nodes, edges = edges))
}
