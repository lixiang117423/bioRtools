#' Convert microbiome_net Result to ggNetView tbl_graph
#'
#' Converts the output of \code{\link{microbiome_net}} for a specified group
#' into a \code{tbl_graph} object compatible with \code{ggNetView::ggNetView()}.
#'
#' @param net_result Output from \code{\link{microbiome_net}}.
#' @param group Character string specifying which group to convert.
#'
#' @return A \code{tidygraph::tbl_graph} object with node attributes
#'   (name, Modularity, Degree, Strength) and edge attributes
#'   (correlation, weight, corr_direction).
#' @export
#' @seealso \code{\link{microbiome_net}}
#'
#' @examples
#' \dontrun{
#' res <- microbiome_net(df.asv, df.sample, groupCol = "group")
#' tg <- net2ggnetview(res, group = "T1")
#' ggNetView::ggNetView(tg, layout = "gephi", add_outer = TRUE, label = TRUE)
#' }
net2ggnetview <- function(net_result, group) {
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("Package 'tidygraph' is required. Install with: install.packages('tidygraph')")
  }

  if (!group %in% names(net_result$networks)) {
    stop(sprintf("Group '%s' not found. Available: %s", group,
                 paste(names(net_result$networks), collapse = ", ")))
  }

  g <- net_result$networks[[group]]
  adja <- net_result$adjaMats[[group]]
  props <- net_result$nodeProps[net_result$nodeProps$group == group, ]

  # Node attributes
  nodes <- data.frame(
    name = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Strength = igraph::strength(g),
    stringsAsFactors = FALSE
  )

  # Cluster membership
  membership <- props$cluster[match(nodes$name, props$feature)]
  top_n <- 15
  uniq_mods <- sort(unique(membership[!is.na(membership)]))
  membership_char <- as.character(membership)
  membership_char[!membership_char %in% as.character(uniq_mods[seq_len(min(top_n, length(uniq_mods)))])] <- "Others"
  nodes$Modularity <- factor(membership_char)
  nodes$modularity <- membership
  nodes$modularity2 <- nodes$Modularity
  nodes$modularity3 <- as.character(nodes$Modularity)

  # Edge attributes from igraph directly (weight + sign)
  edges_raw <- igraph::as_data_frame(g, what = "edges")
  if (nrow(edges_raw) > 0) {
    # Use edge sign attribute stored in microbiome_net
    corr <- edges_raw$weight * edges_raw$sign
    edges <- data.frame(
      from = edges_raw$from,
      to = edges_raw$to,
      correlation = corr,
      weight = abs(corr),
      corr_direction = ifelse(corr > 0, "Positive", "Negative"),
      stringsAsFactors = FALSE
    )
  } else {
    edges <- data.frame(
      from = integer(0), to = integer(0),
      correlation = numeric(0), weight = numeric(0),
      corr_direction = character(0),
      stringsAsFactors = FALSE
    )
  }

  tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
}
