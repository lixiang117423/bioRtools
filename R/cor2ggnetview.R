#' Convert a cor_analysis Result to ggNetView tbl_graph
#'
#' Converts the edge-list data frame returned by \code{\link{cor_analysis}} into
#' a \code{tidygraph::tbl_graph} object compatible with \code{ggNetView::ggNetView()}.
#' Reciprocal edges produced by within-dataset correlations (A-B and B-A) are
#' de-duplicated so each unordered pair appears once. Nodes carry degree,
#' strength, and modularity (community) membership; edges carry the signed
#' correlation, a weight column, direction, and p-value.
#'
#' @param cor_result Data frame returned by \code{\link{cor_analysis}}; must
#'   contain columns \code{from}, \code{to}, \code{cor}.
#' @param weight \code{"absolute"} (default) stores \code{|cor|} in the edge
#'   \code{weight} column; \code{"signed"} stores the signed correlation. The
#'   signed value is always preserved in the \code{correlation} column.
#' @param community Community detection method for node modules:
#'   \code{"auto"} (default, uses Louvain), \code{"louvain"},
#'   \code{"fast_greedy"}, or \code{"none"} (all nodes in a single module).
#'   Modules beyond the top 15 by size are collapsed into \code{"Others"}.
#'
#' @return A \code{tidygraph::tbl_graph} object with node attributes
#'   (name, degree, strength, modularity_factor, modularity, modularity2,
#'   modularity3, and the capitalized Modularity/Degree/Strength aliases
#'   required by \code{ggNetView}) and edge attributes (from, to, correlation,
#'   weight, corr_direction, pvalue when present in the input).
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{cor_analysis}}, \code{\link{cor2gephi}},
#'   \code{\link{net2ggnetview}}
#'
#' @examples
#' \dontrun{
#' res <- cor_analysis(iris[, 1:4], cor = 0.6, pvalue = 0.05)
#' tg <- cor2ggnetview(res)
#' ggNetView::ggNetView(tg, label = TRUE)
#' }
cor2ggnetview <- function(cor_result,
                          weight = c("absolute", "signed"),
                          community = c("auto", "louvain", "fast_greedy", "none")) {
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("Package 'tidygraph' is required. Install with: utils::install.packages('tidygraph')")
  }
  weight <- match.arg(weight)
  community <- match.arg(community)

  cg <- cor_result_to_graph(cor_result, weight)
  g <- cg$g

  # --- Node attributes ----------------------------------------------------
  nodes <- data.frame(
    name = igraph::V(g)$name,
    degree = igraph::degree(g),
    strength = igraph::strength(g, weights = abs(igraph::E(g)$weight)),
    stringsAsFactors = FALSE
  )

  membership <- switch(community,
    louvain = igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight)),
    fast_greedy = igraph::cluster_fast_greedy(g, weights = abs(igraph::E(g)$weight)),
    auto = igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight)),
    none = NULL
  )
  # Fall back to a single module when detection fails (e.g. tiny graphs)
  # or when community = "none"; keeps the node schema stable for ggNetView.
  membership <- if (is.null(membership)) {
    rep(1L, nrow(nodes))
  } else {
    tryCatch(as.integer(igraph::membership(membership)),
             error = function(e) rep(1L, nrow(nodes)))
  }

  # Keep only the top 15 modules by size; larger networks would otherwise
  # produce unreadable color legends in ggNetView.
  top_n <- 15
  uniq_mods <- sort(unique(membership[!is.na(membership)]))
  membership_char <- as.character(membership)
  membership_char[!membership_char %in% as.character(uniq_mods[seq_len(min(top_n, length(uniq_mods)))])] <- "Others"
  nodes$modularity_factor <- factor(membership_char)
  nodes$modularity <- membership
  nodes$modularity2 <- nodes$modularity_factor
  nodes$modularity3 <- as.character(nodes$modularity_factor)

  # ggNetView 0.1.0 的 layout = "gephi" 路径 arrange(modularity3, desc(Degree))、
  # 绘图用 size = Degree / group.by = "Modularity"，要求大写列（与其
  # build_graph_from_* 构建器输出一致），小写列会触发 "object 'Degree' not found"。
  nodes$Modularity <- nodes$modularity2
  nodes$Degree <- nodes$degree
  nodes$Strength <- nodes$strength

  # --- Edge attributes ----------------------------------------------------
  edges_raw <- igraph::as_data_frame(g, what = "edges")
  if (nrow(edges_raw) > 0) {
    edges <- data.frame(
      from = edges_raw$from,
      to = edges_raw$to,
      correlation = edges_raw$cor,
      weight = edges_raw$weight,
      corr_direction = ifelse(edges_raw$cor > 0, "Positive", "Negative"),
      stringsAsFactors = FALSE
    )
    if ("pvalue" %in% names(edges_raw)) edges$pvalue <- edges_raw$pvalue
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
