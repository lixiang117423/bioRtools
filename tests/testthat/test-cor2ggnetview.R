# cor2ggnetview(): cor_analysis edge list -> ggNetView-compatible tbl_graph.

make_cor_input <- function() {
  # Symmetric within-dataset correlations -> reciprocal edges (A-B and B-A).
  X <- data.frame(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(2, 4, 6, 8, 10, 12),   # a,b perfectly correlated
    c = c(6, 5, 4, 3, 2, 1)     # a,c perfectly anti-correlated
  )
  cor_analysis(X, method = "pearson", cor = 0.8, pvalue = 0.05)
}

test_that("cor2ggnetview validates its input", {
  expect_error(cor2ggnetview(data.frame()))
})

test_that("returns a tbl_graph with ggNetView node/edge schema", {
  tg <- cor2ggnetview(make_cor_input())
  expect_s3_class(tg, "tbl_graph")

  nodes <- as.data.frame(igraph::as_data_frame(tg, what = "vertices"))
  edges <- as.data.frame(igraph::as_data_frame(tg, what = "edges"))

  expect_true(all(c("name", "degree", "strength", "modularity_factor",
                    "modularity", "modularity2", "modularity3") %in% names(nodes)))
  expect_true(all(c("from", "to", "correlation", "weight", "corr_direction",
                    "pvalue") %in% names(edges)))
})

test_that("nodes carry capitalized Degree/Strength/Modularity for ggNetView layouts", {
  # ggNetView(layout = "gephi") internally calls arrange(modularity3, desc(Degree))
  # and maps size = Degree, group.by = "Modularity" — lowercase-only columns
  # make ggNetView fail with "object 'Degree' not found".
  tg <- cor2ggnetview(make_cor_input())
  nodes <- as.data.frame(igraph::as_data_frame(tg, what = "vertices"))

  expect_true(all(c("Modularity", "Degree", "Strength") %in% names(nodes)))
  expect_equal(nodes$Degree, nodes$degree)
  expect_equal(nodes$Strength, nodes$strength)
})

test_that("dedup removes reciprocal edges", {
  tg <- cor2ggnetview(make_cor_input())
  edges <- as.data.frame(igraph::as_data_frame(tg, what = "edges"))

  # Each unordered pair appears exactly once.
  key <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "|")
  expect_false(any(duplicated(key)))
  expect_setequal(key, c("a|b", "a|c", "b|c"))
})

test_that("edge attributes are signed cor + abs weight + direction", {
  tg <- cor2ggnetview(make_cor_input())
  edges <- as.data.frame(igraph::as_data_frame(tg, what = "edges"))

  expect_true(all(edges$weight == abs(edges$correlation)))
  expect_true(all(edges$corr_direction == ifelse(edges$correlation > 0,
                                                 "Positive", "Negative")))
  expect_true(all(edges$pvalue <= 0.05))
})

test_that("weight = signed puts signed cor in weight", {
  tg <- cor2ggnetview(make_cor_input(), weight = "signed")
  edges <- as.data.frame(igraph::as_data_frame(tg, what = "edges"))
  expect_true(all(edges$weight == edges$correlation))
})

test_that("community options control module columns", {
  inp <- make_cor_input()

  tg_auto <- cor2ggnetview(inp, community = "auto")
  nodes_auto <- as.data.frame(igraph::as_data_frame(tg_auto, what = "vertices"))
  expect_true(all(!is.na(nodes_auto$modularity)))

  tg_none <- cor2ggnetview(inp, community = "none")
  nodes_none <- as.data.frame(igraph::as_data_frame(tg_none, what = "vertices"))
  expect_length(unique(nodes_none$modularity), 1)
})

test_that("cross-dataset result (no reciprocal edges) round-trips", {
  X1 <- data.frame(a = c(1, 2, 3, 4, 5, 6), b = c(2, 4, 6, 8, 10, 12))
  X2 <- data.frame(x = c(6, 5, 4, 3, 2, 1), y = c(1, 3, 2, 5, 4, 6))
  res <- suppressMessages(cor_analysis(X1, X2, cor = 0.8, pvalue = 0.05))

  tg <- cor2ggnetview(res)
  edges <- as.data.frame(igraph::as_data_frame(tg, what = "edges"))
  expect_true(nrow(edges) > 0)
  expect_true(all(edges$from %in% c("a", "b")))
  expect_true(all(edges$to %in% c("x", "y")))
})
