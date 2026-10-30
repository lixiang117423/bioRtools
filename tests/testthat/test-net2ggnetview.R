if (!can_load_namespace("tidygraph")) {
  test_that("optional dependency integration tests are skipped", {
    testthat::skip("tidygraph is not available")
  })
} else {
# net2ggnetview(): microbiome_net result -> ggNetView-compatible tbl_graph.

# Synthetic microbiome_net result, following the fixture pattern of
# test-find-hub-nodes.R (no SpiecEasi needed).
make_net_result <- function() {
  g <- igraph::make_graph(c("A", "B", "A", "C", "A", "D", "B", "C"),
                          directed = FALSE)
  igraph::E(g)$weight <- c(0.9, 0.8, 0.7, 0.6)   # absolute weights
  igraph::E(g)$sign <- c(1, 1, -1, 1)            # correlation signs
  list(
    networks = list(groupA = g),
    adja_mats = list(groupA = NULL),
    node_props = data.frame(
      feature = c("A", "B", "C", "D"),
      group = "groupA",
      cluster = c(1, 1, 2, 2),
      stringsAsFactors = FALSE
    )
  )
}

test_that("returns a tbl_graph with node/edge schema", {
  tg <- net2ggnetview(make_net_result(), group = "groupA")
  expect_s3_class(tg, "tbl_graph")

  nodes <- as.data.frame(igraph::as_data_frame(tg, what = "vertices"))
  edges <- as.data.frame(igraph::as_data_frame(tg, what = "edges"))

  expect_true(all(c("name", "degree", "strength", "modularity_factor",
                    "modularity", "modularity2", "modularity3") %in% names(nodes)))
  expect_true(all(c("from", "to", "correlation", "weight",
                    "corr_direction") %in% names(edges)))
  # correlation is signed; weight is its absolute value
  expect_true(all(edges$weight == abs(edges$correlation)))
  expect_setequal(sign(edges$correlation[edges$from == "A" & edges$to == "D"]), -1)
})

test_that("nodes carry capitalized Degree/Strength/Modularity for ggNetView layouts", {
  # ggNetView(layout = "gephi") internally calls arrange(modularity3, desc(Degree))
  # and maps size = Degree, group.by = "Modularity" — lowercase-only columns
  # make ggNetView fail with "object 'Degree' not found".
  tg <- net2ggnetview(make_net_result(), group = "groupA")
  nodes <- as.data.frame(igraph::as_data_frame(tg, what = "vertices"))

  expect_true(all(c("Modularity", "Degree", "Strength") %in% names(nodes)))
  expect_equal(nodes$Degree, nodes$degree)
  expect_equal(nodes$Strength, nodes$strength)
})
}
