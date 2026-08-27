# Migrated from tests/test_find_hub_nodes.R (bare script).
# find_hub_nodes(): cor_analysis edge list + microbiome_net result.

# Known 5-node edge list. Absolute strengths:
# A = .9+.8+.7+.3 = 2.7, B = .9, C = .8, D = .7, E = .3
edges <- data.frame(
  from = c("A", "A", "A", "A"),
  to   = c("B", "C", "D", "E"),
  cor  = c(0.9, 0.8, 0.7, 0.3),
  pvalue = c(0.001, 0.001, 0.001, 0.05),
  stringsAsFactors = FALSE
)

test_that("edge-list path centrality + top_percent", {
  res <- find_hub_nodes(edges, method = "top_percent", top_percent = 0.1,
                        weight = "absolute")

  expect_true(all(c("node", "degree", "strength", "betweenness", "closeness",
                    "eigenvector", "is_hub") %in% names(res)))
  expect_equal(nrow(res), 5)
  expect_equal(res$node[1], "A")
  expect_true(round(res$strength[res$node == "A"], 6) == 2.7)
  expect_equal(sum(res$is_hub), 1)   # ceiling(5 * 0.1) = 1
  expect_true(res$is_hub[res$node == "A"])
})

test_that("top_n mode", {
  res2 <- find_hub_nodes(edges, method = "top_n", top_n = 2)

  expect_equal(sum(res2$is_hub), 2)
  expect_setequal(res2$node[res2$is_hub], c("A", "B"))
})

test_that("threshold mode", {
  res3 <- find_hub_nodes(edges, method = "threshold", strength_threshold = 0.75)

  expect_setequal(res3$node[res3$is_hub], c("A", "B", "C"))   # 2.7, 0.9, 0.8 >= 0.75
})

test_that("tie at cutoff is included", {
  edges_tie <- data.frame(
    from = c("A", "A", "B", "B"), to = c("B", "C", "C", "D"),
    cor = c(0.5, 0.5, 0.5, 0.5), pvalue = rep(0.001, 4),
    stringsAsFactors = FALSE
  )
  # strengths: A=1.0, B=1.5, C=1.0, D=0.5  -> ranks: B=1, A&C=2 (tie), D=4
  res4 <- find_hub_nodes(edges_tie, method = "top_n", top_n = 2)

  expect_setequal(res4$node[res4$is_hub], c("B", "A", "C"))   # tie at rank 2 included
})

test_that("microbiome_net result path", {
  # Synthetic network, no SpiecEasi needed.
  g_syn <- igraph::make_graph(c(1, 2, 1, 3, 1, 4, 2, 3), directed = FALSE)
  igraph::V(g_syn)$name <- c("A", "B", "C", "D")
  igraph::E(g_syn)$weight <- c(0.9, 0.8, 0.7, 0.6)   # absolute weights
  igraph::E(g_syn)$sign <- c(1, 1, -1, 1)            # correlation signs
  net_syn <- list(networks = list(groupA = g_syn))

  # absolute strengths: A=2.4, B=1.5, C=1.4, D=0.7
  res5a <- find_hub_nodes(net_syn, group = "groupA", method = "top_percent",
                          top_percent = 0.5, weight = "absolute")
  expect_equal(nrow(res5a), 4)
  expect_equal(res5a$node[1], "A")

  # signed strengths: B=1.5, C=1.4, A=1.0, D=-0.7  -> differs from absolute
  res5b <- find_hub_nodes(net_syn, group = "groupA", method = "top_percent",
                          top_percent = 0.5, weight = "signed")
  expect_equal(nrow(res5b), 4)
  expect_false(identical(res5a$strength, res5b$strength))
  expect_equal(res5b$node[1], "B")
})

test_that("input validation", {
  g_syn <- igraph::make_graph(c(1, 2, 1, 3, 1, 4, 2, 3), directed = FALSE)
  igraph::V(g_syn)$name <- c("A", "B", "C", "D")
  igraph::E(g_syn)$weight <- c(0.9, 0.8, 0.7, 0.6)
  net_syn <- list(networks = list(groupA = g_syn))

  expect_error(find_hub_nodes(123))
  expect_error(find_hub_nodes(data.frame()))
  expect_error(find_hub_nodes(net_syn, group = NULL))
  expect_error(find_hub_nodes(net_syn, group = "missing"))
  expect_error(find_hub_nodes(edges, method = "threshold"))

  # edge-less graph errors cleanly
  g_empty <- igraph::make_empty_graph(n = 2, directed = FALSE)
  igraph::V(g_empty)$name <- c("A", "B")
  net_empty <- list(networks = list(g = g_empty))
  expect_error(find_hub_nodes(net_empty, group = "g"))

  # bad criterion parameters
  expect_error(find_hub_nodes(edges, method = "top_percent", top_percent = 1.5))
  expect_error(find_hub_nodes(edges, method = "top_n", top_n = 0))
})
