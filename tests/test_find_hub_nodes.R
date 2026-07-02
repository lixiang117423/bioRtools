# Tests for find_hub_nodes(): cor_analysis edge list + microbiome_net result.
source("R/cor_result_to_graph.R")
source("R/find_hub_nodes.R")

suppressPackageStartupMessages({ library(igraph) })

# Known 5-node edge list. Absolute strengths:
# A = .9+.8+.7+.3 = 2.7, B = .9, C = .8, D = .7, E = .3
edges <- data.frame(
  from = c("A", "A", "A", "A"),
  to   = c("B", "C", "D", "E"),
  cor  = c(0.9, 0.8, 0.7, 0.3),
  pvalue = c(0.001, 0.001, 0.001, 0.05),
  stringsAsFactors = FALSE
)

# --- Test 1: edge-list path centrality + top_percent ----------------------
cat("Test 1: edge-list path centrality + top_percent\n")
res <- find_hub_nodes(edges, method = "top_percent", top_percent = 0.1, weight = "absolute")
stopifnot(all(c("node", "degree", "strength", "betweenness", "closeness",
                "eigenvector", "is_hub") %in% names(res)))
stopifnot(nrow(res) == 5)
stopifnot(res$node[1] == "A")
stopifnot(round(res$strength[res$node == "A"], 6) == 2.7)
stopifnot(sum(res$is_hub) == 1)                # ceiling(5 * 0.1) = 1
stopifnot(res$is_hub[res$node == "A"] == TRUE)
cat("✓ Test 1 passed\n\n")

# --- Test 2: top_n mode ---------------------------------------------------
cat("Test 2: top_n mode\n")
res2 <- find_hub_nodes(edges, method = "top_n", top_n = 2)
stopifnot(sum(res2$is_hub) == 2)
stopifnot(setequal(res2$node[res2$is_hub], c("A", "B")))
cat("✓ Test 2 passed\n\n")

# --- Test 3: threshold mode ----------------------------------------------
cat("Test 3: threshold mode\n")
res3 <- find_hub_nodes(edges, method = "threshold", strength_threshold = 0.75)
stopifnot(setequal(res3$node[res3$is_hub], c("A", "B", "C")))   # 2.7, 0.9, 0.8 >= 0.75
cat("✓ Test 3 passed\n\n")

# --- Test 4: tie at cutoff is included -----------------------------------
cat("Test 4: tie at cutoff included\n")
edges_tie <- data.frame(
  from = c("A", "A", "B", "B"), to = c("B", "C", "C", "D"),
  cor = c(0.5, 0.5, 0.5, 0.5), pvalue = rep(0.001, 4),
  stringsAsFactors = FALSE
)
# strengths: A=1.0, B=1.5, C=1.0, D=0.5  -> ranks: B=1, A&C=2 (tie), D=4
res4 <- find_hub_nodes(edges_tie, method = "top_n", top_n = 2)
stopifnot(setequal(res4$node[res4$is_hub], c("B", "A", "C")))   # tie at rank 2 included
cat("✓ Test 4 passed\n\n")

# --- Test 5: microbiome_net path (synthetic, no SpiecEasi) ----------------
cat("Test 5: microbiome_net result path\n")
g_syn <- igraph::make_graph(c(1, 2, 1, 3, 1, 4, 2, 3), directed = FALSE)
igraph::V(g_syn)$name <- c("A", "B", "C", "D")
igraph::E(g_syn)$weight <- c(0.9, 0.8, 0.7, 0.6)   # absolute weights
igraph::E(g_syn)$sign <- c(1, 1, -1, 1)            # correlation signs
net_syn <- list(networks = list(groupA = g_syn))
# absolute strengths: A=2.4, B=1.5, C=1.4, D=0.7
res5a <- find_hub_nodes(net_syn, group = "groupA", method = "top_percent",
                        top_percent = 0.5, weight = "absolute")
stopifnot(nrow(res5a) == 4)
stopifnot(res5a$node[1] == "A")
# signed strengths: B=1.5, C=1.4, A=1.0, D=-0.7  -> differs from absolute
res5b <- find_hub_nodes(net_syn, group = "groupA", method = "top_percent",
                        top_percent = 0.5, weight = "signed")
stopifnot(nrow(res5b) == 4)
stopifnot(!identical(res5a$strength, res5b$strength))
stopifnot(res5b$node[1] == "B")
cat("✓ Test 5 passed\n\n")

# --- Test 6: validation ---------------------------------------------------
cat("Test 6: input validation\n")
stopifnot(inherits(try(find_hub_nodes(123), silent = TRUE), "try-error"))
stopifnot(inherits(try(find_hub_nodes(data.frame()), silent = TRUE), "try-error"))
stopifnot(inherits(try(find_hub_nodes(net_syn, group = NULL), silent = TRUE), "try-error"))
stopifnot(inherits(try(find_hub_nodes(net_syn, group = "missing"), silent = TRUE), "try-error"))
stopifnot(inherits(try(find_hub_nodes(edges, method = "threshold"), silent = TRUE), "try-error"))
cat("✓ Test 6 passed\n\n")

cat("=====================================\n")
cat("All find_hub_nodes tests passed! ✓\n")
cat("=====================================\n")
