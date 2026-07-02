# Tests for cor2gephi(): cor_analysis edge list -> Gephi nodes/edges CSV.
source("R/cor_result_to_graph.R")
source("R/cor2gephi.R")
source("R/cor_analysis.R")

suppressPackageStartupMessages({
  library(igraph); library(dplyr)
})

# Symmetric within-dataset correlations -> reciprocal edges (A-B and B-A).
X <- data.frame(
  a = c(1, 2, 3, 4, 5, 6),
  b = c(2, 4, 6, 8, 10, 12),   # a,b perfectly correlated
  c = c(6, 5, 4, 3, 2, 1),    # a,c perfectly anti-correlated
  d = c(1, 1, 1, 1, 1, 1)     # zero variance -> dropped by cor_analysis
)

cor_res <- cor_analysis(X, method = "pearson", cor = 0.8, pvalue = 0.05)
tmp <- tempfile()
on.exit(unlink(paste0(tmp, c("-nodes.csv", "-edges.csv"))), add = TRUE)

# --- Test B1a: validation rejects bad input -------------------------------
cat("Test B1a: cor2gephi validates its input\n")
ok <- try(cor2gephi(data.frame()), silent = TRUE)
stopifnot(inherits(ok, "try-error"))
cat("✓ Test B1a passed\n\n")

# --- Test B1b: reciprocal edges are deduplicated --------------------------
cat("Test B1b: dedup removes reciprocal edges\n")
out <- cor2gephi(cor_res, prefix = tmp, enrich = FALSE)
edges <- out$edges
# Each unordered pair appears exactly once.
key <- paste(pmin(edges$Source, edges$Target), pmax(edges$Source, edges$Target))
stopifnot(!any(duplicated(key)))
stopifnot(all(c("Source", "Target", "Weight", "Type", "Correlation",
                "Direction") %in% names(edges)))
stopifnot(all(edges$Type == "undirected"))
cat("✓ Test B1b passed\n\n")

# --- Test B1c: files are written -----------------------------------------
cat("Test B1c: nodes/edges CSV files written\n")
stopifnot(file.exists(paste0(tmp, "-nodes.csv")))
stopifnot(file.exists(paste0(tmp, "-edges.csv")))
nodes <- out$nodes
stopifnot(all(c("Id", "Label", "Degree", "Strength") %in% names(nodes)))
stopifnot(!"Eigenvector" %in% names(nodes))   # enrich = FALSE
cat("✓ Test B1c passed\n\n")

# --- Test B1d: weight = "signed" puts signed cor in Weight ---------------
cat("Test B1d: weight = signed\n")
out_s <- cor2gephi(cor_res, prefix = tmp, weight = "signed", enrich = FALSE)
stopifnot(all(out_s$edges$Weight == out_s$edges$Correlation))
cat("✓ Test B1d passed\n\n")

# --- Test B2: enrich = TRUE adds topology columns -------------------------
cat("Test B2: enrich columns + community options\n")
out_e <- cor2gephi(cor_res, prefix = tmp, enrich = TRUE, community = "auto")
stopifnot(all(c("Eigenvector", "Modularity", "IsHub") %in% names(out_e$nodes)))
stopifnot(sum(out_e$nodes$IsHub) == length(unique(out_e$nodes$Modularity)))

out_n <- cor2gephi(cor_res, prefix = tmp, enrich = TRUE, community = "none")
stopifnot(!"Modularity" %in% names(out_n$nodes))
stopifnot("Eigenvector" %in% names(out_n$nodes))
cat("✓ Test B2 passed\n\n")

cat("=====================================\n")
cat("All cor2gephi tests passed! ✓\n")
cat("=====================================\n")
