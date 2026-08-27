# Migrated from tests/test_cor2gephi.R (bare script).
# cor2gephi(): cor_analysis edge list -> Gephi nodes/edges CSV.

make_cor_input <- function() {
  # Symmetric within-dataset correlations -> reciprocal edges (A-B and B-A).
  X <- data.frame(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(2, 4, 6, 8, 10, 12),   # a,b perfectly correlated
    c = c(6, 5, 4, 3, 2, 1),    # a,c perfectly anti-correlated
    d = c(1, 1, 1, 1, 1, 1)     # zero variance -> dropped by cor_analysis
  )
  cor_analysis(X, method = "pearson", cor = 0.8, pvalue = 0.05)
}

test_that("cor2gephi validates its input", {
  expect_error(cor2gephi(data.frame()))
})

test_that("dedup removes reciprocal edges", {
  tmp <- tempfile()
  on.exit(unlink(paste0(tmp, c("-nodes.csv", "-edges.csv"))), add = TRUE)

  out <- cor2gephi(make_cor_input(), prefix = tmp, enrich = FALSE)
  edges <- out$edges
  # Each unordered pair appears exactly once.
  key <- paste(pmin(edges$Source, edges$Target),
               pmax(edges$Source, edges$Target))
  expect_false(any(duplicated(key)))
  expect_true(all(c("Source", "Target", "Weight", "Type", "Correlation",
                    "Direction") %in% names(edges)))
  expect_true(all(edges$Type == "undirected"))
})

test_that("nodes/edges CSV files written", {
  tmp <- tempfile()
  on.exit(unlink(paste0(tmp, c("-nodes.csv", "-edges.csv"))), add = TRUE)

  out <- cor2gephi(make_cor_input(), prefix = tmp, enrich = FALSE)

  expect_true(file.exists(paste0(tmp, "-nodes.csv")))
  expect_true(file.exists(paste0(tmp, "-edges.csv")))
  nodes <- out$nodes
  expect_true(all(c("Id", "Label", "Degree", "Strength") %in% names(nodes)))
  expect_false("Eigenvector" %in% names(nodes))   # enrich = FALSE
})

test_that("weight = signed puts signed cor in Weight", {
  tmp <- tempfile()
  on.exit(unlink(paste0(tmp, c("-nodes.csv", "-edges.csv"))), add = TRUE)

  out_s <- cor2gephi(make_cor_input(), prefix = tmp, weight = "signed",
                     enrich = FALSE)
  expect_true(all(out_s$edges$Weight == out_s$edges$Correlation))
})

test_that("enrich columns + community options", {
  tmp <- tempfile()
  on.exit(unlink(paste0(tmp, c("-nodes.csv", "-edges.csv"))), add = TRUE)

  out_e <- cor2gephi(make_cor_input(), prefix = tmp, enrich = TRUE,
                     community = "auto")
  expect_true(all(c("Eigenvector", "Modularity", "IsHub") %in% names(out_e$nodes)))
  expect_equal(sum(out_e$nodes$IsHub),
               length(unique(out_e$nodes$Modularity)))

  out_n <- cor2gephi(make_cor_input(), prefix = tmp, enrich = TRUE,
                     community = "none")
  expect_false("Modularity" %in% names(out_n$nodes))
  expect_true("Eigenvector" %in% names(out_n$nodes))
})
