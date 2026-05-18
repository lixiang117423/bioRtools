#!/usr/bin/env Rscript
# Test extract_pheatmap_clusters function

source("../R/extract_pheatmap_clusters.R")

# Create test data
set.seed(123)
mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(mat) <- paste0("Gene", 1:10)
colnames(mat) <- paste0("Sample", 1:10)

# Add some structure to make clustering meaningful
mat[1:3, ] <- mat[1:3, ] + 2
mat[4:6, ] <- mat[4:6, ] - 2

# Generate pheatmap
library(pheatmap)
ph_result <- pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                      main = "Test Heatmap")

# Test 1: Extract 3 clusters
cat("\n=== Test 1: Extract 3 clusters ===\n")
clusters <- extract_pheatmap_clusters(ph_result, cluster_n = 3)
print(clusters)
cat("\nClass:", class(clusters), "\n")

# Test 2: Filter specific cluster
cat("\n=== Test 2: Filter Cluster1 ===\n")
cluster1_items <- clusters[clusters$cluster == "Cluster1", ]
print(cluster1_items)

# Test 3: Show cluster sizes
cat("\n=== Test 3: Cluster sizes ===\n")
print(table(clusters$cluster))

# Test 4: Test with different cluster number
cat("\n=== Test 4: Extract 2 clusters ===\n")
clusters2 <- extract_pheatmap_clusters(ph_result, cluster_n = 2)
print(clusters2)

cat("\nAll tests passed!\n")
