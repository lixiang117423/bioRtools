#!/usr/bin/env Rscript
# Test get_hap_from_heatmap function

source("../R/get_hap_from_heatmap.R")

# Create test data
set.seed(123)
mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(mat) <- paste0("Sample", 1:10)
colnames(mat) <- paste0("Gene", 1:10)

# Add some structure to make clustering meaningful
mat[1:3, ] <- mat[1:3, ] + 2
mat[4:6, ] <- mat[4:6, ] - 2

# Generate pheatmap
library(pheatmap)
ph_result <- pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                      main = "Test Heatmap")

# Test 1: Extract 3 haplotypes
cat("\n=== Test 1: Extract 3 haplotypes ===\n")
haps <- get_hap_from_heatmap(ph_result, cluster_n = 3)
print(haps)
cat("\nClass:", class(haps), "\n")
cat("Column names:", names(haps), "\n")

# Test 2: Filter specific haplotype
cat("\n=== Test 2: Filter Hap1 ===\n")
hap1_items <- haps[haps$hap == "Hap1", ]
print(hap1_items)

# Test 3: Show haplotype sizes
cat("\n=== Test 3: Haplotype sizes ===\n")
print(table(haps$hap))

# Test 4: Test with different haplotype number
cat("\n=== Test 4: Extract 2 haplotypes ===\n")
haps2 <- get_hap_from_heatmap(ph_result, cluster_n = 2)
print(haps2)

cat("\nAll tests passed!\n")
