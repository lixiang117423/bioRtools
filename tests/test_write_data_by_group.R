# Test script for write_data_by_group function
# Writes one file per group, reads each back, and compares.

source("R/read_write_data.R")

library(dplyr)

# Build a small grouped data frame
df <- tibble::tibble(
  CHROM = c(rep("chr1", 3), rep("chr2", 2)),
  pos = c(100, 200, 300, 400, 500),
  p_value = c(0.1, 0.02, 0.5, 0.001, 0.9)
)

out_dir <- file.path(tempdir(), "write_data_by_group_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Data loaded. Rows:", nrow(df), " Columns:",
    paste(names(df), collapse = ", "), "\n\n")

# Test 1: single-column group, xlsx
cat("Test 1: single-column group, xlsx\n")
df %>%
  dplyr::group_by(CHROM) %>%
  write_data_by_group(file.path(out_dir, "result_{CHROM}.xlsx"))
f1 <- file.path(out_dir, "result_chr1.xlsx")
f2 <- file.path(out_dir, "result_chr2.xlsx")
stopifnot(file.exists(f1), file.exists(f2))
stopifnot(nrow(read_data(f1)) == 3)
stopifnot(nrow(read_data(f2)) == 2)
cat("✓ Test 1 passed\n\n")

# Test 2: multi-column group + csv passthrough
cat("Test 2: multi-column group, csv\n")
df$TYPE <- c("A", "A", "B", "B", "B")
df %>%
  dplyr::group_by(CHROM, TYPE) %>%
  write_data_by_group(file.path(out_dir, "{CHROM}_{TYPE}.csv"))
stopifnot(file.exists(file.path(out_dir, "chr1_A.csv")))
stopifnot(file.exists(file.path(out_dir, "chr2_B.csv")))
cat("✓ Test 2 passed\n\n")

# Test 3: returns input invisibly (pipeable)
cat("Test 3: returns input invisibly\n")
res <- df %>%
  dplyr::group_by(CHROM) %>%
  write_data_by_group(file.path(out_dir, "p_{CHROM}.xlsx"))
stopifnot(dplyr::is_grouped_df(res))
cat("✓ Test 3 passed\n\n")

# Test 4: error on ungrouped data
cat("Test 4: error on ungrouped data\n")
tryCatch({
  write_data_by_group(df, file.path(out_dir, "x_{CHROM}.xlsx"))
  cat("✗ Test 4 failed: should have thrown an error\n")
}, error = function(e) {
  cat("✓ Test 4 passed: correctly caught error -", conditionMessage(e), "\n\n")
})

cat("=====================================\n")
cat("All tests completed successfully! ✓\n")
cat("=====================================\n")
