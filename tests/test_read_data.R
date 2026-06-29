# Test script for read_data()
# Regression guard: a .xls/.xlsx whose content is actually delimited text
# (no ZIP/OLE2 magic bytes) must be read as text, not handed to readxl.

source("R/read_write_data.R")

library(readr)

out_dir <- file.path(tempdir(), "read_data_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df_in <- data.frame(
  `#ID` = c("geneA", "geneB", "geneC"),
  T1R_1 = c(0.46, 4.42, 1.80),
  T1R_2 = c(0.70, 4.27, 3.08),
  check.names = FALSE
)

# --- Test 1: tab-delimited text disguised as .xls (the reported bug) ----------
cat("Test 1: tab-delimited text with .xls extension\n")
f_tsv <- file.path(out_dir, "All_gene_fpkm.xls")
write_tsv(df_in, f_tsv)
res1 <- read_data(f_tsv)
stopifnot(nrow(res1) == 3)
stopifnot(ncol(res1) == 3)
stopifnot(identical(attr(res1, "raw_names"), c("#ID", "T1R_1", "T1R_2")))
cat("✓ Test 1 passed — read", nrow(res1), "rows x", ncol(res1), "cols\n\n")

# --- Test 2: comma-delimited text disguised as .xls ---------------------------
cat("Test 2: comma-delimited text with .xls extension\n")
f_csv <- file.path(out_dir, "table.xls")
write_csv(df_in, f_csv)
res2 <- read_data(f_csv)
stopifnot(nrow(res2) == 3)
stopifnot(ncol(res2) == 3)
cat("✓ Test 2 passed — comma delimiter detected\n\n")

# --- Test 3: magic-byte detector ---------------------------------------------
cat("Test 3: is_excel_file magic-byte detection\n")
stopifnot(is_excel_file(f_tsv) == FALSE)   # text -> not Excel
stopifnot(is_excel_file(f_csv) == FALSE)
if (requireNamespace("writexl", quietly = TRUE)) {
  f_xlsx <- file.path(out_dir, "real.xlsx")
  writexl::write_xlsx(df_in, f_xlsx)
  stopifnot(is_excel_file(f_xlsx) == TRUE)  # ZIP magic -> Excel
  res3 <- read_data(f_xlsx)                 # genuine path still works
  stopifnot(nrow(res3) == 3)
  cat("✓ Test 3 passed — genuine .xlsx detected and read\n\n")
} else {
  cat("✓ Test 3 skipped — writexl not installed (genuine-xlsx path)\n\n")
}

# --- Test 4: error on genuinely unsupported format ---------------------------
cat("Test 4: unsupported extension errors\n")
file.create(file.path(out_dir, "nope.xyz"))  # must exist to reach ext switch
tryCatch({
  read_data(file.path(out_dir, "nope.xyz"))
  cat("✗ Test 4 failed: should have thrown an error\n")
}, error = function(e) {
  stopifnot(grepl("Unsupported format", conditionMessage(e)))
  cat("✓ Test 4 passed: correctly caught error -", conditionMessage(e), "\n\n")
})

cat("=====================================\n")
cat("All read_data tests passed! ✓\n")
cat("=====================================\n")
