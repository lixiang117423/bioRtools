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

# --- Test 5: non-integer skip rejected with a clear error --------------------
cat("Test 5: non-integer skip raises an actionable error\n")
f_skip <- file.path(out_dir, "skip.tsv")
write_tsv(df_in, f_skip)
tryCatch({
  read_data(f_skip, skip = "--")
  cat("✗ Test 5 failed: should have thrown on skip = \"--\"\n")
}, error = function(e) {
  msg <- conditionMessage(e)
  stopifnot(grepl("single non-negative integer", msg))
  cat("✓ Test 5 passed: skip guard caught -", msg, "\n\n")
})
res5 <- read_data(f_skip, skip = 1)
stopifnot(nrow(res5) == 2)
cat("✓ Test 5b passed: skip = 1 still reads cleanly\n\n")

# --- Test 6: GFF3 parse ------------------------------------------------------
# Regression guard: readr's comment = "#" truncates any "#" inside an attribute
# value, so comment lines must be filtered by line-start and data lines left
# untouched; "." must become NA; types must be integer/double as appropriate.
cat("Test 6: GFF3 (.gff3) read into a 9-column data frame\n")
gff_lines <- c(
  "##gff-version 3",
  "##sequence-region Chr01 1 10000",
  "Chr01\ttransdecoder\tgene\t4124\t6282\t.\t+\t.\tID=gene1;Name=gene1",
  "###",
  "Chr01\ttransdecoder\tCDS\t4271\t4358\t5.0\t+\t0\tID=cds1;Note=foo#bar"
)
f_gff <- file.path(out_dir, "test.gff3")
writeLines(gff_lines, f_gff)
res6 <- read_data(f_gff)
stopifnot(nrow(res6) == 2)
stopifnot(ncol(res6) == 9)
stopifnot(identical(
  names(res6),
  c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
))
stopifnot(is.integer(res6$start), is.integer(res6$end))
stopifnot(is.na(res6$score[1]), res6$score[2] == 5.0)   # "." -> NA, 5.0 kept
stopifnot(is.na(res6$phase[1]), res6$phase[2] == 0L)    # "." -> NA, CDS phase kept
stopifnot(is.na(res6$strand[1]) || res6$strand[1] == "+")  # "." strand -> NA
stopifnot(res6$attributes[2] == "ID=cds1;Note=foo#bar") # "#" preserved, not truncated
cat("✓ Test 6 passed — read", nrow(res6), "rows; # preserved; . -> NA\n\n")

# .gff shorthand maps to the same reader
f_gff2 <- file.path(out_dir, "test.gff")
file.copy(f_gff, f_gff2, overwrite = TRUE)
res6b <- read_data(f_gff2)
stopifnot(identical(dim(res6b), c(2L, 9L)))
cat("✓ Test 6b passed: .gff shorthand reads identically\n\n")

cat("=====================================\n")
cat("All read_data tests passed! ✓\n")
cat("=====================================\n")
