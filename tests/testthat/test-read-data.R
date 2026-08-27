# Migrated from tests/test_read_data.R (bare script).
# Regression guard: a .xls/.xlsx whose content is actually delimited text
# (no ZIP/OLE2 magic bytes) must be read as text, not handed to readxl.

is_excel_file <- bioRtools:::is_excel_file

out_dir <- file.path(tempdir(), "read_data_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df_in <- data.frame(
  `#ID` = c("geneA", "geneB", "geneC"),
  T1R_1 = c(0.46, 4.42, 1.80),
  T1R_2 = c(0.70, 4.27, 3.08),
  check.names = FALSE
)

gff_lines <- c(
  "##gff-version 3",
  "##sequence-region Chr01 1 10000",
  "Chr01\ttransdecoder\tgene\t4124\t6282\t.\t+\t.\tID=gene1;Name=gene1",
  "###",
  "Chr01\ttransdecoder\tCDS\t4271\t4358\t5.0\t+\t0\tID=cds1;Note=foo#bar"
)

test_that("tab-delimited text disguised as .xls is read as text", {
  f_tsv <- file.path(out_dir, "All_gene_fpkm.xls")
  readr::write_tsv(df_in, f_tsv)
  res1 <- read_data(f_tsv)
  expect_equal(nrow(res1), 3)
  expect_equal(ncol(res1), 3)
  expect_identical(attr(res1, "raw_names"), c("#ID", "T1R_1", "T1R_2"))
})

test_that("comma-delimited text disguised as .xls is read as text", {
  f_csv <- file.path(out_dir, "table.xls")
  readr::write_csv(df_in, f_csv)
  res2 <- read_data(f_csv)
  expect_equal(nrow(res2), 3)
  expect_equal(ncol(res2), 3)
})

test_that("is_excel_file magic-byte detection", {
  f_tsv <- file.path(out_dir, "All_gene_fpkm.xls")
  f_csv <- file.path(out_dir, "table.xls")
  expect_false(is_excel_file(f_tsv))   # text -> not Excel
  expect_false(is_excel_file(f_csv))
})

test_that("genuine .xlsx is detected and read", {
  skip_if_not_installed("writexl")
  f_xlsx <- file.path(out_dir, "real.xlsx")
  writexl::write_xlsx(df_in, f_xlsx)
  expect_true(is_excel_file(f_xlsx))   # ZIP magic -> Excel
  res3 <- read_data(f_xlsx)            # genuine path still works
  expect_equal(nrow(res3), 3)
})

test_that("unsupported extension errors", {
  file.create(file.path(out_dir, "nope.xyz"))  # must exist to reach ext switch
  expect_error(read_data(file.path(out_dir, "nope.xyz")), "Unsupported format")
})

test_that("non-integer skip is rejected; skip = 1 still reads cleanly", {
  f_skip <- file.path(out_dir, "skip.tsv")
  readr::write_tsv(df_in, f_skip)
  expect_error(read_data(f_skip, skip = "--"), "single non-negative integer")
  res5 <- read_data(f_skip, skip = 1)
  expect_equal(nrow(res5), 2)
})

test_that("GFF3 (.gff3) read into a 9-column data frame", {
  f_gff <- file.path(out_dir, "test.gff3")
  writeLines(gff_lines, f_gff)
  res6 <- read_data(f_gff)
  expect_equal(nrow(res6), 2)
  expect_equal(ncol(res6), 9)
  expect_identical(
    names(res6),
    c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  )
  expect_true(is.integer(res6$start))
  expect_true(is.integer(res6$end))
  expect_true(is.na(res6$score[1]))        # "." -> NA
  expect_equal(res6$score[2], 5.0)         # 5.0 kept
  expect_true(is.na(res6$phase[1]))        # "." -> NA
  expect_equal(res6$phase[2], 0L)          # CDS phase kept
  expect_true(is.na(res6$strand[1]) || res6$strand[1] == "+")  # "." strand -> NA
  expect_equal(res6$attributes[2], "ID=cds1;Note=foo#bar")     # "#" preserved, not truncated
})

test_that(".gff shorthand reads identically", {
  f_gff2 <- file.path(out_dir, "test.gff")
  writeLines(gff_lines, f_gff2)
  res6b <- read_data(f_gff2)
  expect_identical(dim(res6b), c(2L, 9L))
})
