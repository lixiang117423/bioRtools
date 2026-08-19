test_that("read_many row-binds files with an id column", {
  dir <- file.path(tempdir(), "bioRtools_read_many")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f1 <- file.path(dir, "S1.txt")
  f2 <- file.path(dir, "S2.txt")
  writeLines(c("gene	count", "a	1", "b	2"), f1)
  writeLines(c("gene	count", "a	10", "b	20"), f2)

  out <- suppressMessages(read_many(c(f1, f2), id_col = "sample"))
  expect_equal(names(out), c("sample", "gene", "count"))
  expect_equal(out$sample, c("S1.txt", "S1.txt", "S2.txt", "S2.txt"))
  expect_equal(out$count, c(1, 2, 10, 20))

  out_path <- suppressMessages(
    read_many(c(f1, f2), id_col = "sample", id_from = "path")
  )
  expect_equal(out_path$sample, c(f1, f1, f2, f2))
})

test_that("read_many forwards extra arguments to read_delim", {
  dir <- file.path(tempdir(), "bioRtools_read_many_cols")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f1 <- file.path(dir, "a.txt")
  f2 <- file.path(dir, "b.txt")
  writeLines(c("g1	1", "g2	2"), f1)
  writeLines(c("g1	10", "g2	20"), f2)

  out <- suppressMessages(
    read_many(c(f1, f2), id_col = "file", col_names = c("gene", "count"))
  )
  expect_equal(names(out), c("file", "gene", "count"))
  expect_equal(out$count, c(1, 2, 10, 20))
})

test_that("read_many validates inputs", {
  expect_error(read_many(character(0)), "non-empty character vector")
  expect_error(read_many("no_such_file.txt"), "not found")
  expect_error(read_many("a.txt", id_col = ""), "single non-empty string")
  expect_error(read_many("a.txt", id_from = "nope"), "should be one of")
  expect_error(read_many("a.txt", delim = ""), "single non-empty string")
})
