test_that("split_pick picks the requested segment", {
  paths <- c("rawdata/S1/clean.fq", "rawdata/S2/clean.fq")
  expect_equal(split_pick(paths, index = 2), c("S1", "S2"))
  expect_equal(split_pick(paths), c("clean.fq", "clean.fq"))
})

test_that("split_pick works with any separator and position", {
  expect_equal(split_pick(c("a_b_c", "d_e_f"), sep = "_", index = 3),
               c("c", "f"))
  expect_equal(split_pick(c("chr1:100-200", "chr2:300-400"),
                          sep = ":", index = 1),
               c("chr1", "chr2"))
})

test_that("split_pick returns NA for too-short strings", {
  expect_equal(split_pick("a/b", index = 5), NA_character_)
  expect_equal(split_pick(c("a/b", "a/b/c"), index = 3), c(NA, "c"))
})

test_that("split_pick separator is literal, not a regular expression", {
  # A regex "." would match every character and return "" here
  expect_equal(split_pick("a.b.c", sep = ".", index = 1), "a")
})

test_that("split_pick optional second split keeps a sub-segment", {
  paths <- c("rawdata/S1/clean.fq", "rawdata/S2/clean.fq")
  expect_equal(split_pick(paths, split2 = ".", index2 = 1),
               c("clean", "clean"))
  expect_equal(split_pick(paths, split2 = ".", index2 = 2), c("fq", "fq"))
})

test_that("split_pick validates inputs", {
  expect_error(split_pick(1:3), "'x' must be a character vector")
  expect_error(split_pick("a/b", sep = c("/", "\\")), "single non-empty string")
  expect_error(split_pick("a/b", index = 0), "positive integer")
  expect_error(split_pick("a/b", index = 1.5), "positive integer")
  expect_error(split_pick("a/b", split2 = ""), "NULL or a single non-empty string")
  expect_error(split_pick("a/b", index2 = 0), "positive integer")
})

test_that("split_pick handles NA and empty input", {
  expect_equal(split_pick(c("a/b/c", NA)), c("c", NA))
  expect_equal(split_pick(character(0)), character(0))
})
