make_admix_dir <- function() {
  dir <- tempfile("admix")
  dir.create(dir)
  writeLines(c(
    "F1\tS1\t0\t0\t0\t1",
    "F2\tS2\t0\t0\t0\t1",
    "F3\tS3\t0\t0\t0\t1",
    "F4\tS4\t0\t0\t0\t1",
    "F5\tS5\t0\t0\t0\t1"
  ), file.path(dir, "test.fam"))
  writeLines(c("0.9\t0.1", "0.8\t0.2", "0.7\t0.3", "0.1\t0.9", "0.2\t0.8"),
             file.path(dir, "test.2.Q"))
  writeLines(c("0.8\t0.1\t0.1", "0.7\t0.2\t0.1", "0.6\t0.3\t0.1",
               "0.1\t0.2\t0.7", "0.2\t0.1\t0.7"),
             file.path(dir, "test.3.Q"))
  writeLines("(('S1':1,'S2':1):1,('S3':1,'S4':1):1,'S5':2);",
             file.path(dir, "tree.nwk"))
  dir
}

test_that("mode 1 (no tree) uses IID column and returns admixture plot only", {
  dir <- make_admix_dir()
  res <- admixture_phylo_analysis(dir)

  expect_s3_class(res, "admixture_phylo_results")
  expect_equal(res$samples, c("S1", "S2", "S3", "S4", "S5"))
  expect_equal(res$k_values, c(2, 3))
  expect_s3_class(res$admixture_plot, "ggplot")
  expect_null(res$tree_plot)
  expect_null(res$combined_plot)
})

test_that("fam_id_col = 1 takes the FID column", {
  dir <- make_admix_dir()
  res <- admixture_phylo_analysis(dir, fam_id_col = 1)
  expect_equal(res$samples, c("F1", "F2", "F3", "F4", "F5"))
})

test_that("mode 2 accepts a phylo object and combines plots", {
  dir <- make_admix_dir()
  tr <- ape::read.tree(file.path(dir, "tree.nwk"))

  res <- suppressWarnings(admixture_phylo_analysis(dir, tree_file = tr))

  expect_equal(res$matched_samples, c("S1", "S2", "S3", "S4", "S5"))
  expect_s3_class(res$tree_plot, "ggtree")
  expect_false(is.null(res$combined_plot))
})

test_that("mode 2 accepts a tree file path and k_range filter", {
  dir <- make_admix_dir()
  res <- suppressWarnings(admixture_phylo_analysis(
    dir, tree_file = file.path(dir, "tree.nwk"), k_range = 2
  ))
  expect_equal(res$k_values, 2)
  expect_false(is.null(res$combined_plot))
})

test_that("sample_group colors tree tips (data frame input)", {
  dir <- make_admix_dir()
  sg <- data.frame(label = c("S1", "S2", "S3", "S4", "S5"),
                   group = c("A", "A", "B", "B", "C"))
  res <- suppressWarnings(admixture_phylo_analysis(
    dir, tree_file = file.path(dir, "tree.nwk"), sample_group = sg
  ))
  expect_s3_class(res$tree_plot, "ggtree")
})

test_that("save writes the main plot to disk", {
  dir <- make_admix_dir()
  out <- file.path(dir, "out", "combined.png")
  suppressWarnings(admixture_phylo_analysis(dir, save = out))
  expect_true(file.exists(out))
})

test_that("partial sample overlap raises warnings, not silence", {
  dir <- make_admix_dir()
  tr <- ape::read.tree(text = "((S1:1,S2:1):1,(S3:1,S4:1):1,X5:2);")
  expect_warning(
    admixture_phylo_analysis(dir, tree_file = tr),
    "in the tree but not in .Q"
  )
})

test_that("zero sample overlap stops with both sample lists", {
  dir <- make_admix_dir()
  tr <- ape::read.tree(text = "((X1:1,X2:1):1,(X3:1,X4:1):1,X5:2);")
  expect_error(admixture_phylo_analysis(dir, tree_file = tr),
               "No samples match")
  expect_error(admixture_phylo_analysis(dir, tree_file = tr),
               "Tree samples")
})

test_that("Q/fam row mismatch stops with the offending file name", {
  dir <- tempfile("admix_bad")
  dir.create(dir)
  writeLines(c("F1\tS1\t0\t0\t0\t1", "F2\tS2\t0\t0\t0\t1"), file.path(dir, "t.fam"))
  writeLines(c("0.9\t0.1", "0.8\t0.2", "0.7\t0.3"), file.path(dir, "t.2.Q"))
  expect_error(admixture_phylo_analysis(dir), "Row mismatch in t.2.Q")
})

test_that("Q file names without K digits stop with a clear message", {
  dir <- tempfile("admix_nok")
  dir.create(dir)
  writeLines("F1\tS1\t0\t0\t0\t1", file.path(dir, "t.fam"))
  writeLines("0.9\t0.1", file.path(dir, "t.Q"))
  expect_error(admixture_phylo_analysis(dir), "Cannot extract K")
})

test_that("k_range with unavailable values warns and keeps the rest", {
  dir <- make_admix_dir()
  expect_warning(
    res <- admixture_phylo_analysis(dir, k_range = c(2, 9)),
    "not found in .Q files"
  )
  expect_equal(res$k_values, 2)
})

test_that("multiple .fam files warn and use the first", {
  dir <- make_admix_dir()
  # extra.fam sorts before test.fam, so it becomes the one used
  writeLines(c(
    "F1\tS1\t0\t0\t0\t1",
    "F2\tS2\t0\t0\t0\t1",
    "F3\tS3\t0\t0\t0\t1",
    "F4\tS4\t0\t0\t0\t1",
    "F5\tS5\t0\t0\t0\t1"
  ), file.path(dir, "extra.fam"))
  expect_warning(admixture_phylo_analysis(dir), "Multiple .fam files")
})

test_that("print method shows a summary", {
  dir <- make_admix_dir()
  res <- admixture_phylo_analysis(dir)
  expect_output(print(res), "Samples")
  expect_output(print(res), "K values")
})

test_that("input validation rejects bad arguments", {
  dir <- make_admix_dir()
  expect_error(admixture_phylo_analysis(file.path(dir, "nope")), "not found")
  expect_error(admixture_phylo_analysis(dir, fam_id_col = 0), "positive integer")
  expect_error(admixture_phylo_analysis(dir, tree_width = 1.5), "between 0 and 1")
  expect_error(admixture_phylo_analysis(dir, tree_linewidth = -1), "positive number")
  expect_error(admixture_phylo_analysis(dir, tip_point_size = 0), "positive number")
  expect_error(admixture_phylo_analysis(dir, color_scale = "nope"), "should be one of")
  expect_error(admixture_phylo_analysis(dir, save = 42), "NULL or a single")
  expect_error(admixture_phylo_analysis(dir, tree_file = "missing.nwk"),
               "phylo object or a path")
})
