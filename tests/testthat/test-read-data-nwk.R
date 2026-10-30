if (!can_load_namespace("treeio")) {
  test_that("optional dependency integration tests are skipped", {
    testthat::skip("treeio is not available")
  })
} else {
test_that("read_data reads Newick trees as phylo objects", {
  dir <- tempfile("nwk")
  dir.create(dir)
  f <- file.path(dir, "tree.nwk")
  writeLines("((A:1,B:1):1,C:2);", f)

  tr <- read_data(f)
  expect_s3_class(tr, "phylo")
  expect_equal(tr$tip.label, c("A", "B", "C"))
  expect_equal(tr$Nnode, 2L)
})

test_that("read_data supports common tree file extensions", {
  dir <- tempfile("nwk2")
  dir.create(dir)
  for (ext in c("newick", "tree", "treefile")) {
    f <- file.path(dir, paste0("tree.", ext))
    writeLines("((X:1,Y:1):1,Z:2);", f)
    tr <- read_data(f)
    expect_s3_class(tr, "phylo")
    expect_equal(tr$tip.label, c("X", "Y", "Z"))
  }
})

test_that("read_data forwards extra arguments to read.tree", {
  dir <- tempfile("nwk3")
  dir.create(dir)
  f <- file.path(dir, "tree.nwk")
  # node.labels stored as node.label
  writeLines("((A:1,B:1)n1:1,C:2);", f)
  tr <- read_data(f)
  # ape assigns "" to the root node, so the labelled internal node is second
  expect_equal(tr$node.label, c("", "n1"))
})
}
