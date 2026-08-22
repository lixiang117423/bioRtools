# Synthetic OcBSA window table: 3 chromosomes, 100 kb sliding windows.
# Chr2/Chr10 pair catches lexicographic (Chr10 < Chr2) ordering bugs.
make_ocbsa_data <- function() {
  chrs <- c("Chr1", "Chr2", "Chr10")
  chrlens <- c(1000000L, 800000L, 600000L)
  do.call(rbind, lapply(seq_along(chrs), function(i) {
    data.frame(
      Chr = chrs[[i]],
      Pos = seq(100000L, chrlens[[i]] - 100000L, by = 100000L),
      ChrLen = chrlens[[i]],
      VirtualIndex = 0.2,
      Pool1Index = 0.4,
      Pool2Index = 0.35,
      DeltaIndex = 0.05,
      MarkerNum = 50L,
      stringsAsFactors = FALSE
    )
  }))
}

layer_classes <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1], character(1))
}

test_that("plot_ocbsa returns a ggplot object", {
  p <- plot_ocbsa(make_ocbsa_data())
  expect_true(ggplot2::is_ggplot(p))
})

test_that("plot layers match the OcBSA style: boundaries, three index lines, delta points", {
  p <- plot_ocbsa(make_ocbsa_data())
  classes <- layer_classes(p)
  expect_equal(sum(classes == "GeomVline"), 1)
  expect_equal(sum(classes == "GeomLine"), 3)
  expect_equal(sum(classes == "GeomPoint"), 1)
})

test_that("chromosomes are ordered naturally, not lexicographically", {
  p <- plot_ocbsa(make_ocbsa_data())
  labels <- p$scales$get_scales("x")$labels
  expect_equal(as.character(labels), c("Chr1", "Chr2", "Chr10"))
})

test_that("x-axis breaks sit at chromosome midpoints on the cumulative axis", {
  p <- plot_ocbsa(make_ocbsa_data())
  breaks <- p$scales$get_scales("x")$breaks
  # mid = cum_start + chrlen / 2; chr1: 0 + 0.5M, chr2: 1M + 0.4M, chr10: 1.8M + 0.3M
  expect_equal(breaks, c(500000, 1400000, 2100000))
})

test_that("chr_sep strips the accession prefix from labels", {
  df <- make_ocbsa_data()
  df$Chr <- paste0("So120_", df$Chr)
  p <- plot_ocbsa(df, chr_sep = "_", chr_index = 2)
  labels <- p$scales$get_scales("x")$labels
  expect_equal(as.character(labels), c("Chr1", "Chr2", "Chr10"))
})

test_that("plot_ocbsa accepts a file path", {
  df <- make_ocbsa_data()
  path <- tempfile(fileext = ".tsv")
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  p <- plot_ocbsa(path)
  expect_true(ggplot2::is_ggplot(p))
})

test_that("plot_ocbsa validates inputs", {
  df <- make_ocbsa_data()

  expect_error(plot_ocbsa(list(a = 1)),
               "must be a file path \\(character\\) or data frame")
  expect_error(plot_ocbsa(df[, -which(names(df) == "DeltaIndex")]),
               "Required columns missing.*DeltaIndex")
  expect_error(plot_ocbsa(df, chr_sep = ""),
               "single non-empty string")
})
