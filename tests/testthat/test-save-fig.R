# Migrated from tests/test_save_fig.R (bare script).

# raw magic-byte check
has_magic <- function(f, magic) {
  if (!file.exists(f)) return(FALSE)
  con <- file(f, "rb")
  on.exit(close(con), add = TRUE)
  identical(readBin(con, "raw", nchar(magic, type = "bytes")), charToRaw(magic))
}

make_plot <- function() {
  ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width, colour = Species)) +
    ggplot2::geom_point()
}

base <- file.path(tempdir(), "bioRtools_savefig_test")
cleanup <- function() unlink(c(
  paste0(base, c(".pdf", ".png", ".tiff", ".svg", ".eps")),
  tempfile("sf_", fileext = ".pdf")
))

test_that("single .pdf from extension", {
  cleanup()
  p <- make_plot()
  f <- tempfile("sf_single", fileext = ".pdf")
  save_fig(p, f)
  expect_true(file.exists(f))
  expect_true(has_magic(f, "%PDF"))
})

test_that("dual format (pdf+png) from bare base name", {
  cleanup()
  p <- make_plot()
  save_fig(p, base, formats = c("pdf", "png"))
  expect_true(file.exists(paste0(base, ".pdf")))
  expect_true(file.exists(paste0(base, ".png")))
  expect_true(has_magic(paste0(base, ".pdf"), "%PDF"))
  expect_true(has_magic(paste0(base, ".png"), "\x89PNG"))
})

test_that("filename with extension + formats strips extension", {
  cleanup()
  p <- make_plot()
  save_fig(p, paste0(base, ".pdf"), formats = c("pdf", "png"))
  expect_true(file.exists(paste0(base, ".pdf")))
  expect_true(file.exists(paste0(base, ".png")))
})

test_that("unsupported format errors", {
  cleanup()
  p <- make_plot()
  expect_error(save_fig(p, tempfile("sf_bad", fileext = ".xyz")))
})

test_that("no extension and formats = NULL errors", {
  cleanup()
  p <- make_plot()
  expect_error(save_fig(p, file.path(tempdir(), "sf_noext")))
})
