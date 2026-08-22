# Dimensions are read straight from the PNG IHDR header so the test needs
# no image package: width/height live in bytes 17-24 of the file.
png_dims <- function(path) {
  raw <- readBin(path, "raw", n = 24)
  c(
    width  = sum(as.integer(raw[17:20]) * 256^(3:0)),
    height = sum(as.integer(raw[21:24]) * 256^(3:0))
  )
}

expected_px <- function(cm) cm / 2.54 * 600

test_that("write_data() derives the missing image side from the golden ratio", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
  f <- tempfile(fileext = ".png")

  expect_px <- function(actual, cm) {
    expect_true(abs(actual - expected_px(cm)) <= 2,
      info = paste0(actual, "px vs expected ", expected_px(cm), "px (", cm, "cm)"))
  }

  # width only -> height = width / 1.618
  write_data(p, f, width = 10)
  d <- png_dims(f)
  expect_px(d[["width"]], 10)
  expect_px(d[["height"]], 10 / 1.618)

  # height only -> width = height * 1.618
  write_data(p, f, height = 5)
  d <- png_dims(f)
  expect_px(d[["width"]], 5 * 1.618)
  expect_px(d[["height"]], 5)

  # neither side -> 8.5 x 6 cm default
  write_data(p, f)
  d <- png_dims(f)
  expect_px(d[["width"]], 8.5)
  expect_px(d[["height"]], 6)

  # both sides -> used as-is (not golden)
  write_data(p, f, width = 4, height = 4)
  d <- png_dims(f)
  expect_px(d[["width"]], 4)
  expect_px(d[["height"]], 4)
})

test_that("write_data() table branches ignore image dimensions", {
  df <- data.frame(a = 1:3, b = c("x", "y", "z"))
  f <- tempfile(fileext = ".csv")
  write_data(df, f, width = 10)
  expect_equal(read.csv(f), df)
})
