# Migrated from tests/test_wong_palette.R (bare script).

create_discrete_scale <- bioRtools:::create_discrete_scale

wong_expected <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

test_that("wong_palette() returns the 8 Wong colours in order", {
  expect_identical(wong_palette(), wong_expected)
})

test_that("subsetting (n<=8) and recycling (n>8)", {
  expect_identical(wong_palette(3), wong_expected[1:3])
  expect_identical(wong_palette(0), character(0))
  rec <- wong_palette(10)
  expect_length(rec, 10)
  expect_identical(rec, rep(wong_expected, length.out = 10))
})

test_that("scale_*_wong() palette resolves to wong colours", {
  sc <- scale_fill_wong()
  expect_s3_class(sc, "ScaleDiscrete")
  expect_identical(sc$palette(3), wong_palette(3))
  expect_identical(scale_colour_wong()$palette(4), wong_palette(4))
})

test_that("render a 3-group plot with scale_fill_wong", {
  p <- ggplot2::ggplot(iris, ggplot2::aes(Species, Sepal.Length, fill = Species)) +
    ggplot2::geom_boxplot() +
    scale_fill_wong()
  gb <- ggplot2::ggplot_build(p)
  expect_setequal(unique(gb$data[[1]]$fill), wong_palette(3))
  grob <- ggplot2::ggplotGrob(p)
  expect_s3_class(grob, "gtable")
})

test_that("create_discrete_scale(n, 'wong') uses wong colours", {
  s <- create_discrete_scale(5, "wong")
  expect_s3_class(s, "ScaleDiscrete")
  expect_setequal(s$palette(5), wong_palette(5))
  s2 <- create_discrete_scale(3, "wong", aesthetic = "fill")
  expect_setequal(s2$palette(3), wong_palette(3))
})
