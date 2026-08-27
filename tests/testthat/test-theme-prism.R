# Migrated from tests/test_theme_prism.R (bare script).

gray_fill <- ggplot2::theme_gray()$panel.background$fill

test_that("gray_background theme elements", {
  th <- theme_prism(gray_background = TRUE)
  expect_identical(as.character(th$panel.background$fill), as.character(gray_fill))
  expect_identical(as.character(th$plot.background$fill), "white")
  expect_identical(as.character(th$panel.grid.major$colour), "white")
  expect_s3_class(th$panel.grid.minor, "element_blank")
  expect_s3_class(th$panel.border, "element_blank")
})

test_that("default (gray_background = FALSE) unchanged", {
  th0 <- theme_prism()
  expect_true(is.na(th0$panel.background$fill))
  expect_s3_class(th0$panel.grid, "element_blank")
})

test_that("invalid gray_background errors", {
  expect_error(theme_prism(gray_background = "nope"))
})

test_that("render a plot with gray_background", {
  p <- ggplot2::ggplot(iris, ggplot2::aes(Species, Sepal.Length, fill = Species)) +
    ggplot2::geom_boxplot() +
    theme_prism(gray_background = TRUE)
  grob <- ggplot2::ggplotGrob(p)
  expect_s3_class(grob, "gtable")
})

test_that("light = TRUE (default) — plain text, thinner lines", {
  th_light <- theme_prism()
  expect_identical(th_light$text$face, "plain")
  expect_equal(th_light$line$linewidth, 14 / 22)
  expect_equal(th_light$rect$linewidth, 14 / 22)
})

test_that("light = FALSE — reproduces old bold + base_size/14", {
  th_bold <- theme_prism(light = FALSE)
  expect_identical(th_bold$text$face, "bold")
  expect_equal(th_bold$line$linewidth, 14 / 14)
  expect_equal(th_bold$rect$linewidth, 14 / 14)
})

test_that("plot.tag face bold in both light modes", {
  th_light <- theme_prism()
  th_bold <- theme_prism(light = FALSE)
  expect_identical(th_light$plot.tag$face, "bold")
  expect_identical(th_bold$plot.tag$face, "bold")
})

test_that("strip.text face bold in both modes", {
  th_light <- theme_prism()
  th_bold <- theme_prism(light = FALSE)
  expect_identical(th_light$strip.text$face, "bold")
  expect_identical(th_bold$strip.text$face, "bold")
})

test_that("invalid light errors", {
  expect_error(theme_prism(light = "nope"))
})
