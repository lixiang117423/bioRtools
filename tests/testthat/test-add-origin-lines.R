# Migrated from tests/test_add_origin_lines.R (bare script).
# add_origin_lines() reference-line helper.

base_plot <- ggplot2::ggplot(datasets::iris,
                             ggplot2::aes(Sepal.Length, Sepal.Width)) +
  ggplot2::geom_point()

test_that("returns a list of two layers, vertical then horizontal", {
  layers <- add_origin_lines()

  expect_type(layers, "list")
  expect_length(layers, 2)
  expect_s3_class(layers[[1]]$geom, "GeomVline")
  expect_s3_class(layers[[2]]$geom, "GeomHline")
})

test_that("defaults render identically to bare geom_vline + geom_hline", {
  p_ref <- base_plot +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  gb_ref <- ggplot2::ggplot_build(p_ref)

  p <- base_plot + add_origin_lines()
  expect_length(p$layers, 3)

  gb <- ggplot2::ggplot_build(p)
  compare_cols <- c("xintercept", "yintercept", "colour", "linewidth",
                    "linetype", "alpha")
  for (col in compare_cols) {
    expect_equal(gb$data[[2]][[col]], gb_ref$data[[2]][[col]],
                 tolerance = 1.5e-8)
    expect_equal(gb$data[[3]][[col]], gb_ref$data[[3]][[col]],
                 tolerance = 1.5e-8)
  }
})

test_that("custom intercepts and style overrides", {
  p2 <- base_plot + add_origin_lines(
    x_intercept = 1, y_intercept = -2,
    color = "red", linewidth = 2
  )
  gb2 <- ggplot2::ggplot_build(p2)

  expect_true(all(gb2$data[[2]]$xintercept == 1))
  expect_true(all(gb2$data[[3]]$yintercept == -2))
  expect_true(all(gb2$data[[2]]$colour == "red"))
  expect_true(all(gb2$data[[3]]$linewidth == 2))
})

test_that("extra arguments via ... reach both layers", {
  layers4 <- add_origin_lines(show.legend = FALSE)

  expect_false(layers4[[1]]$show.legend)
  expect_false(layers4[[2]]$show.legend)
})

test_that("input validation rejects bad arguments", {
  expect_error(add_origin_lines(x_intercept = "a"))
  expect_error(add_origin_lines(y_intercept = c(0, 1)))
  expect_error(add_origin_lines(alpha = 2))
  expect_error(add_origin_lines(linewidth = -1))
  expect_error(add_origin_lines(linetype = c("dashed", "solid")))
  expect_error(add_origin_lines(color = NA_character_))
  expect_error(add_origin_lines(alpha = 2), "'alpha'")
})

test_that("plot renders", {
  p <- base_plot + add_origin_lines()
  grob <- ggplot2::ggplotGrob(p)

  expect_s3_class(grob, "gtable")
})
