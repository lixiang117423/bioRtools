# Migrated from tests/test_plot_composition_alluvium.R (bare script).
# plot_composition_alluvium() with synthetic wide-format composition data, so
# it runs in R CMD check without external files.

# Synthetic wide-format composition data (taxon column + sample columns)
set.seed(1)
test_data <- data.frame(
  "Viral Taxonomy (family)" = paste0("Family_", 1:18),
  Sample1 = sample(0:100, 18, replace = TRUE),
  Sample2 = sample(0:100, 18, replace = TRUE),
  Sample3 = sample(0:100, 18, replace = TRUE),
  Sample4 = sample(0:100, 18, replace = TRUE),
  check.names = FALSE
)

make_alluvium_plot <- function(...) {
  plot_composition_alluvium(test_data,
                            taxon_column = "Viral Taxonomy (family)", ...)
}

test_that("basic usage returns a ggplot object", {
  p1 <- make_alluvium_plot()
  expect_true(ggplot2::is_ggplot(p1))
})

test_that("custom color palette", {
  p2 <- make_alluvium_plot(palette = "Set2")
  expect_true(ggplot2::is_ggplot(p2))
})

test_that("custom alpha", {
  p3 <- make_alluvium_plot(alpha = 0.9)
  expect_true(ggplot2::is_ggplot(p3))
})

test_that("custom y-axis formatting", {
  p4 <- make_alluvium_plot(y_accuracy = 0.1)
  expect_true(ggplot2::is_ggplot(p4))
})

test_that("custom legend position", {
  p5 <- make_alluvium_plot(legend_position = "bottom")
  expect_true(ggplot2::is_ggplot(p5))
})

test_that("id_column is used when taxon_column is not found", {
  p6 <- plot_composition_alluvium(
    test_data,
    taxon_column = "family",
    id_column = "Viral Taxonomy (family)"
  )
  expect_true(ggplot2::is_ggplot(p6))
})

test_that("plot can be saved to a file", {
  p1 <- make_alluvium_plot()
  out_file <- tempfile(fileext = ".png")
  # RColorBrewer warns when 18 taxa exceed the palette size while rendering
  suppressWarnings(
    ggplot2::ggsave(out_file, p1, width = 8, height = 6, dpi = 100)
  )
  expect_true(file.exists(out_file))
})
