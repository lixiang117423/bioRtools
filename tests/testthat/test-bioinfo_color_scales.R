# Anchor colors mirrored from the palette definitions in R/bioinfo_color_scales.R.
# Categorical assertions are exact (direct subsetting, no interpolation).
# Ramp endpoints use a tolerance: colorRampPalette(space = "Lab") round-trips
# anchors through Lab space, perturbing endpoints by a few 1/255 steps.
bioinfo_deep_first3 <- c("#497EAD", "#B85449", "#479E8D")
bioinfo_deep_last <- "#855D49"
bioinfo_light_11th <- "#B1D3F0"
bioinfo_light_last <- "#CCB0A3"

near_color <- function(actual, expected, tol = 2) {
  delta <- max(abs(grDevices::col2rgb(actual) - grDevices::col2rgb(expected)))
  if (delta > tol) {
    paste0(actual, " vs ", expected, " (max channel delta ", delta, ")")
  } else {
    TRUE
  }
}

test_that("get_cat_colors() returns deep tones first for n <= 10", {
  expect_identical(get_cat_colors(3), bioinfo_deep_first3)
  expect_identical(get_cat_colors(1), "#497EAD")
  cols10 <- get_cat_colors(10)
  expect_length(cols10, 10)
  expect_identical(cols10[10], bioinfo_deep_last)
})

test_that("get_cat_colors() interleaves light tones for 10 < n <= 20", {
  cols11 <- get_cat_colors(11)
  expect_length(cols11, 11)
  expect_identical(cols11[1], bioinfo_deep_first3[1])
  # deep/light pairs stay adjacent per hue family: slot 2 is blue's light tone
  expect_identical(cols11[2], bioinfo_light_11th)
  expect_identical(cols11[11], "#D17997")
  cols20 <- get_cat_colors(20)
  expect_length(cols20, 20)
  expect_identical(cols20[20], bioinfo_light_last)
  expect_true(!anyDuplicated(cols20))
})

test_that("get_cat_colors() interpolates beyond 20", {
  cols25 <- get_cat_colors(25)
  expect_length(cols25, 25)
  expect_true(near_color(cols25[1], bioinfo_deep_first3[1]))
  expect_true(near_color(cols25[25], bioinfo_light_last))
})

test_that("get_cat_colors() handles edge and invalid n", {
  expect_identical(get_cat_colors(0), character(0))
  expect_error(get_cat_colors("a"), "'n' must be")
  expect_error(get_cat_colors(2.5), "'n' must be")
  expect_error(get_cat_colors(c(1, 2)), "'n' must be")
})

test_that("get_cat_colors() output is unnamed", {
  expect_null(names(get_cat_colors(5)))
})

test_that("get_seq_colors() ramps between the family anchors", {
  cols <- get_seq_colors("blue")
  expect_length(cols, 100)
  expect_true(near_color(cols[1], "#EBF3FA"))
  expect_true(near_color(cols[100], "#2F5A80"))
  two <- get_seq_colors("red", 2)
  expect_true(near_color(two[1], "#FAECEB"))
  expect_true(near_color(two[2], "#8A3930"))
})

test_that("get_seq_colors() supports all ten families", {
  families <- c("blue", "red", "teal", "purple", "orange", "pink",
                "gold", "navy", "olive", "brown")
  for (fam in families) {
    expect_length(get_seq_colors(fam, 10), 10)
  }
})

test_that("get_seq_colors() rejects unknown family names", {
  expect_error(get_seq_colors("chartreuse"), "blue, red, teal")
})

test_that("ramp getters reject invalid n", {
  expect_error(get_seq_colors("blue", 0), "'n' must be")
  expect_error(get_seq_colors("blue", 2.5), "'n' must be")
  expect_error(get_div_colors(1), "'n' must be")
  expect_error(get_multiseq_colors("ten"), "'n' must be")
})

test_that("get_div_colors() returns the requested variant", {
  soft <- get_div_colors(100)
  expect_length(soft, 100)
  expect_true(near_color(soft[1], "#2C5F8A"))
  expect_true(near_color(soft[100], "#B23A3A"))
  strong <- get_div_colors(5, palette = "strong")
  expect_true(near_color(strong[1], "#08519C"))
  expect_true(near_color(strong[5], "#B2182B"))
  expect_error(get_div_colors(10, palette = "neon"), "soft")
})

test_that("get_multiseq_colors() ramps through the multi-hue anchors", {
  cols <- get_multiseq_colors()
  expect_length(cols, 100)
  expect_true(near_color(cols[1], "#3A1F5D"))
  expect_true(near_color(cols[100], "#F2E94E"))
  four <- get_multiseq_colors(4)
  for (i in 1:4) {
    expect_true(near_color(four[i], c("#3A1F5D", "#2E6E8E", "#3FA796", "#F2E94E")[i]))
  }
})

# --- ggplot2 scale wiring ----------------------------------------------------

make_grp_plot <- function() {
  df <- data.frame(x = 1:6, y = 1:6, grp = rep(c("a", "b", "c"), each = 2))
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = grp, fill = grp)) +
    ggplot2::geom_point(shape = 21, size = 3)
}

make_con_plot <- function() {
  df <- data.frame(x = 1:10, y = 1:10, z = seq(-1, 1, length.out = 10))
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = z, fill = z)) +
    ggplot2::geom_point(shape = 21, size = 3)
}

test_that("scale_color_bioinfo() maps groups onto the categorical palette", {
  s <- scale_color_bioinfo()
  expect_s3_class(s, "ScaleDiscrete")
  expect_identical(s$palette(3), get_cat_colors(3))
  expect_identical(s$palette(12), get_cat_colors(12))

  built <- ggplot2::ggplot_build(make_grp_plot() + scale_color_bioinfo())
  expect_length(unique(built$data[[1]]$colour), 3)
})

test_that("scale_colour_bioinfo() is an alias of scale_color_bioinfo()", {
  expect_identical(scale_colour_bioinfo, scale_color_bioinfo)
})

test_that("scale_fill_bioinfo() targets the fill aesthetic", {
  s <- scale_fill_bioinfo()
  expect_identical(s$aesthetics, "fill")
  expect_identical(s$palette(2), get_cat_colors(2))
})

test_that("discrete bioinfo scales apply alpha", {
  s <- scale_fill_bioinfo(alpha = 0.5)
  cols <- s$palette(1)
  expect_match(cols, "^#[0-9A-F]{6}80$")
})

test_that("scale_color_bioinfo_seq() maps [0,1] onto the family ramp", {
  s <- scale_color_bioinfo_seq(name = "red")
  expect_s3_class(s, "ScaleContinuous")
  expect_identical(s$palette(0), get_seq_colors("red", 256)[1])
  expect_identical(s$palette(1), get_seq_colors("red", 256)[256])

  rev_s <- scale_color_bioinfo_seq(name = "red", reverse = TRUE)
  expect_identical(rev_s$palette(0), get_seq_colors("red", 256)[256])
})

test_that("continuous bioinfo scales build real plots", {
  expect_no_error(ggplot2::ggplot_build(
    make_con_plot() +
      scale_color_bioinfo_seq(name = "teal") +
      scale_fill_bioinfo_seq(name = "teal", alpha = 0.6)
  ))
  expect_error(ggplot2::ggplot_build(
    make_con_plot() + scale_color_bioinfo_seq(name = "nope")
  ), "blue, red, teal")
})

test_that("scale_fill_bioinfo_div() maps [0,1] onto the diverging ramp", {
  s <- scale_fill_bioinfo_div()
  expect_s3_class(s, "ScaleContinuous")
  expect_identical(s$palette(0), get_div_colors(256)[1])
  expect_identical(s$palette(1), get_div_colors(256)[256])

  strong_s <- scale_fill_bioinfo_div(palette = "strong")
  expect_identical(strong_s$palette(0), get_div_colors(256L, "strong")[1])

  expect_no_error(ggplot2::ggplot_build(
    make_con_plot() +
      scale_color_bioinfo_div() +
      scale_fill_bioinfo_div(palette = "strong", reverse = TRUE)
  ))
})

test_that("bioinfo scales reject alpha outside (0, 1]", {
  expect_error(scale_fill_bioinfo(alpha = 0), "'alpha' must be")
  expect_error(scale_color_bioinfo(alpha = 1.5), "'alpha' must be")
  expect_error(scale_color_bioinfo_seq(alpha = -1), "'alpha' must be")
  expect_error(scale_fill_bioinfo_div(alpha = 2), "'alpha' must be")
})

test_that("scale_colour_bioinfo_seq() and scale_colour_bioinfo_div() are aliases", {
  expect_identical(scale_colour_bioinfo_seq, scale_color_bioinfo_seq)
  expect_identical(scale_colour_bioinfo_div, scale_color_bioinfo_div)
})

test_that("preview_bioinfo_palettes() draws without opening a new device", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(preview_bioinfo_palettes())
  expect_no_error(preview_bioinfo_palettes(family = "seq"))
})
