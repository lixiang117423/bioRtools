# Migrated from tests/test_plot_ideogram.R (bare script).
# Single-species idiogram (RIdeogram port): validates return type, that every
# label mode runs, and that the chromosome layout matches the original
# RIdeogram formulas to full precision.

ideo_layout_chr <- bioRtools:::ideo_layout_chr

is_ggplot <- function(x) inherits(x, "ggplot")

test_that("karyotype-only returns a ggplot", {
  p <- plot_ideogram(df.ideo.karyotype)
  expect_true(is_ggplot(p))
})

test_that("overlaid heatmap + marker + heatmap/line/polygon all run", {
  skip_if_namespace_unavailable("ggnewscale")

  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density)))
  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density,
                                      label = df.ideo.rna_marker, label_type = "marker")))
  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density,
                                      label = df.ideo.ltr_density, label_type = "heatmap")))
  ld <- aggregate(Value ~ Chr + Start + End, df.ideo.gene_density, mean)
  ld$color <- "4575b4"
  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, label = ld, label_type = "line")))
  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, label = ld, label_type = "polygon")))
})

test_that("chromosome layout matches original RIdeogram formulas", {
  # human chr1 is the longest (End == max), at index 1. From RIdeogram c1.svg:
  #   x_left = 70.86614, x_right = 80.51938, y_top = 93.40930, cap r = 4.82662
  lay <- ideo_layout_chr(df.ideo.karyotype, width = 170)
  r1 <- lay[lay$Chr == 1, ]
  expect_equal(r1$x_left, 70.86614, tolerance = 1e-4)
  expect_equal(r1$x_right, 80.51938, tolerance = 1e-4)
  expect_equal(r1$y_top, 93.409295, tolerance = 1e-4)
  expect_equal(r1$chr_width / 2, 4.82662, tolerance = 1e-4)
})

test_that("3-column karyotype (no centromere) also works", {
  kat3 <- df.ideo.karyotype[, c("Chr", "Start", "End")]
  expect_true(is_ggplot(plot_ideogram(kat3, overlaid = df.ideo.gene_density)))
  expect_null(ideo_layout_chr(kat3, 170)$y_ce_start)
})

test_that("input validation", {
  expect_error(plot_ideogram(data.frame()))
  expect_error(plot_ideogram(data.frame(Chr = "1")))
})

test_that("line / polygon dual-series (7-column Value_1/Value_2)", {
  ld2 <- aggregate(Value ~ Chr + Start + End, df.ideo.gene_density, mean)
  ld2$Value_2 <- rev(ld2$Value)                       # a second series
  ld2$Value_1 <- ld2$Value
  ld2$color_1 <- "4575b4"; ld2$color_2 <- "d73027"
  ld2$Value <- NULL
  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, label = ld2, label_type = "line")))
  expect_true(is_ggplot(plot_ideogram(df.ideo.karyotype, label = ld2, label_type = "polygon")))
})

test_that("Chr-mismatch warning fires for unknown chromosomes", {
  bad_overlaid <- df.ideo.gene_density
  bad_overlaid$Chr[1] <- "ZZ"
  w <- tryCatch(plot_ideogram(df.ideo.karyotype, overlaid = bad_overlaid),
                warning = function(w) w$message)
  expect_true(is.character(w) && grepl("not in 'karyotype\\$Chr'", w))
})
