# Migrated from tests/test_plot_ideogram_synteny.R (bare script).
# Dual-species synteny (RIdeogram port).

build_synteny_ribbons <- bioRtools:::build_synteny_ribbons
ideo_synteny_row_layout <- bioRtools:::ideo_synteny_row_layout
ideo_cm_px <- bioRtools:::ideo_cm_px

is_ggplot <- function(x) inherits(x, "ggplot")

test_that("synteny plot returns a ggplot", {
  p <- plot_ideogram_synteny(df.ideo.synteny_karyotype, df.ideo.synteny)
  expect_true(is_ggplot(p))
})

test_that("chromosome row layout matches original RIdeogram formulas", {
  # From RIdeogram c6_synteny.svg, chr1 of species 1 (Grape):
  #   x = 124.015745, y = 88.582675, width = 26.240873
  total_cm <- 170 / 10
  wt <- function(n) (total_cm - 1.5 - (n - 1) * 0.1) * ideo_cm_px
  k1 <- df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Grape", ]
  rl <- ideo_synteny_row_layout(setNames(k1$End, k1$Chr), wt(nrow(k1)))
  expect_equal(unname(rl$starts[1]), 124.015745, tolerance = 1e-4)
  expect_equal(unname(rl$widths[1]), 26.240873, tolerance = 1e-4)
  expect_equal(2.5 * ideo_cm_px, 88.582675, tolerance = 1e-4)   # top-row y
})

test_that("species-label y matches original (Grape 103.299, Populus 298.181)", {
  h_outer <- 17.71654
  sp1_y <- 2.5 * ideo_cm_px + h_outer - 3
  sp2_y <- 8 * ideo_cm_px + h_outer - 3
  expect_equal(sp1_y, 103.29921, tolerance = 1e-4)
  expect_equal(sp2_y, 298.18110, tolerance = 1e-4)
})

test_that("ribbon polygon count equals synteny rows", {
  total_cm <- 170 / 10
  wt <- function(n) (total_cm - 1.5 - (n - 1) * 0.1) * ideo_cm_px
  k1 <- df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Grape", ]
  k2 <- df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Populus", ]
  ribbons <- build_synteny_ribbons(df.ideo.synteny,
    ideo_synteny_row_layout(setNames(k1$End, k1$Chr), wt(nrow(k1))),
    ideo_synteny_row_layout(setNames(k2$End, k2$Chr), wt(nrow(k2))),
    (2.5 + 0.5) * ideo_cm_px, (2.5 + 0.5 + 5) * ideo_cm_px)
  expect_equal(length(unique(ribbons$ribbon_id)), nrow(df.ideo.synteny))
})

test_that("validation — wrong species count / bad columns", {
  bad <- df.ideo.synteny_karyotype
  bad$species <- as.factor("Grape")
  expect_error(plot_ideogram_synteny(bad, df.ideo.synteny))
  expect_error(plot_ideogram_synteny(df.ideo.synteny_karyotype[, 1:3], df.ideo.synteny))
})
