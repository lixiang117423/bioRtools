# Migrated from tests/test_plot_ideogram_ternary.R (bare script).
# Three-species triangular synteny (RIdeogram port).

ideo_ternary_params <- bioRtools:::ideo_ternary_params
ideo_ternary_species <- bioRtools:::ideo_ternary_species
ideo_ternary_anchor <- bioRtools:::ideo_ternary_anchor

is_ggplot <- function(x) inherits(x, "ggplot")

test_that("ternary plot returns a ggplot", {
  p <- plot_ideogram_ternary(df.ideo.ternary_karyotype, df.ideo.ternary)
  expect_true(is_ggplot(p))
})

test_that("species-1 chromosome geometry matches original RIdeogram", {
  # From RIdeogram c7_ternary.svg: chr1 species 1 x=106.29921, width=7.40977,
  # baseline y=487.49987.
  params <- ideo_ternary_params()
  kat1 <- df.ideo.ternary_karyotype[df.ideo.ternary_karyotype$species == "Amborella", ]
  lay1 <- ideo_ternary_species(kat1)
  expect_equal(unname(lay1$starts[1]), 106.29921, tolerance = 1e-4)
  expect_equal(unname(lay1$widths[1]), 7.40977, tolerance = 1e-4)
  expect_equal(params$baseline, 487.49987, tolerance = 1e-4)
})

test_that("species-2/3 anchors are rotated relative to the local frame", {
  params <- ideo_ternary_params()
  layouts <- lapply(unique(df.ideo.ternary_karyotype$species), function(s)
    ideo_ternary_species(df.ideo.ternary_karyotype[df.ideo.ternary_karyotype$species == s, ]))
  a2 <- ideo_ternary_anchor(2, 1, 0, layouts, params)   # species 2, chr1, pos 0
  a3 <- ideo_ternary_anchor(3, 1, 0, layouts, params)
  local2 <- c(layouts[[2]]$starts[1], params$baseline)
  expect_false(isTRUE(all.equal(a2[1, ], local2, tolerance = 1e-3)))   # rotated -> different
  expect_lt(a2[1, 2], params$baseline)                                 # rotated upward
  expect_lt(a3[1, 2], params$baseline)
})

test_that("gradient ribbons (fill = 'gradient') render without error", {
  sg <- df.ideo.ternary
  sg$fill[1:20] <- "gradient"
  expect_true(is_ggplot(plot_ideogram_ternary(df.ideo.ternary_karyotype, sg)))
})

test_that("validation — wrong species count", {
  bad2 <- droplevels(df.ideo.ternary_karyotype[df.ideo.ternary_karyotype$species != "Liriodendron", ])
  bad2$species <- as.character(bad2$species)
  expect_error(plot_ideogram_ternary(bad2, df.ideo.ternary))
})
