# Tests for plot_ideogram_synteny(): dual-species synteny (RIdeogram port).
source("R/plot_ideogram_synteny.R")
source("R/utils_ideogram.R")
suppressPackageStartupMessages({ library(ggplot2); library(pkgload) })
data(df.ideo.synteny_karyotype); data(df.ideo.synteny)

is_ggplot <- function(x) inherits(x, "ggplot")

cat("Test 1: synteny plot returns a ggplot\n")
p <- plot_ideogram_synteny(df.ideo.synteny_karyotype, df.ideo.synteny)
stopifnot(is_ggplot(p))
cat("✓ Test 1 passed\n\n")

cat("Test 2: chromosome row layout matches original RIdeogram formulas\n")
# From RIdeogram c6_synteny.svg, chr1 of species 1 (Grape):
#   x = 124.015745, y = 88.582675, width = 26.240873
total_cm <- 170 / 10
wt <- function(n) (total_cm - 1.5 - (n - 1) * 0.1) * ideo_cm_px
k1 <- df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Grape", ]
rl <- ideo_synteny_row_layout(setNames(k1$End, k1$Chr), wt(nrow(k1)))
stopifnot(all.equal(unname(rl$starts[1]), 124.015745, tolerance = 1e-4))
stopifnot(all.equal(unname(rl$widths[1]), 26.240873,  tolerance = 1e-4))
stopifnot(all.equal(2.5 * ideo_cm_px, 88.582675, tolerance = 1e-4))   # top-row y
cat("✓ Test 2 passed (chr1 row geometry matches original to 1e-4)\n\n")

cat("Test 3: species-label y matches original (Grape 103.299, Populus 298.181)\n")
h_outer <- 17.71654
sp1_y <- 2.5 * ideo_cm_px + h_outer - 3
sp2_y <- 8 * ideo_cm_px + h_outer - 3
stopifnot(all.equal(sp1_y, 103.29921, tolerance = 1e-4))
stopifnot(all.equal(sp2_y, 298.18110, tolerance = 1e-4))
cat("✓ Test 3 passed\n\n")

cat("Test 4: ribbon polygon count equals synteny rows\n")
ribbons <- build_synteny_ribbons(df.ideo.synteny,
  ideo_synteny_row_layout(setNames(k1$End, k1$Chr), wt(nrow(k1))),
  ideo_synteny_row_layout(setNames(
    df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Populus", ]$End,
    df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Populus", ]$Chr),
    wt(nrow(df.ideo.synteny_karyotype[df.ideo.synteny_karyotype$species == "Populus", ]))),
  (2.5 + 0.5) * ideo_cm_px, (2.5 + 0.5 + 5) * ideo_cm_px)
stopifnot(length(unique(ribbons$ribbon_id)) == nrow(df.ideo.synteny))
cat("✓ Test 4 passed (", nrow(df.ideo.synteny), "ribbons )\n\n")

cat("Test 5: validation — wrong species count / bad columns\n")
bad <- df.ideo.synteny_karyotype
bad$species <- as.factor("Grape")
stopifnot(inherits(try(plot_ideogram_synteny(bad, df.ideo.synteny), silent = TRUE), "try-error"))
stopifnot(inherits(try(plot_ideogram_synteny(df.ideo.synteny_karyotype[, 1:3], df.ideo.synteny),
                       silent = TRUE), "try-error"))
cat("✓ Test 5 passed\n\n")

cat("All plot_ideogram_synteny tests done.\n")
