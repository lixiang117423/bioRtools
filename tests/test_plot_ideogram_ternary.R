# Tests for plot_ideogram_ternary(): three-species triangular synteny (RIdeogram port).
source("R/utils_ideogram.R")
source("R/plot_ideogram_ternary.R")
suppressPackageStartupMessages({ library(ggplot2); library(pkgload) })
data(df.ideo.ternary_karyotype); data(df.ideo.ternary)

is_ggplot <- function(x) inherits(x, "ggplot")

cat("Test 1: ternary plot returns a ggplot\n")
p <- plot_ideogram_ternary(df.ideo.ternary_karyotype, df.ideo.ternary)
stopifnot(is_ggplot(p))
cat("✓ Test 1 passed\n\n")

cat("Test 2: species-1 chromosome geometry matches original RIdeogram\n")
# From RIdeogram c7_ternary.svg: chr1 species 1 x=106.29921, width=7.40977,
# baseline y=487.49987.
params <- ideo_ternary_params()
kat1 <- df.ideo.ternary_karyotype[df.ideo.ternary_karyotype$species == "Amborella", ]
lay1 <- ideo_ternary_species(kat1)
stopifnot(all.equal(unname(lay1$starts[1]), 106.29921, tolerance = 1e-4))
stopifnot(all.equal(unname(lay1$widths[1]),   7.40977, tolerance = 1e-4))
stopifnot(all.equal(params$baseline,          487.49987, tolerance = 1e-4))
cat("✓ Test 2 passed (species-1 chr1 matches original to 1e-4)\n\n")

cat("Test 3: species-2/3 anchors are rotated relative to the local frame\n")
layouts <- lapply(unique(df.ideo.ternary_karyotype$species), function(s)
  ideo_ternary_species(df.ideo.ternary_karyotype[df.ideo.ternary_karyotype$species == s, ]))
a2 <- ideo_ternary_anchor(2, 1, 0, layouts, params)   # species 2, chr1, pos 0
a3 <- ideo_ternary_anchor(3, 1, 0, layouts, params)
local2 <- c(layouts[[2]]$starts[1], params$baseline)
stopifnot(!isTRUE(all.equal(a2[1, ], local2, tolerance = 1e-3)))   # rotated -> different
stopifnot(a2[1, 2] < params$baseline)                              # rotated upward
stopifnot(a3[1, 2] < params$baseline)
cat("✓ Test 3 passed (species 2/3 anchors rotated up from baseline)\n\n")

cat("Test 4: gradient ribbons (fill = \"gradient\") render without error\n")
sg <- df.ideo.ternary; sg$fill[1:20] <- "gradient"
stopifnot(is_ggplot(plot_ideogram_ternary(df.ideo.ternary_karyotype, sg)))
cat("✓ Test 4 passed\n\n")

cat("Test 5: validation — wrong species count\n")
bad2 <- droplevels(df.ideo.ternary_karyotype[df.ideo.ternary_karyotype$species != "Liriodendron", ])
bad2$species <- as.character(bad2$species)
stopifnot(inherits(try(plot_ideogram_ternary(bad2, df.ideo.ternary), silent = TRUE), "try-error"))
cat("✓ Test 5 passed\n\n")

cat("All plot_ideogram_ternary tests done.\n")
