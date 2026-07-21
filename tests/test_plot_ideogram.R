# Tests for plot_ideogram(): single-species idiogram (RIdeogram port).
# Validates return type, that every label mode runs, and that the chromosome
# layout matches the original RIdeogram formulas to full precision.
source("R/plot_ideogram.R")
source("R/utils_ideogram.R")
suppressPackageStartupMessages({ library(ggplot2); library(pkgload) })
data(df.ideo.karyotype); data(df.ideo.gene_density)
data(df.ideo.ltr_density); data(df.ideo.rna_marker)

is_ggplot <- function(x) inherits(x, "ggplot")

cat("Test 1: karyotype-only returns a ggplot\n")
p <- plot_ideogram(df.ideo.karyotype)
stopifnot(is_ggplot(p))
cat("✓ Test 1 passed\n\n")

cat("Test 2: overlaid heatmap + marker + heatmap/line/polygon all run\n")
stopifnot(is_ggplot(plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density)))
stopifnot(is_ggplot(plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density,
                                  label = df.ideo.rna_marker, label_type = "marker")))
stopifnot(is_ggplot(plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density,
                                  label = df.ideo.ltr_density, label_type = "heatmap")))
ld <- aggregate(Value ~ Chr + Start + End, df.ideo.gene_density, mean); ld$color <- "4575b4"
stopifnot(is_ggplot(plot_ideogram(df.ideo.karyotype, label = ld, label_type = "line")))
stopifnot(is_ggplot(plot_ideogram(df.ideo.karyotype, label = ld, label_type = "polygon")))
cat("✓ Test 2 passed\n\n")

cat("Test 3: chromosome layout matches original RIdeogram formulas\n")
# human chr1 is the longest (End == max), at index 1. From RIdeogram c1.svg:
#   x_left = 70.86614, x_right = 80.51938, y_top = 93.40930, cap r = 4.82662
lay <- ideo_layout_chr(df.ideo.karyotype, width = 170)
r1 <- lay[lay$Chr == 1, ]
stopifnot(all.equal(r1$x_left,   70.86614,  tolerance = 1e-4))
stopifnot(all.equal(r1$x_right,  80.51938,  tolerance = 1e-4))
stopifnot(all.equal(r1$y_top,    93.409295, tolerance = 1e-4))
stopifnot(all.equal(r1$chr_width / 2, 4.82662, tolerance = 1e-4))
cat("✓ Test 3 passed (chr1 geometry matches original to 1e-4)\n\n")

cat("Test 4: 3-column karyotype (no centromere) also works\n")
kat3 <- df.ideo.karyotype[, c("Chr", "Start", "End")]
stopifnot(is_ggplot(plot_ideogram(kat3, overlaid = df.ideo.gene_density)))
stopifnot(is.null(ideo_layout_chr(kat3, 170)$y_ce_start))
cat("✓ Test 4 passed\n\n")

cat("Test 5: input validation\n")
stopifnot(inherits(try(plot_ideogram(data.frame()), silent = TRUE), "try-error"))
stopifnot(inherits(try(plot_ideogram(data.frame(Chr = "1")), silent = TRUE), "try-error"))
cat("✓ Test 5 passed\n\n")

cat("All plot_ideogram tests done.\n")
