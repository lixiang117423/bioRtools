# Test script for the Wong (Okabe-Ito) colour-blind palette.

source("R/wong_palette.R")
source("R/utils_plotting.R")   # for create_discrete_scale()

suppressPackageStartupMessages({
  library(ggplot2)
})

wong_expected <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# --- Test 1: full palette --------------------------------------------------
cat("Test 1: wong_palette() returns the 8 Wong colours in order\n")
stopifnot(identical(wong_palette(), wong_expected))
cat("✓ Test 1 passed — 8 hex codes in Wong order\n\n")

# --- Test 2: subsetting + recycling ----------------------------------------
cat("Test 2: subsetting (n<=8) and recycling (n>8)\n")
stopifnot(identical(wong_palette(3), wong_expected[1:3]))
stopifnot(identical(wong_palette(0), character(0)))
rec <- wong_palette(10)
stopifnot(length(rec) == 10)
stopifnot(identical(rec, rep(wong_expected, length.out = 10)))
cat("✓ Test 2 passed — n=3 subset, n=0 empty, n=10 recycled\n\n")

# --- Test 3: scale palette resolves to wong colours ------------------------
cat("Test 3: scale_*_wong() palette resolves to wong colours\n")
sc <- scale_fill_wong()
stopifnot(inherits(sc, "ScaleDiscrete"))
stopifnot(identical(sc$palette(3), wong_palette(3)))
stopifnot(identical(scale_colour_wong()$palette(4), wong_palette(4)))
cat("✓ Test 3 passed — fill + colour scales map to wong_palette\n\n")

# --- Test 4: real plot renders + built fills match -------------------------
cat("Test 4: render a 3-group plot with scale_fill_wong\n")
p <- ggplot(iris, aes(Species, Sepal.Length, fill = Species)) +
  geom_boxplot() +
  scale_fill_wong()
gb <- ggplot_build(p)
stopifnot(setequal(unique(gb$data[[1]]$fill), wong_palette(3)))
grob <- ggplotGrob(p)
stopifnot(inherits(grob, "gtable"))
cat("✓ Test 4 passed — 3 boxplot fills are exactly wong_palette(3)\n\n")

# --- Test 5: create_discrete_scale supports palette = "wong" ---------------
cat("Test 5: create_discrete_scale(n, 'wong') uses wong colours\n")
s <- create_discrete_scale(5, "wong")
stopifnot(inherits(s, "ScaleDiscrete"))
stopifnot(setequal(s$palette(5), wong_palette(5)))
s2 <- create_discrete_scale(3, "wong", aesthetic = "fill")
stopifnot(setequal(s2$palette(3), wong_palette(3)))
cat("✓ Test 5 passed — internal helper branches to wong_palette\n\n")

cat("=====================================\n")
cat("All wong_palette tests passed! ✓\n")
cat("=====================================\n")
