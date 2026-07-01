# Test script for theme_prism() gray_background option.

source("R/theme_prism.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

gray_fill <- ggplot2::theme_gray()$panel.background$fill

# --- Test 1: gray_background = TRUE uses theme_gray() colors ----------------
cat("Test 1: gray_background theme elements\n")
th <- theme_prism(gray_background = TRUE)
stopifnot(identical(as.character(th$panel.background$fill), as.character(gray_fill)))
stopifnot(identical(as.character(th$plot.background$fill), "white"))
stopifnot(identical(as.character(th$panel.grid.major$colour), "white"))
stopifnot(inherits(th$panel.grid.minor, "element_blank"))
stopifnot(inherits(th$panel.border, "element_blank"))
cat("✓ Test 1 passed — grey panel (#EBEBEB), white grid lines, no border\n\n")

# --- Test 2: default look unchanged -----------------------------------------
cat("Test 2: default (gray_background = FALSE) unchanged\n")
th0 <- theme_prism()
stopifnot(is.na(th0$panel.background$fill))
stopifnot(inherits(th0$panel.grid, "element_blank"))
cat("✓ Test 2 passed — clean Prism panel preserved\n\n")

# --- Test 3: invalid value errors -------------------------------------------
cat("Test 3: invalid gray_background errors\n")
tryCatch({
  theme_prism(gray_background = "nope")
  cat("✗ Test 3 failed: should have errored\n")
}, error = function(e) {
  cat("✓ Test 3 passed:", conditionMessage(e), "\n\n")
})

# --- Test 4: renders on a real ggplot ---------------------------------------
cat("Test 4: render a plot with gray_background\n")
p <- ggplot(iris, aes(Species, Sepal.Length, fill = Species)) +
  geom_boxplot() +
  theme_prism(gray_background = TRUE)
grob <- ggplot2::ggplotGrob(p)   # errors if the theme is malformed
stopifnot(inherits(grob, "gtable"))
cat("✓ Test 4 passed — plot builds without error\n\n")

cat("=====================================\n")
cat("All theme_prism tests passed! ✓\n")
cat("=====================================\n")
