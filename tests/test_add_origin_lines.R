# Test script for add_origin_lines() reference-line helper.

source("R/add_origin_lines.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

base_plot <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
  geom_point()

# --- Test 1: returns two layers, vertical then horizontal ----------------
cat("Test 1: returns a list of two layers, vertical then horizontal\n")
layers <- add_origin_lines()
stopifnot(is.list(layers), length(layers) == 2)
stopifnot(inherits(layers[[1]]$geom, "GeomVline"))
stopifnot(inherits(layers[[2]]$geom, "GeomHline"))
cat("✓ Test 1 passed — GeomVline + GeomHline, in order\n\n")

# --- Test 2: defaults render exactly like bare geom_vline + geom_hline -----
cat("Test 2: defaults render identically to bare geom_vline + geom_hline\n")
p_ref <- base_plot +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
gb_ref <- ggplot2::ggplot_build(p_ref)
p <- base_plot + add_origin_lines()
stopifnot(length(p$layers) == 3)
gb <- ggplot2::ggplot_build(p)
compare_cols <- c("xintercept", "yintercept", "colour", "linewidth", "linetype", "alpha")
for (col in compare_cols) {
  stopifnot(isTRUE(all.equal(gb$data[[2]][[col]], gb_ref$data[[2]][[col]])))
  stopifnot(isTRUE(all.equal(gb$data[[3]][[col]], gb_ref$data[[3]][[col]])))
}
cat("✓ Test 2 passed — built layers identical to bare geoms\n\n")

# --- Test 3: custom intercepts and overrides -------------------------------
cat("Test 3: custom intercepts and style overrides\n")
p2 <- base_plot + add_origin_lines(
  x_intercept = 1, y_intercept = -2,
  color = "red", linewidth = 2
)
gb2 <- ggplot_build(p2)
stopifnot(isTRUE(all(gb2$data[[2]]$xintercept == 1)))
stopifnot(isTRUE(all(gb2$data[[3]]$yintercept == -2)))
stopifnot(isTRUE(all(gb2$data[[2]]$colour == "red")))
stopifnot(isTRUE(all(gb2$data[[3]]$linewidth == 2)))
cat("✓ Test 3 passed — intercepts 1/-2, colour red, linewidth 2\n\n")

# --- Test 4: ... passes extra arguments to both geoms ----------------------
cat("Test 4: extra arguments via ... reach both layers\n")
layers4 <- add_origin_lines(show.legend = FALSE)
stopifnot(isTRUE(layers4[[1]]$show.legend == FALSE))
stopifnot(isTRUE(layers4[[2]]$show.legend == FALSE))
cat("✓ Test 4 passed — show.legend = FALSE forwarded to both geoms\n\n")

# --- Test 5: input validation -----------------------------------------------
cat("Test 5: input validation rejects bad arguments\n")
stopifnot(inherits(try(add_origin_lines(x_intercept = "a"), silent = TRUE), "try-error"))
stopifnot(inherits(try(add_origin_lines(y_intercept = c(0, 1)), silent = TRUE), "try-error"))
stopifnot(inherits(try(add_origin_lines(alpha = 2), silent = TRUE), "try-error"))
stopifnot(inherits(try(add_origin_lines(linewidth = -1), silent = TRUE), "try-error"))
stopifnot(inherits(try(add_origin_lines(linetype = c("dashed", "solid")), silent = TRUE), "try-error"))
stopifnot(inherits(try(add_origin_lines(color = NA_character_), silent = TRUE), "try-error"))
err_msg <- tryCatch(add_origin_lines(alpha = 2), error = function(e) conditionMessage(e))
stopifnot(grepl("'alpha'", err_msg))
cat("✓ Test 5 passed — bad intercepts/style rejected with actionable message\n\n")

# --- Test 6: renders to a grob ----------------------------------------------
cat("Test 6: plot renders\n")
grob <- ggplotGrob(p)
stopifnot(inherits(grob, "gtable"))
cat("✓ Test 6 passed — ggplotGrob returns a gtable\n\n")

cat("=====================================\n")
cat("All add_origin_lines tests passed! ✓\n")
cat("=====================================\n")
