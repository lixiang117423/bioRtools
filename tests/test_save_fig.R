# Test script for save_fig() journal export helper.

source("R/save_fig.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

p <- ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
  geom_point()

# raw magic-byte check
has_magic <- function(f, magic) {
  if (!file.exists(f)) return(FALSE)
  con <- file(f, "rb"); on.exit(close(con), add = TRUE)
  identical(readBin(con, "raw", nchar(magic, type = "bytes")), charToRaw(magic))
}

base <- file.path(tempdir(), "bioRtools_savefig_test")
cleanup <- function() unlink(c(
  paste0(base, c(".pdf", ".png", ".tiff", ".svg", ".eps")),
  tempfile("sf_", fileext = ".pdf")
))
cleanup()

# --- Test 1: single format inferred from extension ------------------------
cat("Test 1: single .pdf from extension\n")
f <- tempfile("sf_single", fileext = ".pdf")
save_fig(p, f)
stopifnot(file.exists(f))
stopifnot(has_magic(f, "%PDF"))
cat("✓ Test 1 passed — valid PDF written\n\n")

# --- Test 2: dual format, bare base name -----------------------------------
cat("Test 2: dual format (pdf+png) from bare base name\n")
cleanup()
save_fig(p, base, formats = c("pdf", "png"))
stopifnot(file.exists(paste0(base, ".pdf")))
stopifnot(file.exists(paste0(base, ".png")))
stopifnot(has_magic(paste0(base, ".pdf"), "%PDF"))
stopifnot(has_magic(paste0(base, ".png"), "\x89PNG"))
cat("✓ Test 2 passed — base.pdf + base.png, correct magic\n\n")

# --- Test 3: filename WITH extension + formats -> extension stripped -------
cat("Test 3: filename with extension + formats strips extension\n")
cleanup()
save_fig(p, paste0(base, ".pdf"), formats = c("pdf", "png"))
stopifnot(file.exists(paste0(base, ".pdf")))
stopifnot(file.exists(paste0(base, ".png")))
cat("✓ Test 3 passed — extension stripped, both formats written\n\n")

# --- Test 4: unsupported format errors -------------------------------------
cat("Test 4: unsupported format errors\n")
res4 <- tryCatch({
  save_fig(p, tempfile("sf_bad", fileext = ".xyz"))
  "no-error"
}, error = function(e) conditionMessage(e))
stopifnot(!identical(res4, "no-error"))
cat("✓ Test 4 passed:", res4, "\n\n")

# --- Test 5: no extension + formats = NULL errors --------------------------
cat("Test 5: no extension and formats = NULL errors\n")
res5 <- tryCatch({
  save_fig(p, file.path(tempdir(), "sf_noext"))
  "no-error"
}, error = function(e) conditionMessage(e))
stopifnot(!identical(res5, "no-error"))
cat("✓ Test 5 passed:", res5, "\n\n")

cleanup()
cat("=====================================\n")
cat("All save_fig tests passed! ✓\n")
cat("=====================================\n")
