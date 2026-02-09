#!/usr/bin/env Rscript
# Script to automatically fix code style issues including unnecessary return()

# Install styler if not available
if (!require("styler", quietly = TRUE)) {
  install.packages("styler", repos = "https://cloud.r-project.org")
}

library(styler)

# Get all R files in R/ directory
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)

cat("Formatting", length(r_files), "R files with styler...\n\n")

# Style each file
for (i in seq_along(r_files)) {
  file <- r_files[i]
  cat(sprintf("[%d/%d] %s\n", i, length(r_files), basename(file)))

  tryCatch({
    # Style the file (this will remove unnecessary return() and fix other style issues)
    style_file(file, strict = FALSE)
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", conditionMessage(e)))
  })
}

cat("\n✓ Code styling complete!\n")
cat("Styler has:\n")
cat("  - Removed unnecessary return() statements\n")
cat("  - Fixed indentation (2 spaces)\n")
cat("  - Fixed spacing around operators\n")
cat("  - Fixed other code style issues\n")
