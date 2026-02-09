#!/usr/bin/env Rscript
# Script to find examples that need library() calls

# Find all @examples
examples <- list.files("R", pattern = "\\.R$", full.names = TRUE)

cat("Checking examples for library() calls...\n\n")

files_needing_library <- c()

for (file in examples) {
  lines <- readLines(file, warn = FALSE)
  example_start <- grep("^#' @examples", lines)

  if (length(example_start) > 0) {
    # Get next 15 lines after @examples
    example_block <- lines[min(example_start + 1) : min(example_start + 15, length(lines))]

    # Check if library() is present
    has_library <- any(grepl("^#' library\\(", example_block))

    if (!has_library) {
      files_needing_library <- c(files_needing_library, basename(file))
      cat("Missing library():", basename(file), "\n")
    }
  }
}

cat("\n=== Summary ===\n")
cat("Total files needing library():", length(files_needing_library), "\n")

if (length(files_needing_library) > 0) {
  cat("\nFiles to update:\n")
  for (f in files_needing_library) {
    cat("  -", f, "\n")
  }
}
