#!/usr/bin/env Rscript
# Script to add @author tags to all functions missing them

author_line <- "#' @author Xiang LI \\email{lixiang117423@@foxmail.com}"

# Files that need @author tags
files_to_update <- c(
  "R/run_wgcna.R",
  "R/plot_multi_volcano.R",
  "R/identify_core_microbiome.R",
  "R/get_lm_stats.R",
  "R/academic_color_scales.R"
)

cat("Processing files...\n\n")

for (filepath in files_to_update) {
  if (!file.exists(filepath)) {
    cat("Skipping:", filepath, "(not found)\n")
    next
  }

  cat("Processing:", basename(filepath), "\n")

  # Read file
  lines <- readLines(filepath, warn = FALSE)

  # Find @export tags
  export_lines <- grep("^#' @export", lines)

  if (length(export_lines) == 0) {
    cat("  No @export found, skipping...\n")
    next
  }

  modified <- FALSE
  new_lines <- lines

  # Process each @export (go backwards to maintain line numbers)
  for (i in length(export_lines):1) {
    export_line <- export_lines[i]

    # Look back 50 lines for @author
    search_start <- max(1, export_line - 50)
    search_block <- new_lines[search_start:export_line]

    has_author <- any(grepl("^#' @author", search_block))

    if (!has_author) {
      # Look for @export and insert @author before it
      # Find the last roxygen comment before @export
      insert_pos <- export_line - 1

      # Skip empty lines
      while (insert_pos > 0 && (new_lines[insert_pos] == "" || !grepl("^#'", new_lines[insert_pos]))) {
        insert_pos <- insert_pos - 1
      }

      if (insert_pos > 0) {
        # Insert @author after the last roxygen comment
        insert_pos <- insert_pos + 1

        # Check if @author already exists at this position
        if (insert_pos <= length(new_lines) && !grepl("^#' @author", new_lines[insert_pos])) {
          new_lines <- append(new_lines, author_line, after = insert_pos - 1)
          modified <- TRUE
          cat("  Added @author at line", insert_pos, "\n")
        }
      }
    }
  }

  if (modified) {
    writeLines(new_lines, filepath)
    cat("  ✓ File updated\n")
  } else {
    cat("  No changes needed\n")
  }

  cat("\n")
}

cat("✓ Complete!\n")
