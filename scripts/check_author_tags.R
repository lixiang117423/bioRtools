#!/usr/bin/env Rscript
# Script to find functions missing @author tags

files <- list.files("R", pattern = "\\.R$", full.names = TRUE)

cat("Checking for missing @author tags...\n\n")

files_missing_author <- c()

for (file in files) {
  lines <- readLines(file, warn = FALSE)

  # Find @export tags
  export_lines <- grep("^#' @export", lines)

  if (length(export_lines) == 0) {
    next
  }

  # For each @export, check if there's an @author before it
  for (export_line in export_lines) {
    # Look back 100 lines for @author
    search_start <- max(1, export_line - 100)
    search_block <- lines[search_start:export_line]

    has_author <- any(grepl("^#' @author", search_block))

    if (!has_author) {
      # Try to find function name
      func_match <- grep("<- function\\(", lines[min(export_line + 1, length(lines))])
      if (length(func_match) > 0) {
        func_line <- export_line + 1
        func_name <- gsub("^\\s*(\\w+)\\s*<-.*", "\\1", lines[func_line])
        func_name <- gsub("^.*\\b(\\w+)\\s*<-.*", "\\1", func_name) # Fallback

        files_missing_author <- rbind(
          data.frame(
            file = basename(file),
            function_name = func_name,
            export_line = export_line,
            stringsAsFactors = FALSE
          ),
          if (exists("files_missing_author")) files_missing_author else data.frame()
        )
      }
    }
  }
}

cat("=== Summary ===\n")
cat("Functions missing @author:", nrow(files_missing_author), "\n\n")

if (nrow(files_missing_author) > 0) {
  cat("Functions to update:\n")
  for (i in 1:nrow(files_missing_author)) {
    cat(sprintf("  %s::%s() (line %d)\n",
                files_missing_author$file[i],
                files_missing_author$function_name[i],
                files_missing_author$export_line[i]))
  }
}

# Save to file for processing
if (nrow(files_missing_author) > 0) {
  write.csv(files_missing_author, "missing_author_tags.csv", row.names = FALSE)
  cat("\nSaved list to: missing_author_tags.csv\n")
}
