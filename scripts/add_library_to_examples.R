#!/usr/bin/env Rscript
# Script to automatically add library() calls to examples

files_needing_library <- c(
  "admixture_phylo_analysis.R",
  "calc_expression_delta_ct.R",
  "calc_expression_delta_delta_ct.R",
  "calc_expression_qpcr_efficiency.R",
  "calc_expression_standard_curve.R",
  "calc_standard_curve.R",
  "data.R",
  "df2fasta.R",
  "extract_tree_hierarchy.R",
  "fasta2df.R",
  "get_lm_stats.R",
  "get_motif_from_meme.R",
  "identify_core_microbiome.R",
  "label_signif.R",
  "ormalize_int.R",
  "pav_gwas.R",
  "plot_LDheatmap.R",
  "plot_motif_location.R",
  "plot_pfam.R",
  "row_stat.R",
  "run_wgcna.R",
  "scale01.R"
)

cat("Processing", length(files_needing_library), "files...\n\n")

for (filename in files_needing_library) {
  filepath <- file.path("R", filename)

  cat("Processing:", filename, "\n")

  # Read file
  lines <- readLines(filepath, warn = FALSE)

  # Find @examples line
  example_idx <- grep("^#' @examples", lines)

  if (length(example_idx) == 0) {
    cat("  No @examples found, skipping...\n")
    next
  }

  # Get the example block (next 20 lines)
  example_block <- lines[min(example_idx + 1) : min(example_idx + 20, length(lines))]

  # Check if already has library()
  if (any(grepl("^#' library\\(", example_block))) {
    cat("  Already has library(), skipping...\n")
    next
  }

  # Extract packages used in examples (look for ::)
  package_calls <- grep("^#'.*::", example_block, value = TRUE)

  # Extract package names
  packages <- unique(
    gsub("^#'.*?([a-zA-Z][a-zA-Z0-9._]*)::.*", "\\1", package_calls)
  )

  # Remove common base R packages
  base_packages <- c("stats", "base", "utils", "graphics", "methods", "grDevices")
  packages <- setdiff(packages, base_packages)

  # Add bioRtools if not present
  if (!"bioRtools" %in% packages) {
    packages <- c("bioRtools", packages)
  }

  # Build library() lines
  if (length(packages) > 0) {
    library_lines <- paste0("#' library(", packages, ")")
    library_text <- paste(library_lines, collapse = "\n")

    # Insert library calls after @examples
    insert_pos <- example_idx[1] + 1

    # Check if next line is empty or starts with \dontrun
    next_line <- lines[insert_pos]

    if (grepl("^#'(\\\\dontrun|\\s*$)", next_line)) {
      # Insert before the empty line or \dontrun
      new_lines <- c(
        lines[1:(example_idx[1])],
        library_lines,
        lines[(example_idx[1] + 1):length(lines)]
      )
    } else {
      # Insert directly after @examples
      new_lines <- c(
        lines[1:example_idx[1]],
        library_lines,
        lines[(example_idx[1] + 1):length(lines)]
      )
    }

    # Write back
    writeLines(new_lines, filepath)
    cat("  Added library() for:", paste(packages, collapse = ", "), "\n")
  } else {
    cat("  No packages to add\n")
  }

  cat("\n")
}

cat("✓ Complete!\n")
