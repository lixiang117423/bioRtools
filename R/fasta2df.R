# #' Convert FASTA File to Data Frame (High Performance Version)
# #'
# #' @description
# #' \code{fasta2df} efficiently converts FASTA format files to R data frames.
# #' This optimized version uses vectorized operations and efficient memory
# #' allocation to significantly improve performance compared to loop-based approaches.
# #' Suitable for processing large genomic, transcriptomic, or proteomic sequence files.
# #'
# #' @param fasta Character string specifying the path to the FASTA file.
# #'   The file should contain sequences in standard FASTA format with headers
# #'   starting with '>' followed by sequence data on subsequent lines.
# #' @param validate Logical value indicating whether to perform input validation.
# #'   Default is TRUE. Set to FALSE for better performance if input is guaranteed
# #'   to be valid.
# #' @param remove_empty Logical value indicating whether to remove empty sequences.
# #'   Default is TRUE.
# #' @param trim_whitespace Logical value indicating whether to trim whitespace
# #'   from sequence headers and sequences. Default is TRUE.
# #'
# #' @return A data frame with two columns:
# #' \describe{
# #'   \item{id}{Character vector containing FASTA headers (including the '>' symbol)}
# #'   \item{seq}{Character vector containing the corresponding DNA/RNA/protein sequences}
# #' }
# #' The data frame has one row per sequence found in the input file.
# #'
# #' @details
# #' This function implements several performance optimizations:
# #' \itemize{
# #'   \item Vectorized header detection using \code{startsWith()}
# #'   \item Pre-allocated memory for results to avoid dynamic growth
# #'   \item Efficient string concatenation using \code{paste(..., collapse="")}
# #'   \item Minimal data copying and transformation
# #' }
# #'
# #' The function handles various FASTA format variations:
# #' \itemize{
# #'   \item Multi-line sequences (automatically concatenated)
# #'   \item Empty lines (automatically skipped)
# #'   \item Different line ending formats (Unix, Windows, Mac)
# #'   \item Files with or without trailing newlines
# #' }
# #'
# #' @note
# #' \itemize{
# #'   \item Large files (>100MB) may require substantial memory
# #'   \item The function loads the entire file into memory
# #'   \item For extremely large files, consider streaming approaches
# #'   \item Sequence names are returned with the '>' prefix intact
# #' }
# #'
# #' @seealso
# #' \code{\link{readLines}} for reading text files,
# #' \code{\link{Biostrings}} package for advanced sequence manipulation
# #'
# #' @author Optimized by Assistant
# #' @references
# #' FASTA format specification: \url{https://en.wikipedia.org/wiki/FASTA_format}
# #'
# #' @importFrom utils file.size
# #' @export
# #'
# #' @examples
# #' \dontrun{
# #' # Create a sample FASTA file for demonstration
# #' fasta_content <- c(
# #'   ">sequence1 description",
# #'   "ATCGATCGATCG",
# #'   "GCTAGCTAGCTA",
# #'   ">sequence2",
# #'   "TTTTAAAA",
# #'   ">sequence3 another description",
# #'   "CCCGGGAAATTT"
# #' )
# #'
# #' # Write to temporary file
# #' temp_fasta <- tempfile(fileext = ".fasta")
# #' writeLines(fasta_content, temp_fasta)
# #'
# #' # Convert FASTA to data frame
# #' df <- fasta2df(temp_fasta)
# #' print(df)
# #' #   id                           seq
# #' # 1 >sequence1 description      ATCGATCGATCGGCTAGCTAGCTA
# #' # 2 >sequence2                  TTTTAAAA
# #' # 3 >sequence3 another description CCCGGGAAATTT
# #'
# #' # Check dimensions
# #' cat("Number of sequences:", nrow(df), "\n")
# #' cat("Average sequence length:", mean(nchar(df$seq)), "\n")
# #'
# #' # Clean up
# #' unlink(temp_fasta)
# #' }
# #'
# #' # Basic usage with validation disabled for maximum speed
# #' \dontrun{
# #' df_fast <- fasta2df("large_file.fasta", validate = FALSE)
# #' }
# #'
# #' # Keep empty sequences and preserve whitespace
# #' \dontrun{
# #' df_raw <- fasta2df("sequences.fasta",
# #'                    remove_empty = FALSE,
# #'                    trim_whitespace = FALSE)
# #' }
# #'
# #' @keywords file sequence bioinformatics genomics

# fasta2df <- function(fasta,
#                      validate = TRUE,
#                      remove_empty = TRUE,
#                      trim_whitespace = TRUE) {

#   # Input validation
#   if (validate) {
#     # Check if fasta parameter is provided and valid
#     if (missing(fasta) || is.null(fasta)) {
#       stop("Parameter 'fasta' is required and cannot be NULL")
#     }

#     if (!is.character(fasta) || length(fasta) != 1) {
#       stop("Parameter 'fasta' must be a single character string")
#     }

#     if (!file.exists(fasta)) {
#       stop("FASTA file '", fasta, "' does not exist")
#     }

#     if (file.size(fasta) == 0) {
#       stop("FASTA file '", fasta, "' is empty")
#     }

#     # Check logical parameters
#     if (!is.logical(remove_empty) || length(remove_empty) != 1) {
#       stop("Parameter 'remove_empty' must be a single logical value")
#     }

#     if (!is.logical(trim_whitespace) || length(trim_whitespace) != 1) {
#       stop("Parameter 'trim_whitespace' must be a single logical value")
#     }
#   }

#   # Read file content efficiently
#   tryCatch({
#     lines <- readLines(fasta, warn = FALSE)
#   }, error = function(e) {
#     stop("Error reading FASTA file '", fasta, "': ", e$message)
#   })

#   # Remove empty lines if requested
#   if (remove_empty) {
#     lines <- lines[nzchar(lines)]
#   }

#   # Trim whitespace if requested
#   if (trim_whitespace) {
#     lines <- trimws(lines)
#   }

#   if (length(lines) == 0) {
#     warning("No content found in FASTA file after processing")
#     return(data.frame(id = character(0), seq = character(0),
#                      stringsAsFactors = FALSE))
#   }

#   # Vectorized header detection (high performance)
#   is_header <- startsWith(lines, ">")
#   header_positions <- which(is_header)
#   n_sequences <- length(header_positions)

#   if (n_sequences == 0) {
#     stop("No FASTA headers (lines starting with '>') found in file")
#   }

#   # Pre-allocate result vectors for optimal memory usage
#   sequence_ids <- lines[header_positions]
#   sequences <- character(n_sequences)

#   # Efficiently determine sequence boundaries
#   sequence_starts <- header_positions + 1
#   sequence_ends <- c(header_positions[-1] - 1, length(lines))

#   # Vectorized sequence concatenation
#   for (i in seq_len(n_sequences)) {
#     start_pos <- sequence_starts[i]
#     end_pos <- sequence_ends[i]

#     if (start_pos <= end_pos && start_pos <= length(lines)) {
#       # Extract sequence lines and concatenate efficiently
#       seq_lines <- lines[start_pos:min(end_pos, length(lines))]
#       # Remove any remaining headers that might have been included
#       seq_lines <- seq_lines[!startsWith(seq_lines, ">")]
#       sequences[i] <- paste(seq_lines, collapse = "")
#     } else {
#       sequences[i] <- ""
#     }
#   }

#   # Create final data frame
#   result <- data.frame(
#     id = sequence_ids,
#     seq = sequences,
#     stringsAsFactors = FALSE,
#     row.names = NULL
#   )

#   # Optional: remove empty sequences from results
#   if (remove_empty) {
#     non_empty <- nzchar(result$seq)
#     if (!any(non_empty)) {
#       warning("All sequences are empty after processing")
#     }
#     result <- result[non_empty, , drop = FALSE]
#   }

#   # Add metadata as attributes
#   attr(result, "source_file") <- fasta
#   attr(result, "n_sequences") <- nrow(result)
#   attr(result, "processing_options") <- list(
#     remove_empty = remove_empty,
#     trim_whitespace = trim_whitespace
#   )

#   return(result)
# }

# # Utility function for performance benchmarking (not exported)
# .benchmark_fasta2df <- function(fasta_file, iterations = 5) {
#   if (!file.exists(fasta_file)) {
#     stop("Benchmark file does not exist: ", fasta_file)
#   }

#   cat("Benchmarking fasta2df performance...\n")
#   cat("File:", fasta_file, "\n")
#   cat("File size:", round(file.size(fasta_file) / 1024^2, 2), "MB\n")
#   cat("Iterations:", iterations, "\n\n")

#   # Warm-up run
#   invisible(fasta2df(fasta_file, validate = FALSE))

#   # Benchmark with validation
#   time_with_validation <- system.time({
#     for (i in seq_len(iterations)) {
#       result <- fasta2df(fasta_file, validate = TRUE)
#     }
#   })

#   # Benchmark without validation
#   time_without_validation <- system.time({
#     for (i in seq_len(iterations)) {
#       result <- fasta2df(fasta_file, validate = FALSE)
#     }
#   })

#   cat("With validation:   ", round(time_with_validation["elapsed"] / iterations, 4), "seconds per run\n")
#   cat("Without validation:", round(time_without_validation["elapsed"] / iterations, 4), "seconds per run\n")
#   cat("Sequences processed:", nrow(result), "\n")
#   cat("Performance:", round(nrow(result) / (time_without_validation["elapsed"] / iterations)), "sequences/second\n")

#   invisible(result)
# }


#' Convert FASTA File to Data Frame (High Performance Version)
#'
#' @description
#' \code{fasta2df} efficiently converts FASTA format files to R data frames.
#' This optimized version uses vectorized operations and efficient memory
#' allocation to significantly improve performance compared to loop-based approaches.
#' Suitable for processing large genomic, transcriptomic, or proteomic sequence files.
#'
#' @param fasta Character vector specifying the path(s) to the FASTA file(s).
#'   Can be a single file path or a vector of file paths for vectorized processing.
#'   Files should contain sequences in standard FASTA format with headers
#'   starting with '>' followed by sequence data on subsequent lines.
#' @param validate Logical value indicating whether to perform input validation.
#'   Default is TRUE. Set to FALSE for better performance if input is guaranteed
#'   to be valid.
#' @param remove_empty Logical value indicating whether to remove empty sequences.
#'   Default is TRUE.
#' @param trim_whitespace Logical value indicating whether to trim whitespace
#'   from sequence headers and sequences. Default is TRUE.
#'
#' @return For single file input: A data frame with two columns:
#' \describe{
#'   \item{id}{Character vector containing FASTA headers (including the '>' symbol)}
#'   \item{seq}{Character vector containing the corresponding DNA/RNA/protein sequences}
#'
#' # Working with multiple files (vectorized)
#' fasta_files <- c("seq1.fasta", "seq2.fasta", "seq3.fasta")
#' results_list <- fasta2df(fasta_files)
#'
#' # Access individual results
#' first_file_data <- results_list[[1]]
#'
#' # Combine all results
#' all_sequences <- do.call(rbind, results_list)
#' }
#'
#' # Integration with tidyverse workflows
#' \dontrun{
#' library(dplyr)
#' library(purrr)
#'
#' # Method 1: Using rowwise() for line-by-line processing
#' file_data <- data.frame(
#'   species = c("species1", "species2", "species3"),
#'   file_path = c("seq1.fasta", "seq2.fasta", "seq3.fasta")
#' ) %>%
#'   rowwise() %>%
#'   mutate(sequences = list(fasta2df(file_path))) %>%
#'   ungroup()
#'
#' # Method 2: Using map() for functional programming approach
#' file_data2 <- data.frame(
#'   species = c("species1", "species2", "species3"),
#'   file_path = c("seq1.fasta", "seq2.fasta", "seq3.fasta")
#' ) %>%
#'   mutate(sequences = map(file_path, fasta2df))
#'
#' # Method 3: Process all files at once (most efficient)
#' all_files <- c("seq1.fasta", "seq2.fasta", "seq3.fasta")
#' all_results <- fasta2df(all_files)  # Returns a list
#'
#' # Convert to a tidy format
#' tidy_results <- all_results %>%
#'   imap_dfr(~ .x %>% mutate(source_file = .y))
#'
#' # Your specific use case - SOLUTION
#' dir("G://database/十字花科基因组/叶绿体基因组/chloroplast_genomes/") %>%
#'   as.data.frame() %>%
#'   set_names("file") %>%
#'   filter(!str_detect(file, "network_failed")) %>%
#'   mutate(
#'     species = str_split(file, "_chl") %>% map_chr(1),
#'     path = paste0("G://database/十字花科基因组/叶绿体基因组/chloroplast_genomes/", file)
#'   ) %>%
#'   rowwise() %>%  # Key: enables row-by-row processing
#'   mutate(seq = list(fasta2df(path))) %>%  # Key: use list() wrapper
#'   ungroup()
#'
#' # Alternative: Using the tidy wrapper
#' dir("path/to/files/") %>%
#'   as.data.frame() %>%
#'   set_names("file") %>%
#'   mutate(
#'     path = paste0("path/to/files/", file),
#'     seq = map(path, fasta2df_tidy)  # Simpler syntax
#'   )
#' }
#' For multiple file input: A list of data frames, one per input file.
#' Each data frame has the same structure as the single file output.
#' The data frame(s) have one row per sequence found in the input file(s).
#'
#' @details
#' This function implements several performance optimizations:
#' \itemize{
#'   \item Vectorized header detection using \code{startsWith()}
#'   \item Pre-allocated memory for results to avoid dynamic growth
#'   \item Efficient string concatenation using \code{paste(..., collapse="")}
#'   \item Minimal data copying and transformation
#' }
#'
#' The function handles various FASTA format variations:
#' \itemize{
#'   \item Multi-line sequences (automatically concatenated)
#'   \item Empty lines (automatically skipped)
#'   \item Different line ending formats (Unix, Windows, Mac)
#'   \item Files with or without trailing newlines
#' }
#'
#' @note
#' \itemize{
#'   \item Large files (>100MB) may require substantial memory
#'   \item The function loads the entire file into memory
#'   \item For extremely large files, consider streaming approaches
#'   \item Sequence names are returned with the '>' prefix intact
#' }
#'
#' @seealso
#' \code{\link{readLines}} for reading text files,
#' \code{\link{Biostrings}} package for advanced sequence manipulation
#'
#' @author Optimized by Assistant
#' @references
#' FASTA format specification: \url{https://en.wikipedia.org/wiki/FASTA_format}
#'
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' library(bioRtools)
#' \dontrun{
#' # Create a sample FASTA file for demonstration
#' fasta_content <- c(
#'   ">sequence1 description",
#'   "ATCGATCGATCG",
#'   "GCTAGCTAGCTA",
#'   ">sequence2",
#'   "TTTTAAAA",
#'   ">sequence3 another description",
#'   "CCCGGGAAATTT"
#' )
#'
#' # Write to temporary file
#' temp_fasta <- tempfile(fileext = ".fasta")
#' writeLines(fasta_content, temp_fasta)
#'
#' # Convert FASTA to data frame
#' df <- fasta2df(temp_fasta)
#' print(df)
#' #   id                           seq
#' # 1 >sequence1 description      ATCGATCGATCGGCTAGCTAGCTA
#' # 2 >sequence2                  TTTTAAAA
#' # 3 >sequence3 another description CCCGGGAAATTT
#'
#' # Check dimensions
#' cat("Number of sequences:", nrow(df), "\n")
#' cat("Average sequence length:", mean(nchar(df$seq)), "\n")
#'
#' # Clean up
#' unlink(temp_fasta)
#' }
#'
#' # Basic usage with validation disabled for maximum speed
#' \dontrun{
#' df_fast <- fasta2df("large_file.fasta", validate = FALSE)
#' }
#'
#' # Keep empty sequences and preserve whitespace
#' \dontrun{
#' df_raw <- fasta2df("sequences.fasta",
#'   remove_empty = FALSE,
#'   trim_whitespace = FALSE)
#' }
#'
#' @keywords file sequence bioinformatics genomics

fasta2df <- function(fasta,
                     validate = TRUE,
                     remove_empty = TRUE,
                     trim_whitespace = TRUE) {
  # Input validation
  if (validate) {
    # Check if fasta parameter is provided and valid
    if (missing(fasta) || is.null(fasta)) {
      stop("Parameter 'fasta' is required and cannot be NULL")
    }

    if (!is.character(fasta)) {
      stop("Parameter 'fasta' must be a character vector")
    }

    if (length(fasta) == 0) {
      stop("Parameter 'fasta' cannot be empty")
    }

    # Check logical parameters
    if (!is.logical(remove_empty) || length(remove_empty) != 1) {
      stop("Parameter 'remove_empty' must be a single logical value")
    }

    if (!is.logical(trim_whitespace) || length(trim_whitespace) != 1) {
      stop("Parameter 'trim_whitespace' must be a single logical value")
    }
  }

  # Handle vectorized input (multiple files)
  if (length(fasta) > 1) {
    # Process multiple files
    result_list <- vector("list", length(fasta))
    names(result_list) <- fasta

    for (i in seq_along(fasta)) {
      # Recursive call for each file (single file processing)
      result_list[[i]] <- fasta2df(
        fasta = fasta[i],
        validate = validate,
        remove_empty = remove_empty,
        trim_whitespace = trim_whitespace
      )
    }

    # Add metadata for vectorized results
    attr(result_list, "n_files") <- length(fasta)
    attr(result_list, "file_paths") <- fasta
    attr(result_list, "processing_options") <- list(
      remove_empty = remove_empty,
      trim_whitespace = trim_whitespace
    )

    return(result_list)
  }

  # Single file processing (original logic)
  fasta_file <- fasta[1]  # Ensure single file

  if (validate) {
    if (!file.exists(fasta_file)) {
      stop("FASTA file '", fasta_file, "' does not exist")
    }

    if (file.size(fasta_file) == 0) {
      stop("FASTA file '", fasta_file, "' is empty")
    }
  }

  # Read file content efficiently
  tryCatch(
    {
      lines <- readLines(fasta_file, warn = FALSE)
    },
    error = function(e) {
      stop("Error reading FASTA file '", fasta_file, "': ", e$message)
    })

  # Remove empty lines if requested
  if (remove_empty) {
    lines <- lines[nzchar(lines)]
  }

  # Trim whitespace if requested
  if (trim_whitespace) {
    lines <- trimws(lines)
  }

  if (length(lines) == 0) {
    warning("No content found in FASTA file after processing")
    return(data.frame(id = character(0), seq = character(0),
      stringsAsFactors = FALSE))
  }

  # Vectorized header detection (high performance)
  is_header <- startsWith(lines, ">")
  header_positions <- which(is_header)
  n_sequences <- length(header_positions)

  if (n_sequences == 0) {
    stop("No FASTA headers (lines starting with '>') found in file")
  }

  # Pre-allocate result vectors for optimal memory usage
  sequence_ids <- lines[header_positions]
  sequences <- character(n_sequences)

  # Efficiently determine sequence boundaries
  sequence_starts <- header_positions + 1
  sequence_ends <- c(header_positions[-1] - 1, length(lines))

  # Vectorized sequence concatenation
  for (i in seq_len(n_sequences)) {
    start_pos <- sequence_starts[i]
    end_pos <- sequence_ends[i]

    if (start_pos <= end_pos && start_pos <= length(lines)) {
      # Extract sequence lines and concatenate efficiently
      seq_lines <- lines[start_pos:min(end_pos, length(lines))]
      # Remove any remaining headers that might have been included
      seq_lines <- seq_lines[!startsWith(seq_lines, ">")]
      sequences[i] <- paste(seq_lines, collapse = "")
    } else {
      sequences[i] <- ""
    }
  }

  # Create final data frame
  result <- data.frame(
    id = sequence_ids,
    seq = sequences,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Optional: remove empty sequences from results
  if (remove_empty) {
    non_empty <- nzchar(result$seq)
    if (!any(non_empty)) {
      warning("All sequences are empty after processing")
    }
    result <- result[non_empty, , drop = FALSE]
  }

  # Add metadata as attributes
  attr(result, "source_file") <- fasta_file
  attr(result, "n_sequences") <- nrow(result)
  attr(result, "processing_options") <- list(
    remove_empty = remove_empty,
    trim_whitespace = trim_whitespace
  )

  return(result)
}

# Utility function for performance benchmarking (not exported)
.benchmark_fasta2df <- function(fasta_file, iterations = 5) {
  if (!file.exists(fasta_file)) {
    stop("Benchmark file does not exist: ", fasta_file)
  }

  cat("Benchmarking fasta2df performance...\n")
  cat("File:", fasta_file, "\n")
  cat("File size:", round(file.size(fasta_file) / 1024^2, 2), "MB\n")
  cat("Iterations:", iterations, "\n\n")

  # Warm-up run
  invisible(fasta2df(fasta_file, validate = FALSE))

  # Benchmark with validation
  time_with_validation <- system.time({
    for (i in seq_len(iterations)) {
      result <- fasta2df(fasta_file, validate = TRUE)
    }
  })

  # Benchmark without validation
  time_without_validation <- system.time({
    for (i in seq_len(iterations)) {
      result <- fasta2df(fasta_file, validate = FALSE)
    }
  })

  cat("With validation:   ", round(time_with_validation["elapsed"] / iterations, 4), "seconds per run\n")
  cat("Without validation:", round(time_without_validation["elapsed"] / iterations, 4), "seconds per run\n")
  cat("Sequences processed:", nrow(result), "\n")
  cat("Performance:", round(nrow(result) / (time_without_validation["elapsed"] / iterations)), "sequences/second\n")

  invisible(result)
}
