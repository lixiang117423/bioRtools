#' Extract and Visualize Motif Information from MEME Software
#'
#' @description
#' Extract motif information from MEME software results in either TXT or XML format.
#' This function parses MEME output files and returns structured motif data suitable
#' for downstream analysis and visualization.
#'
#' @param data Character string. Path to the MEME output file (either .txt or .xml).
#' @param format Character string. Format of the MEME results file. Must be either
#'   "txt" or "xml". Default is "txt".
#'
#' @return A data frame (tibble) containing motif information:
#'   \itemize{
#'     \item For TXT format: motif number, sequence name, and motif sequence
#'     \item For XML format: detailed motif information including positions, p-values,
#'       and statistical measures
#'   }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename mutate filter select left_join
#' @importFrom stringr str_sub str_split str_trim
#' @importFrom XML xmlParse xmlRoot xmlSize xmlToList
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Load TXT format
#' filepath_txt <- system.file("examples", "meme.txt", package = "ggmotif")
#' motif_info <- get_motif_from_meme(data = filepath_txt, format = "txt")
#'
#' # Load XML format
#' filepath_xml <- system.file("examples", "meme.xml", package = "ggmotif")
#' motif_info <- get_motif_from_meme(data = filepath_xml, format = "xml")
#' }
#'
#' @export
get_motif_from_meme <- function(data, format = "txt") {
  # Input validation
  if (!file.exists(data)) {
    stop("File not found: ", data, call. = FALSE)
  }
  
  format <- tolower(format)
  if (!format %in% c("txt", "xml")) {
    stop("Format must be either 'txt' or 'xml', not '", format, "'", call. = FALSE)
  }
  
  # Dispatch to appropriate parser
  if (format == "txt") {
    parse_meme_txt(data)
  } else {
    parse_meme_xml(data)
  }
}


#' Parse MEME TXT Format
#'
#' @param file_path Path to MEME TXT file
#' @return A tibble with motif information
#' @keywords internal
#' @noRd
parse_meme_txt <- function(file_path) {
  # Read raw data
  raw_data <- readLines(file_path, warn = FALSE) %>%
    tibble::tibble(raw = .)
  
  # Find motif boundaries
  motif_boundaries <- raw_data %>%
    dplyr::mutate(
      row_num = dplyr::row_number(),
      line_length = nchar(.data$raw),
      suffix_length = 16,
      end_position = .data$line_length,
      start_position = .data$line_length - .data$suffix_length,
      last_chars = stringr::str_sub(.data$raw, .data$start_position, .data$end_position),
      first_chars = stringr::str_sub(.data$raw, 1, 2)
    ) %>%
    dplyr::filter(
      .data$last_chars == " in BLOCKS format" | .data$first_chars == "//"
    ) %>%
    dplyr::mutate(
      row_num = ifelse(
        .data$raw == "//",
        .data$row_num - 1,
        .data$row_num + 3
      )
    ) %>%
    dplyr::select(.data$raw, .data$row_num)
  
  # Extract motifs
  n_motifs <- nrow(motif_boundaries) / 2
  
  motif_data <- purrr::map_dfr(
    seq_len(n_motifs),
    function(i) {
      start_idx <- (i - 1) * 2 + 1
      end_idx <- start_idx + 1
      
      start_row <- motif_boundaries$row_num[start_idx]
      end_row <- motif_boundaries$row_num[end_idx]
      
      extract_motif_block(raw_data$raw, start_row, end_row, motif_num = i)
    }
  )
  
  return(motif_data)
}


#' Extract a Single Motif Block
#'
#' @param raw_lines Character vector of raw lines
#' @param start_row Starting row number
#' @param end_row Ending row number
#' @param motif_num Motif number
#' @return A tibble with extracted motif information
#' @keywords internal
#' @noRd
extract_motif_block <- function(raw_lines, start_row, end_row, motif_num) {
  motif_lines <- raw_lines[start_row:end_row]
  
  tibble::tibble(raw = motif_lines) %>%
    dplyr::mutate(
      motif_num = paste0("Motif.", motif_num),
      split_line = stringr::str_split(.data$raw, "\\s+")
    ) %>%
    dplyr::mutate(
      input_seq_name = purrr::map_chr(.data$split_line, ~.x[1]),
      input_seq_motif = purrr::map_chr(
        .data$split_line,
        ~.x[length(.x) - 3]
      )
    ) %>%
    dplyr::select(
      .data$raw,
      .data$motif_num,
      .data$input_seq_name,
      .data$input_seq_motif
    )
}


#' Parse MEME XML Format
#'
#' @param file_path Path to MEME XML file
#' @return A tibble with comprehensive motif information
#' @keywords internal
#' @noRd
parse_meme_xml <- function(file_path) {
  # Parse XML
  xml_doc <- XML::xmlParse(file = file_path)
  xml_root <- XML::xmlRoot(xml_doc)
  
  # Extract sequence information
  seq_info <- extract_sequence_info(xml_root[[1]])
  
  # Extract motif information
  motif_info <- extract_motif_info(xml_root[[3]])
  
  # Extract gene-motif mapping
  gene_motif_info <- extract_gene_motif_info(xml_root[[4]])
  
  # Combine all information
  combine_motif_data(gene_motif_info, motif_info, seq_info)
}


#' Extract Sequence Information from XML
#'
#' @param seq_node XML node containing sequence information
#' @return A tibble with sequence information
#' @keywords internal
#' @noRd
extract_sequence_info <- function(seq_node) {
  n_sequences <- XML::xmlSize(seq_node) - 1
  
  purrr::map_dfr(
    2:n_sequences,
    function(i) {
      seq_node[[i]] %>%
        XML::xmlToList() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::as_tibble()
    }
  ) %>%
    dplyr::rename(
      seq_id = 1,
      input_seq_id = 2
    )
}


#' Extract Motif Information from XML
#'
#' @param motif_node XML node containing motif information
#' @return A tibble with motif information
#' @keywords internal
#' @noRd
extract_motif_info <- function(motif_node) {
  n_motifs <- XML::xmlSize(motif_node)
  
  purrr::map_dfr(
    seq_len(n_motifs),
    function(i) {
      motif_node[[i]] %>%
        XML::xmlToList() %>%
        `[[`(".attrs") %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::as_tibble()
    }
  ) %>%
    dplyr::rename(
      motif_id = 1,
      motif_name = 2
    )
}


#' Extract Gene-Motif Mapping from XML
#'
#' @param gene_node XML node containing gene-motif mapping
#' @return A tibble with gene-motif relationships
#' @keywords internal
#' @noRd
extract_gene_motif_info <- function(gene_node) {
  n_genes <- XML::xmlSize(gene_node)
  
  purrr::map_dfr(
    seq_len(n_genes),
    function(i) {
      gene_data <- XML::xmlToList(gene_node[[i]])
      
      # Skip if character (no motif found)
      if (is.character(gene_data)) {
        return(NULL)
      }
      
      n_elements <- length(gene_data) - 1
      
      purrr::map_dfr(
        seq_len(n_elements),
        function(j) {
          gene_data[[j]] %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            tibble::as_tibble()
        }
      ) %>%
        dplyr::mutate(
          seq_id = gene_data[[".attrs"]][[1]],
          p_value_seq = gene_data[[".attrs"]][[2]],
          num_site4seq = gene_data[[".attrs"]][[3]]
        )
    }
  ) %>%
    dplyr::select(
      .data$seq_id,
      .data$num_site4seq,
      .data$p_value_seq,
      .data$motif_id,
      .data$strand,
      .data$position,
      .data$pvalue
    )
}


#' Combine All Motif Data
#'
#' @param gene_motif_info Gene-motif mapping tibble
#' @param motif_info Motif information tibble
#' @param seq_info Sequence information tibble
#' @return A comprehensive tibble with all motif data
#' @keywords internal
#' @noRd
combine_motif_data <- function(gene_motif_info, motif_info, seq_info) {
  gene_motif_info %>%
    dplyr::left_join(motif_info, by = "motif_id") %>%
    dplyr::left_join(seq_info[, 1:3], by = "seq_id") %>%
    dplyr::select(-.data$seq_id) %>%
    dplyr::mutate(
      position = as.numeric(.data$position),
      width = as.numeric(.data$width),
      start_position = .data$position + 1,
      end_position = .data$position + .data$width
    ) %>%
    dplyr::select(
      .data$input_seq_id,
      .data$length,
      .data$num_site4seq,
      .data$p_value_seq,
      .data$motif_id,
      .data$position,
      .data$width,
      .data$start_position,
      .data$end_position,
      .data$pvalue,
      .data$ic,
      .data$re,
      .data$llr,
      .data$p_value,
      .data$e_value,
      .data$bayes_threshold
    )
}


# Global variables (for R CMD check)
utils::globalVariables(c(
  ".", ".data"
))