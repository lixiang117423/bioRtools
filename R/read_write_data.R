#' Read Data File by Extension
#'
#' Reads a data file using the appropriate function based on file extension.
#' Supports Excel (.xlsx/.xls), CSV (.csv), TSV (.tsv), delimited text
#' (.txt), FASTA (.fasta/.fa), and RDS/RData (.rds/.rdata) files.
#' First row is treated as column headers for tabular formats.
#'
#' For tabular formats (.xlsx/.xls/.csv/.tsv/.txt), the original column names
#' from the file are preserved as attribute \code{raw_names}. Use
#' \code{attr(df, "raw_names")} to retrieve them. Useful when \code{readxl} or
#' \code{readr} sanitizes names (e.g., \code{"Plant height (PH)"} becomes
#' \code{"Plant.height..PH.."}) — \code{raw_names} keeps the human-readable
#' originals for plotting, reporting, or pivoting back to original labels.
#'
#' @param file File path.
#' @param delim Character string used as field separator for \code{.txt} files.
#'   Default is \code{"\\t"} (tab). Ignored for other file types.
#' @param ... Additional arguments passed to the underlying read function.
#'
#' @return A data frame for tabular formats, the deserialized R object for
#'   \code{.rds}/\code{.rdata} files, or a FASTA data frame for
#'   \code{.fasta}/\code{.fa} files. For tabular formats, attribute
#'   \code{raw_names} holds the original column names from the file.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{write_data}}
#'
#' @examples
#' \dontrun{
#' df <- read_data("data/samples.xlsx")
#' attr(df, "raw_names")  # original column names from the file
#'
#' df <- read_data("data/counts.csv")
#' df <- read_data("data/table.tsv")
#' df <- read_data("data/table.txt", delim = ",")
#' }
read_data <- function(file, delim = "\t", ...) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))

  switch(ext,
    xlsx = , xls = {
      df <- readxl::read_excel(file, .name_repair = "universal", ...)
      raw_names <- names(suppressMessages(
        readxl::read_excel(file, n_max = 0, .name_repair = "minimal")
      ))
      attr(df, "raw_names") <- raw_names
      df
    },
    csv = {
      df <- readr::read_csv(file, ...)
      raw_names <- names(readr::read_csv(file, n_max = 0,
                                          show_col_types = FALSE,
                                          name_repair = "minimal"))
      attr(df, "raw_names") <- raw_names
      df
    },
    tsv = {
      df <- readr::read_tsv(file, ...)
      raw_names <- names(readr::read_tsv(file, n_max = 0,
                                          show_col_types = FALSE,
                                          name_repair = "minimal"))
      attr(df, "raw_names") <- raw_names
      df
    },
    txt = {
      df <- readr::read_delim(file, delim = delim, ...)
      raw_names <- names(readr::read_delim(file, delim = delim, n_max = 0,
                                            show_col_types = FALSE,
                                            name_repair = "minimal"))
      attr(df, "raw_names") <- raw_names
      df
    },
    rds  = readRDS(file, ...),
    rdata = , rda = get(load(file, ...)),
    fasta = , fa = fasta2df(file, ...),
    stop("Unsupported format: .", ext, "\n  Supported: .xlsx, .xls, .csv, .tsv, .txt, .fasta, .fa, .rds, .rdata")
  )
}


#' Write Data File by Extension
#'
#' Writes a data frame to file using the appropriate function based on file
#' extension. Supports Excel (.xlsx), CSV (.csv), TSV (.tsv), tab-delimited
#' text (.txt), RDS (.rds), shell script (.sh), and image formats
#' (.pdf/.png/.svg/.tiff/.jpg/.eps) for saving ggplot objects.
#'
#' For .bed files, writes tab-delimited output with no header, following the
#' BED (Browser Extensible Data) format convention.
#'
#' For .sh files, writes the first column (or the \code{col} column) of the
#' data frame as plain text with no header and no quoting, suitable for
#' shell scripts or sample lists.
#'
#' For image formats, calls \code{ggplot2::ggsave} with sensible defaults
#' (width=8, height=6, dpi=600).
#'
#' @param data A data frame (or ggplot object for image formats) to write.
#' @param file File path. The extension determines the output format.
#' @param ... Additional arguments passed to the underlying write function.
#'
#' @return The input \code{data}, invisibly.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{read_data}}
#'
#' @examples
#' \dontrun{
#' write_data(df, "output/table.xlsx")
#' write_data(df, "output/table.csv")
#' write_data(df, "output/table.tsv")
#' }
write_data <- function(data, file, col = 1, col_names = TRUE,
                       width = 8, height = 6, dpi = 600, ...) {
  ext <- tolower(tools::file_ext(file))

  switch(ext,
    xlsx = writexl::write_xlsx(data, path = file, ...),
    csv  = readr::write_csv(data, file = file, col_names = col_names, ...),
    tsv  = readr::write_tsv(data, file = file, col_names = col_names, ...),
    bed  = readr::write_tsv(data, file = file, col_names = FALSE, ...),
    txt  = readr::write_delim(data, file = file, delim = "\t", col_names = col_names, ...),
    rds  = saveRDS(data, file = file, ...),
    sh   = readr::write_delim(data[, col, drop = FALSE], file = file,
              col_names = FALSE, quote = "none", delim = "\t", ...),
    fa = , fasta = if ("fasta_line" %in% names(data)) {
              writeLines(data$fasta_line, file)
            } else {
              df2fasta(data, output_file = file, ...)
            },
    pdf = ggplot2::ggsave(filename = file, plot = data, width = width, height = height,
              dpi = dpi, device = cairo_pdf, ...),
    png = , svg = , tiff = , jpg = , jpeg = , eps =
      ggplot2::ggsave(filename = file, plot = data, width = width, height = height, dpi = dpi, ...),
    stop("Unsupported format: .", ext, "\n  Supported: .xlsx, .csv, .tsv, .bed, .txt, .sh, .fa, .fasta, .pdf, .png, .svg, .tiff, .jpg, .eps, .rds")
  )
  invisible(data)
}


#' Write Data by Group
#'
#' Writes a grouped data frame to multiple files — one per group. Designed to
#' sit after \code{dplyr::group_by()} in a pipe: \code{file} is a path template
#' interpolated per group with \pkg{glue} syntax, where each \code{{column}}
#' placeholder is replaced by that group's value of the column.
#'
#' Each group is handed to \code{\link{write_data}} unchanged, so every format
#' \code{write_data} supports (.xlsx, .csv, .tsv, .bed, .txt, .rds, images, ...)
#' works here too. Multi-column groups are supported — use one placeholder per
#' grouping column, e.g. \code{"{CHROM}_{TYPE}.csv"}.
#'
#' @param data A grouped data frame (from \code{dplyr::group_by()}).
#' @param file File path template. Use \code{{column}} to insert a group's value
#'   of a grouping column, e.g. \code{"~/Downloads/result_{CHROM}.xlsx"}.
#' @param ... Additional arguments forwarded to \code{\link{write_data}} for
#'   every group (e.g. \code{col_names = FALSE}, or image \code{width}/\code{height}).
#'
#' @return The input \code{data}, invisibly, so the pipe can continue.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{write_data}}
#'
#' @examples
#' \dontrun{
#' df %>%
#'   dplyr::group_by(CHROM) %>%
#'   write_data_by_group("~/Downloads/GWAS_result_new_{CHROM}.xlsx")
#'
#' # Multi-column group, csv output
#' df %>%
#'   dplyr::group_by(CHROM, TYPE) %>%
#'   write_data_by_group("~/Downloads/{CHROM}_{TYPE}.csv")
#' }
write_data_by_group <- function(data, file, ...) {
  if (!dplyr::is_grouped_df(data)) {
    stop("'data' must be a grouped data frame; call dplyr::group_by() first")
  }

  dots <- list(...)

  dplyr::group_walk(data, function(.x, .y) {
    out_file <- glue::glue_data(.y, file)
    do.call(write_data, c(list(data = .x, file = out_file), dots))
  })

  invisible(data)
}
