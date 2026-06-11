#' Read Data File by Extension
#'
#' Reads a data file using the appropriate function based on file extension.
#' Supports Excel (.xlsx/.xls), CSV (.csv), TSV (.tsv), delimited text
#' (.txt), FASTA (.fasta/.fa), and RDS/RData (.rds/.rdata) files.
#' First row is treated as column headers for tabular formats.
#'
#' @param file File path.
#' @param delim Character string used as field separator for \code{.txt} files.
#'   Default is \code{"\\t"} (tab). Ignored for other file types.
#' @param ... Additional arguments passed to the underlying read function.
#'
#' @return A data frame for tabular formats, the deserialized R object for
#'   \code{.rds}/\code{.rdata} files, or a FASTA data frame for
#'   \code{.fasta}/\code{.fa} files.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link{write_data}}
#'
#' @examples
#' \dontrun{
#' df <- read_data("data/samples.xlsx")
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
    xlsx = , xls = readxl::read_excel(file, .name_repair = "universal", ...),
    csv  = readr::read_csv(file, ...),
    tsv  = readr::read_tsv(file, ...),
    txt  = readr::read_delim(file, delim = delim, ...),
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
