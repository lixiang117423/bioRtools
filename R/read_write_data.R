#' Read Data File by Extension
#'
#' Reads a data file using the appropriate function based on file extension.
#' Supports Excel (.xlsx/.xls), CSV (.csv), TSV (.tsv), delimited text
#' (.txt), FASTA (.fasta/.fa), GFF3 (.gff/.gff3), and RDS/RData (.rds/.rdata)
#' files.
#' First row is treated as column headers for tabular formats.
#'
#' A forwarded \code{skip} argument (tabular formats) must be a single
#' non-negative integer line count; a non-integer such as \code{skip = "--"}
#' is rejected here with an actionable error instead of being passed to
#' vroom, which would fail with the opaque "C++ error (unknown cause)".
#'
#' A .xls/.xlsx file whose content is actually delimited text (common with
#' bioinformatics exporters — no ZIP/OLE2 magic bytes) is detected and read as
#' delimited text with the delimiter guessed from the first line, instead of
#' being handed to readxl where it would fail.
#'
#' For tabular formats (.xlsx/.xls/.csv/.tsv/.txt), the original column names
#' from the file are preserved as attribute \code{raw_names}. Use
#' \code{attr(df, "raw_names")} to retrieve them. Useful when \code{readxl} or
#' \code{readr} sanitizes names (e.g., \code{"Plant height (PH)"} becomes
#' \code{"Plant.height..PH.."}) — \code{raw_names} keeps the human-readable
#' originals for plotting, reporting, or pivoting back to original labels.
#'
#' For GFF3 (.gff/.gff3), directive/comment lines (those starting with
#' \code{#}) are dropped and the nine standard columns are returned with
#' fixed names (\code{seqid, source, type, start, end, score, strand, phase,
#' attributes}); column 9 is kept as the raw \code{key=value;...} string and
#' "\code{.}" is read as \code{NA}.
#'
#' @param file File path.
#' @param delim Character string used as field separator for \code{.txt} files.
#'   Default is \code{"\\t"} (tab). Ignored for other file types.
#' @param ... Additional arguments passed to the underlying read function.
#'
#' @return A data frame for tabular formats, the deserialized R object for
#'   \code{.rds}/\code{.rdata} files, a FASTA data frame for
#'   \code{.fasta}/\code{.fa} files, or a 9-column GFF3 data frame for
#'   \code{.gff}/\code{.gff3} files. For tabular formats, attribute
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
#' df <- read_data("data/annotation.gff3")  # 9-column GFF3 data frame
#' }
read_data <- function(file, delim = "\t", ...) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # readr/vroom and readxl both treat 'skip' as a non-negative integer line
  # count. A non-integer such as skip = "--" is forwarded unchanged to vroom's
  # C++ backend, which dies with the opaque "C++ error (unknown cause)" instead
  # of an actionable message — validate it here so the failure is self-evident.
  dots <- list(...)
  if (!is.null(dots$skip)) {
    s <- dots$skip
    s_num <- suppressWarnings(as.numeric(s))
    if (length(s) != 1L || is.na(s) || is.na(s_num) ||
        s_num < 0 || s_num != floor(s_num)) {
      stop("'skip' must be a single non-negative integer (number of lines ",
           "to skip from the top of the file), not ", deparse(s))
    }
  }

  ext <- tolower(tools::file_ext(file))

  switch(ext,
    xlsx = , xls = {
      if (is_excel_file(file)) {
        df <- readxl::read_excel(file, .name_repair = "universal", ...)
        raw_names <- names(suppressMessages(
          readxl::read_excel(file, n_max = 0, .name_repair = "minimal")
        ))
      } else {
        # Bioinformatics tools routinely emit tab/comma-delimited text with a
        # .xls extension; readxl rejects those (libxls "Unable to open file"),
        # so detect the real delimiter and read as plain text instead.
        dlm <- guess_delim(file)
        df <- readr::read_delim(file, delim = dlm,
                                name_repair = "universal", ...)
        raw_names <- names(readr::read_delim(file, delim = dlm, n_max = 0,
                                             show_col_types = FALSE,
                                             name_repair = "minimal"))
      }
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
    gff3 = , gff = gff3_to_df(file, ...),
    stop("Unsupported format: .", ext, "\n  Supported: .xlsx, .xls, .csv, .tsv, .txt, .fasta, .fa, .gff, .gff3, .rds, .rdata")
  )
}


#' Detect a genuine Excel workbook by magic bytes
#'
#' Real .xlsx files are ZIP containers (magic bytes \code{PK}); real .xls files
#' are OLE2 compound documents (magic bytes \code{D0 CF 11 E0}). Delimited text
#' exported with a .xls extension matches neither and must be read as text.
#'
#' @param file File path.
#' @return \code{TRUE} if the file is a genuine Excel binary, else \code{FALSE}.
#' @keywords internal
is_excel_file <- function(file) {
  magic <- readBin(file, "raw", n = 4)
  if (length(magic) < 4) {
    return(FALSE)
  }
  zip_magic <- as.raw(c(0x50, 0x4B))
  ole2_magic <- as.raw(c(0xD0, 0xCF, 0x11, 0xE0))
  identical(magic[1:2], zip_magic) || identical(magic, ole2_magic)
}

#' Guess a text table's field delimiter from its first line
#'
#' @param file File path.
#' @return \code{"\\t"} if the first line holds at least as many tabs as commas,
#'   otherwise \code{","}.
#' @keywords internal
guess_delim <- function(file) {
  first <- readLines(file, n = 1, warn = FALSE)
  n_tab <- length(regmatches(first, gregexpr("\t", first, fixed = TRUE))[[1L]])
  n_comma <- length(regmatches(first, gregexpr(",", first, fixed = TRUE))[[1L]])
  if (n_tab >= n_comma) "\t" else ","
}


#' Read a GFF3 file into a 9-column data frame
#'
#' GFF3 carries no column header and uses "#" for directive/comment lines
#' (\code{##gff-version}, \code{##sequence-region}, the \code{###} feature
#' separator, etc.). \code{readr}'s \code{comment = "#"} would also drop any
#' "\code{#}" appearing inside an attribute value (e.g. \code{Note=foo#bar}),
#' silently truncating data, so comment lines are filtered here and only data
#' rows are handed to \code{read_tsv}. Column 9 stays a verbatim
#' \code{key=value;...} string; the "\code{.}" sentinel becomes \code{NA}.
#'
#' @param file Path to a \code{.gff}/\code{.gff3} file.
#' @param ... Passed to \code{readr::read_tsv} (e.g. \code{n_max}, \code{skip}).
#' @return A 9-column data frame: \code{seqid, source, type, start, end, score,
#'   strand, phase, attributes}. \code{start}/\code{end}/\code{phase} are integer,
#'   \code{score} is double, \code{attributes} is the raw attribute string.
#' @keywords internal
gff3_to_df <- function(file, ...) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  lines <- lines[nzchar(trimws(lines))]

  readr::read_tsv(
    I(paste(lines, collapse = "\n")),
    col_names = c("seqid", "source", "type", "start", "end",
                  "score", "strand", "phase", "attributes"),
    col_types = readr::cols(
      seqid      = readr::col_character(),
      source     = readr::col_character(),
      type       = readr::col_character(),
      start      = readr::col_integer(),
      end        = readr::col_integer(),
      score      = readr::col_double(),
      strand     = readr::col_character(),
      phase      = readr::col_integer(),
      attributes = readr::col_character()
    ),
    na = ".",
    show_col_types = FALSE,
    ...
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
