#' Parse the bundled FAPROTAX database
#'
#' Reads the FAPROTAX `.txt` database (shipped in `inst/extdata/`) into a list
#' of function-to-clade mappings, composite set-operations, and metadata.
#' Internal helper for [predict_func()].
#'
#' @param path Path to FAPROTAX.txt (default: the bundled copy).
#'
#' @return A list: `func_tax` (named list, function -> character vector of
#'   `*`-delimited clade expressions), `func_ops` (named list, function ->
#'   list of `c(op, ref)` with op in add/subtract/intersect), `func_meta`
#'   (named list of metadata strings), `ver` (database version).
#'
#' @keywords internal
parse_faprotax <- function(path = system.file("extdata", "FAPROTAX.txt",
                                              package = "bioRtools")) {
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "extdata", "FAPROTAX.txt")   # in-package dev fallback
    if (!file.exists(path)) {
      stop("FAPROTAX.txt not found in inst/extdata; reinstall bioRtools")
    }
  }
  lines <- readLines(path, warn = FALSE)
  func_tax <- list(); func_ops <- list(); func_meta <- list()
  cur <- NULL; ver <- NA_character_
  op_re <- "^(add_group|subtract_group|intersect_group):"

  for (ln in lines) {
    trimmed <- trimws(ln)
    if (!nzchar(trimmed)) next
    if (startsWith(trimmed, "#")) {
      if (startsWith(trimmed, "# Version:")) {
        ver <- trimws(sub("^#\\s*Version:", "", trimmed))
      }
      next
    }
    ln_clean <- trimws(sub("#.*$", "", ln))   # strip trailing DOI/source comment
    if (!nzchar(ln_clean)) next

    if (grepl(op_re, ln_clean)) {
      op <- if (startsWith(ln_clean, "add_group:")) {
        "add"
      } else if (startsWith(ln_clean, "subtract_group:")) {
        "subtract"
      } else {
        "intersect"
      }
      ref <- trimws(sub(op_re, "", ln_clean))
      func_ops[[cur]] <- c(func_ops[[cur]], list(c(op = op, ref = ref)))
      next
    }

    if (startsWith(ln_clean, "*")) {
      func_tax[[cur]] <- c(func_tax[[cur]], ln_clean)
      next
    }

    # function header: "<name>[<ws><metadata...>]"
    parts <- strsplit(ln_clean, "\\s+", perl = TRUE)[[1]]
    cur <- parts[1]
    if (is.null(func_tax[[cur]])) func_tax[[cur]] <- character(0)
    if (is.null(func_ops[[cur]])) func_ops[[cur]] <- list()
    func_meta[[cur]] <- paste(parts[-1], collapse = " ")
  }

  list(func_tax = func_tax, func_ops = func_ops, func_meta = func_meta, ver = ver)
}
