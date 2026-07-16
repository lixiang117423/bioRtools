#' Match a FAPROTAX clade expression against a taxonomy lineage
#'
#' Port of FAPROTAX `collapse_table.py` `find_matches_to_words_expression` in
#' 'words' mode (`valid_word_symbols = "-"`): the `*`-separated words in
#' `expression` must each appear in `candidate` as a complete word, in order.
#' A "complete word" means flanked by characters that are neither alphanumeric
#' nor `-`. Internal helper.
#'
#' @keywords internal
match_words <- function(expression, candidate, vws = "-") {
  words <- strsplit(expression, "*", fixed = TRUE)[[1]]
  words <- words[nzchar(words)]
  if (length(words) == 0) return(FALSE)
  cand_chars <- strsplit(candidate, "", fixed = TRUE)[[1]]
  LC <- length(cand_chars)
  vws_chars <- if (nzchar(vws)) strsplit(vws, "", fixed = TRUE)[[1]] else character(0)
  is_internal <- function(idx) {
    if (idx < 1 || idx > LC) return(FALSE)
    ch <- cand_chars[idx]
    grepl("[[:alnum:]]", ch) || ch %in% vws_chars
  }
  next_start <- 1L
  for (w in words) {
    positions <- gregexpr(w, candidate, fixed = TRUE)[[1]]
    if (positions[1] == -1L) return(FALSE)
    positions <- positions[positions >= next_start]
    matched <- FALSE
    wl <- nchar(w)
    for (pos in positions) {
      end <- pos + wl                       # index just past the match
      before_ok <- (pos == 1L) || !is_internal(pos - 1L)
      after_ok  <- (end > LC)   || !is_internal(end)
      if (before_ok && after_ok) {
        next_start <- end; matched <- TRUE; break
      }
    }
    if (!matched) return(FALSE)
  }
  TRUE
}

#' Predict prokaryote functions via FAPROTAX
#'
#' Assigns each ASV (taxon) to FAPROTAX functional groups by matching its
#' taxonomic lineage against the FAPROTAX clade database (bundled, v1.2.12,
#' BSD-2 license). Reimplements the official `collapse_table.py` 'words' mode
#' clade matching plus composite `add_group`/`subtract_group`/`intersect_group`
#' set operations. Ported from microeco's `trans_func$cal_func`, but using the
#' original BSD-2 database (not microeco's GPL-3 conversion).
#'
#' @param data Community count matrix (rows = taxa/ASVs, cols = samples) or a
#'   data.frame whose first column holds ASV IDs.
#' @param taxonomy Taxonomy table (data.frame). ASV IDs in row names (or first
#'   column); columns named per `tax_cols`. One rank per column.
#' @param db Database; only `"FAPROTAX"` is supported.
#' @param tax_cols Named character vector mapping rank (kingdom, phylum, class,
#'   order, family, genus, species) -> column name in `taxonomy`. Ranks are
#'   joined high-to-low with `;` to form each ASV's lineage.
#' @param case_insensitive Match case-insensitively (default TRUE).
#' @param verbose Print a summary (default TRUE).
#'
#' @return A list:
#'   \item{asv_func}{0/1 matrix, ASV x function.}
#'   \item{func_table}{Abundance matrix, function x sample (sum of ASV counts
#'     per function).}
#'   \item{version}{FAPROTAX database version.}
#'
#' @references Louca, S., et al. (2016). Decoupling function and taxonomy in
#'   the global ocean microbiome. *Science* 353:1272-1277. FAPROTAX:
#'   https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' otu  <- read.delim("asv_table.txt", row.names = 1)
#' taxo <- read.delim("taxonomy.txt", row.names = 1)
#' res  <- predict_func(otu, taxo)
#' res$func_table
#' }
predict_func <- function(data, taxonomy, db = "FAPROTAX",
                         tax_cols = c(kingdom = "Kingdom", phylum = "Phylum",
                                      class = "Class", order = "Order",
                                      family = "Family", genus = "Genus",
                                      species = "Species"),
                         case_insensitive = TRUE, verbose = TRUE) {

  if (db != "FAPROTAX") {
    stop("Only db = 'FAPROTAX' is supported in this version")
  }
  parsed <- parse_faprotax()
  func_tax <- parsed$func_tax
  func_ops <- parsed$func_ops
  func_names <- names(func_tax)

  # --- prepare count matrix (taxa x samples) ---
  data <- as.data.frame(data)
  if (ncol(data) > 0 && !is.numeric(data[[1]])) {
    asv_ids <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    asv_ids <- rownames(data)
    if (is.null(asv_ids)) asv_ids <- paste0("ASV_", seq_len(nrow(data)))
  }
  m <- as.matrix(data)
  rownames(m) <- asv_ids
  n_asv <- length(asv_ids)

  # --- build lineage strings from taxonomy ---
  taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE)
  if (!is.null(rownames(taxonomy))) {
    tax_ids <- rownames(taxonomy)
  } else {
    tax_ids <- as.character(taxonomy[[1]])
    taxonomy <- taxonomy[, -1, drop = FALSE]
  }
  use_cols <- tax_cols[tax_cols %in% names(taxonomy)]
  if (length(use_cols) == 0) {
    stop("None of the tax_cols were found in the taxonomy table")
  }
  tax_idx <- match(asv_ids, tax_ids)
  if (all(is.na(tax_idx))) {
    stop("taxonomy row names must match data row names")
  }

  lineages <- character(n_asv)
  for (i in seq_len(n_asv)) {
    j <- tax_idx[i]
    if (is.na(j)) { lineages[i] <- ""; next }
    vals <- as.character(unlist(taxonomy[j, use_cols, drop = FALSE]))
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    vals <- sub("^[a-z]__", "", vals)          # strip SILVA-style "k__" etc.
    lineages[i] <- paste(vals, collapse = ";")
  }
  if (case_insensitive) lineages <- tolower(lineages)
  if (any(lineages == "")) {
    warning("Some ASVs have empty taxonomy lineages and will match no function")
  }

  # --- primitive assignment: clade matching ---
  asv_func <- matrix(0L, nrow = n_asv, ncol = length(func_names),
                     dimnames = list(asv_ids, func_names))
  primitive <- !vapply(func_names, function(f) length(func_ops[[f]]) > 0, logical(1))
  for (f in func_names) {
    clades <- func_tax[[f]]
    if (length(clades) == 0) next
    for (cl in clades) {
      cl_use <- if (case_insensitive) tolower(cl) else cl
      for (i in seq_len(n_asv)) {
        if (asv_func[i, f] == 0L && nzchar(lineages[i]) &&
            match_words(cl_use, lineages[i])) {
          asv_func[i, f] <- 1L
        }
      }
    }
  }

  # --- composite set operations (file order => referents already finalized) ---
  for (f in func_names) {
    ops <- func_ops[[f]]
    if (length(ops) == 0) next
    idx <- which(asv_func[, f] == 1L)
    for (op in ops) {
      ref_idx <- which(asv_func[, op[["ref"]]] == 1L)
      op_kind <- op[["op"]]
      if (op_kind == "add") {
        idx <- union(idx, ref_idx)
      } else if (op_kind == "subtract") {
        idx <- setdiff(idx, ref_idx)
      } else {
        idx <- intersect(idx, ref_idx)
      }
    }
    asv_func[, f] <- 0L
    asv_func[idx, f] <- 1L
  }

  # --- function x sample abundance ---
  func_table <- t(t(m) %*% asv_func)

  if (verbose) {
    n_funcs_hit <- sum(colSums(asv_func) > 0)
    message(sprintf("FAPROTAX %s: %d ASVs, %d/%d functions assigned (%d primitive).",
                    parsed$ver, n_asv, n_funcs_hit, length(func_names), sum(primitive)))
  }

  list(asv_func = asv_func, func_table = func_table, version = parsed$ver)
}
