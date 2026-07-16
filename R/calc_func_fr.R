#' Functional redundancy (FR)
#'
#' Computes per-sample, per-function functional redundancy from an ASV x
#' function assignment matrix (e.g. from [predict_func()]) and a count matrix.
#' Ported from microeco's `trans_func$cal_func_FR`.
#'
#' For sample s and function f:
#' \itemize{
#'   \item unweighted: \code{(N_taxa_with_f_in_s / N_taxa_in_s) * AF}
#'   \item weighted:   \code{(sum_abund_with_f_in_s / sum_abund_in_s) * AF}
#' }
#' where the optional adjustment factor \code{AF = N_unique_<level>_with_f /
#' N_unique_<level>_total} corrects for taxonomic dispersion (default AF = 1).
#'
#' @param asv_func ASV x function 0/1 matrix (from [predict_func()]).
#' @param abundance ASV x sample count matrix (same row names as `asv_func`).
#' @param abundance_weighted Use abundance-weighted FR (default FALSE).
#' @param adj_tax Apply the taxonomic adjustment factor (default FALSE).
#' @param taxonomy Taxonomy table (required if `adj_tax = TRUE`); ASV IDs in row
#'   names, rank columns per `tax_cols`.
#' @param adj_tax_by Rank for the adjustment factor (default "Genus").
#' @param tax_cols Named rank -> column mapping (default SILVA-style).
#' @param remove_zero Exclude zero-abundance taxa from each sample's totals
#'   (default TRUE).
#'
#' @return Long data.frame: `sample`, `function`, `fr` (plus `adjustment_factor`
#'   if `adj_tax = TRUE`).
#'
#' @references Louca et al. (2016); microeco `trans_func$cal_func_FR`.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' res <- predict_func(otu, taxo)
#' fr  <- calc_func_fr(res$asv_func, otu)
#' }
calc_func_fr <- function(asv_func, abundance, abundance_weighted = FALSE,
                         adj_tax = FALSE, taxonomy = NULL,
                         adj_tax_by = "Genus",
                         tax_cols = c(kingdom = "Kingdom", phylum = "Phylum",
                                      class = "Class", order = "Order",
                                      family = "Family", genus = "Genus",
                                      species = "Species"),
                         remove_zero = TRUE) {

  asv_func <- as.matrix(asv_func)
  abundance <- as.matrix(abundance)
  func_names <- colnames(asv_func)
  samples <- colnames(abundance)
  if (is.null(func_names) || is.null(samples)) {
    stop("asv_func and abundance need row/column names")
  }

  # --- optional adjustment factor (AF) per function ---
  AF <- stats::setNames(rep(1, length(func_names)), func_names)
  af_col <- NULL
  if (adj_tax) {
    if (is.null(taxonomy)) {
      stop("adj_tax = TRUE requires a taxonomy table")
    }
    col <- tax_cols[match(tolower(adj_tax_by), tolower(names(tax_cols)))]
    if (is.na(col)) stop("adj_tax_by rank not found in tax_cols names")
    taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE)
    level <- as.character(taxonomy[rownames(asv_func), col])
    uniq_total <- length(unique(stats::na.omit(level)))
    for (f in func_names) {
      uniq_with_f <- length(unique(stats::na.omit(level[asv_func[, f] == 1L])))
      AF[f] <- if (uniq_total > 0) uniq_with_f / uniq_total else NA_real_
    }
    af_col <- AF
  }

  rows <- vector("list", length(samples) * length(func_names))
  k <- 0L
  for (s in samples) {
    abund_s <- abundance[, s]
    present <- if (remove_zero) abund_s > 0 else rep(TRUE, length(abund_s))
    n_total <- sum(present)
    abund_total <- sum(abund_s[present])
    for (f in func_names) {
      has_f <- asv_func[, f] == 1L
      if (abundance_weighted) {
        num <- sum(abund_s[present & has_f])
        den <- abund_total
      } else {
        num <- sum(present & has_f)
        den <- n_total
      }
      fr <- if (isTRUE(den > 0)) (num / den) * AF[f] else NA_real_
      k <- k + 1L
      rows[[k]] <- data.frame(sample = s, func = f, fr = fr,
                              stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  names(out)[2] <- "function"
  if (adj_tax) {
    out$adjustment_factor <- af_col[as.character(out[["function"]])]
  }
  out
}


#' Community-level functional redundancy
#'
#' Geometric mean of per-function FR across functions, for each sample (zeros
#' floored at `floor_val` before logging). Ported from microeco's
#' `trans_func$cal_func_FR_comm`.
#'
#' @param func_fr Long data.frame from [calc_func_fr()] (`sample`, `function`, `fr`).
#' @param floor_val Small value replacing 0/NA before the log (default 1e-6).
#'
#' @return Data.frame: `sample`, `fr_comm`.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
calc_func_fr_comm <- function(func_fr, floor_val = 1e-6) {
  if (!all(c("sample", "function", "fr") %in% names(func_fr))) {
    stop("func_fr must have columns: sample, function, fr")
  }
  v <- pmax(func_fr$fr, floor_val, na.rm = TRUE)
  v[is.na(func_fr$fr)] <- NA_real_
  func_fr$.lfr <- log(v)
  agg <- stats::aggregate(.lfr ~ sample, data = func_fr, FUN = mean, na.rm = TRUE)
  agg$fr_comm <- exp(agg$.lfr)
  agg[, c("sample", "fr_comm")]
}
