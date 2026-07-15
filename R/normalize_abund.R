#' Normalize microbiome abundance data
#'
#' Unified dispatcher for common microbiome count normalization methods, ported
#' from microeco's `trans_norm`. Operates on a taxa x samples matrix (rows =
#' taxa, cols = samples) to align with [calc_beta_nti()] / [calc_rc_bray()].
#'
#' Methods:
#' \itemize{
#'   \item \code{clr}: centered log-ratio per sample (\code{log(x + pseudocount) - column mean}).
#'   \item \code{tss}: total-sum scaling (relative abundance).
#'   \item \code{hellinger}: square root of relative abundance.
#'   \item \code{css}: cumulative-sum scaling via \pkg{metagenomeSeq} (requires the package).
#'   \item \code{tmm}: trimmed mean of M-values via \pkg{edgeR}.
#'   \item \code{rle}: relative log expression via \pkg{edgeR}.
#'   \item \code{gmpr}: geometric mean of pairwise ratios (own implementation).
#'   \item \code{rarefy}: rarefaction to a fixed depth via \pkg{vegan}.
#' }
#'
#' For \code{css}, \code{tmm}, \code{rle}, and \code{gmpr} the per-sample
#' size/library factor is attached as attribute \code{"size_factor"}.
#'
#' @param data Community matrix (rows = taxa, cols = samples) or a data.frame
#'   whose first column holds taxa IDs.
#' @param method One of \code{"clr"}, \code{"tss"}, \code{"hellinger"},
#'   \code{"css"}, \code{"tmm"}, \code{"rle"}, \code{"gmpr"}, \code{"rarefy"}.
#' @param pseudocount Count added to zeros before \code{clr} (default 0.5).
#' @param depth Rarefaction depth for \code{rarefy}; \code{NULL} = min colSum.
#' @param seed Integer seed for \code{rarefy} reproducibility (default 123).
#' @param verbose Print progress (default TRUE).
#'
#' @return Normalized matrix (taxa x samples) with the same dimnames. For
#'   \code{css}/\code{tmm}/\code{rle}/\code{gmpr}, \code{attr(,"size_factor")}
#'   holds the per-sample factor.
#'
#' @references Chen et al. (2018) GMPR, PeerJ 6:e4600; Paulson et al. (2013) CSS,
#'   Nat Methods 10:1200; Robinson & Oshlack (2010) TMM, Genome Biol 11:R25;
#'   Aitchison (1986) clr.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' otu <- read.delim("asv_table.txt", row.names = 1)
#' clr <- normalize_abund(otu, method = "clr")
#' tss <- normalize_abund(otu, method = "tss")
#' }
normalize_abund <- function(data,
                            method = c("clr", "tss", "hellinger", "css",
                                       "tmm", "rle", "gmpr", "rarefy")[1],
                            pseudocount = 0.5,
                            depth = NULL,
                            seed = 123,
                            verbose = TRUE) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a community matrix (taxa x samples) or data.frame")
  }
  method <- match.arg(tolower(method),
                      c("clr", "tss", "hellinger", "css", "tmm", "rle", "gmpr", "rarefy"))

  # Prepare numeric matrix (taxa x samples)
  data <- as.data.frame(data)
  if (ncol(data) > 0 && !is.numeric(data[[1]])) {
    feat_names <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    feat_names <- rownames(data)
    if (is.null(feat_names)) feat_names <- paste0("Taxa_", seq_len(nrow(data)))
  }
  m <- as.matrix(data)
  rownames(m) <- feat_names

  if (method == "clr") {
    if (any(m == 0)) m <- m + pseudocount
    out <- log(m)
    out <- sweep(out, 2, colMeans(out), "-")
    return(out)
  }

  if (method == "tss") {
    cs <- colSums(m)
    if (any(cs == 0)) warning("Samples with zero total abundance produce NA columns")
    return(sweep(m, 2, cs, "/"))
  }

  if (method == "hellinger") {
    cs <- colSums(m)
    if (any(cs == 0)) warning("Samples with zero total abundance produce NA columns")
    return(sqrt(sweep(m, 2, cs, "/")))
  }

  if (method == "css") {
    if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
      stop("method = 'css' requires the 'metagenomeSeq' package; ",
           "install it with BiocManager::install('metagenomeSeq')")
    }
    obj <- metagenomeSeq::newMRexperiment(m)
    p <- tryCatch(metagenomeSeq::cumNormStatFast(obj), error = function(e) 0.5)
    obj <- metagenomeSeq::cumNorm(obj, p = p)
    out <- metagenomeSeq::MRcounts(obj, norm = TRUE)
    sf <- exp(metagenomeSeq::pData(obj)$normFactors)   # metagenomeSeq stores log factors
    attr(out, "size_factor") <- sf
    return(out)
  }

  if (method %in% c("tmm", "rle")) {
    nf <- edgeR::normLibSizes(m, method = toupper(method))
    effec <- colSums(m) * nf
    ref <- mean(effec)
    out <- sweep(m, 2, effec, "/") * ref
    attr(out, "size_factor") <- effec
    return(out)
  }

  if (method == "gmpr") {
    sf <- gmpr(m)
    out <- sweep(m, 2, sf, "/")
    attr(out, "size_factor") <- sf
    return(out)
  }

  # method == "rarefy"
  if (!isTRUE(all(m == floor(m)))) {
    stop("rarefy requires integer count data")
  }
  if (is.null(depth)) {
    depth <- min(colSums(m))
    if (verbose) message("rarefy depth (min colSum): ", depth)
  }
  set.seed(seed)
  r <- vegan::rrarefy(t(m), sample = depth)   # rrarefy: rows = samples
  t(r)
}
