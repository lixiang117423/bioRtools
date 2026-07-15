#' GMPR size factors (geometric mean of pairwise ratios)
#'
#' Computes per-sample size factors via the GMPR method (Chen et al. 2018,
#' <doi:10.7717/peerj.4600>): for each sample, the geometric mean over peers of
#' the median count ratio on shared non-zero taxa. Ported from microeco's
#' internal `GMPR`. Internal helper for [normalize_abund()].
#'
#' @param comm Community matrix, taxa x samples (features in rows).
#' @param intersect_no Minimum number of shared non-zero taxa required between
#'   a sample pair for it to count toward the geometric mean (default 10).
#' @param ct_min Minimum count threshold; counts below this are treated as 0
#'   (default 1).
#'
#' @return Named numeric vector of per-sample size factors (names = sample
#'   names). Samples that share fewer than `intersect_no` taxa with all peers
#'   get `NA` (with a warning). `attr(,"nss")` gives the number of qualifying
#'   peers per sample.
#'
#' @keywords internal
gmpr <- function(comm, intersect_no = 10, ct_min = 1) {
  comm[comm < ct_min] <- 0
  results <- sapply(seq_len(ncol(comm)), function(i) {
    x <- comm[, i]
    pr <- x / comm                      # taxa x samples ratio matrix vs every sample
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    included_count <- colSums(!is.na(pr))
    pr_median <- suppressWarnings(apply(pr, 2, median, na.rm = TRUE))
    nss_i <- sum(included_count >= intersect_no)
    gmpr_i <- if (nss_i > 1) {
      exp(mean(log(pr_median[included_count >= intersect_no])))
    } else {
      NA_real_
    }
    c(gmpr = gmpr_i, nss = nss_i)
  })
  output <- results["gmpr", ]
  nss <- results["nss", ]
  if (any(is.na(output))) {
    warning("Some samples share fewer than ", intersect_no,
            " common taxa with all peers; their size factors are NA.")
  }
  names(output) <- names(nss) <- colnames(comm)
  attr(output, "nss") <- nss
  output
}
