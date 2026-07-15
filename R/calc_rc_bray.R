#' Calculate Raup-Crick on Bray-Curtis (RCbray)
#'
#' Computes RCbray, the compositional null-model axis used (alongside
#' betaNTI/betaNRI) to infer community-assembly processes (Stegen et al. 2013).
#' RCbray ranges from -1 to 1; |RCbray| > 0.95 with |betaNTI| <= 2 indicates
#' dispersal-related processes. Requires NO phylogenetic tree.
#'
#' RCbray = (fraction of null Bray-Curtis values strictly greater than the
#' observed value, with the observed value counted in the denominator) - 0.5,
#' then doubled. Null communities are generated with
#' \code{picante::randomizeMatrix} (default null model \code{"independentswap"}).
#'
#' @param data Community matrix (rows = taxa, cols = samples) or a data.frame
#'   whose first column holds taxon IDs.
#' @param sample Optional sample metadata data.frame. When provided with
#'   \code{group_col}, RCbray is computed independently within each group.
#' @param group_col Column in \code{sample} holding group labels (default "group").
#' @param sample_id Column in \code{sample} matching \code{colnames(data)}
#'   (default "sample").
#' @param runs Number of null-model randomizations (default 1000).
#' @param null_model Null model for \code{picante::randomizeMatrix}
#'   (default "independentswap").
#' @param verbose Print progress (default TRUE).
#'
#' @return Long data.frame: \code{from, to, rc_bray} (plus \code{group} when
#'   grouped). Pairs use the lower-triangle convention to align with
#'   \code{calc_beta_nti()}.
#'
#' @references Stegen, J. C., et al. (2013). ISME J, 7(11):2069-2079.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' otu <- read.delim("asv_table.txt", row.names = 1)
#' rc <- calc_rc_bray(otu, runs = 1000)
#' }
calc_rc_bray <- function(data, sample = NULL, group_col = "group",
                         sample_id = "sample", runs = 1000,
                         null_model = "independentswap", verbose = TRUE) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a community matrix (taxa x samples) or data.frame")
  }

  data <- as.data.frame(data)
  if (ncol(data) > 0 && !is.numeric(data[[1]])) {
    feat_names <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    feat_names <- rownames(data)
    if (is.null(feat_names)) feat_names <- paste0("Taxa_", seq_len(nrow(data)))
  }
  comm <- t(as.matrix(data))
  colnames(comm) <- feat_names
  rownames(comm) <- colnames(data)

  run_one <- function(comm_sub) {
    n <- nrow(comm_sub)
    obs <- as.matrix(vegan::vegdist(comm_sub, method = "bray"))
    if (verbose) message(sprintf("Randomizing RCbray (%d runs, %d pairs)...",
                                 runs, n * (n - 1) / 2))
    greater <- matrix(0L, nrow = n, ncol = n, dimnames = dimnames(obs))
    for (r in seq_len(runs)) {
      null_d <- as.matrix(vegan::vegdist(
        picante::randomizeMatrix(comm_sub, null.model = null_model), "bray"))
      greater <- greater + (null_d > obs)
      if (verbose && r %% 100 == 0) message(sprintf("  Run %d / %d", r, runs))
    }
    rc <- (greater / (runs + 1) - 0.5) * 2
    rc[upper.tri(rc, diag = TRUE)] <- NA   # keep lower triangle (align w/ calc_beta_nti)
    data.frame(from = rownames(rc), rc, stringsAsFactors = FALSE) %>%
      tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "rc_bray") %>%
      dplyr::filter(!is.na(.data$rc_bray))
  }

  if (is.null(sample)) {
    return(run_one(comm))
  }

  sample <- as.data.frame(sample)
  if (!is.null(sample_id) && sample_id %in% names(sample)) {
    rownames(sample) <- as.character(sample[[sample_id]])
  } else {
    rn <- rownames(sample)
    if (is.null(rn) || identical(rn, as.character(seq_len(nrow(sample))))) {
      for (col in names(sample)) {
        if (is.character(sample[[col]]) && all(rownames(comm) %in% sample[[col]])) {
          rownames(sample) <- sample[[col]]
          break
        }
      }
    }
  }
  if (!group_col %in% names(sample)) {
    stop(sprintf("Column '%s' not found in sample data", group_col))
  }
  groups <- unique(as.character(sample[[group_col]]))
  sample_groups <- stats::setNames(as.character(sample[[group_col]]), rownames(sample))

  results <- list()
  for (grp in groups) {
    grp_samples <- intersect(names(sample_groups[sample_groups == grp]), rownames(comm))
    if (length(grp_samples) < 2) {
      if (verbose) message(sprintf("Group '%s': skipped (< 2 samples)", grp))
      next
    }
    if (verbose) message(sprintf("\n=== Group: %s (%d samples) ===", grp, length(grp_samples)))
    res <- run_one(comm[grp_samples, , drop = FALSE])
    res$group <- grp
    results[[grp]] <- res
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
