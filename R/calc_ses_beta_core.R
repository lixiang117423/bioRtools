#' Shared standardized-effect-size engine for betaNTI/betaNRI
#'
#' Internal workhorse for [calc_beta_nti()] (dist_fun = "mntd") and
#' [calc_beta_nri()] (dist_fun = "mpd"). Handles data preparation, tree
#' matching, optional per-group computation, the taxa-label null-model
#' randomization, and long-format output. Not exported.
#'
#' @param data Community matrix (rows = taxa, cols = samples) or a data.frame
#'   whose first column holds taxon IDs.
#' @param tree A phylogenetic tree (`ape::phylo`).
#' @param dist_fun `"mntd"` (picante::comdistnt, nearest taxon) or
#'   `"mpd"` (picante::comdist, mean pairwise).
#' @param obs_name Name for the observed-metric output column.
#' @param ses_name Name for the standardized-effect-size output column.
#' @inheritParams calc_beta_nti
#'
#' @keywords internal
calc_ses_beta_core <- function(data, tree, beta_reps = 999,
                               abundance_weighted = TRUE, verbose = TRUE,
                               sample = NULL, sample_id = "sample",
                               group_col = "group",
                               dist_fun = c("mntd", "mpd"),
                               obs_name = "beta_mntd",
                               ses_name = "beta_nti") {
  dist_fun <- match.arg(dist_fun)
  dist_fn <- if (dist_fun == "mntd") picante::comdistnt else picante::comdist

  # --- Prepare data ---
  data <- as.data.frame(data)
  if (ncol(data) > 0 && !is.numeric(data[[1]])) {
    feat_names <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    feat_names <- rownames(data)
    if (is.null(feat_names)) feat_names <- paste0("Taxa_", seq_len(nrow(data)))
  }
  data_mat <- as.matrix(data)
  rownames(data_mat) <- feat_names

  common_taxa <- intersect(tree$tip.label, rownames(data_mat))
  if (length(common_taxa) == 0) {
    stop("No matching taxa between tree tips and data row names")
  }
  phylo_matched <- ape::drop.tip(tree, setdiff(tree$tip.label, common_taxa))
  otu_matched <- data_mat[phylo_matched$tip.label, , drop = FALSE]

  if (verbose) message(sprintf("Matched: %d taxa, %d samples", nrow(otu_matched), ncol(otu_matched)))
  if (verbose) message("Computing cophenetic distance matrix...")
  cophenetic_full <- ape::cophenetic.phylo(phylo_matched)
  rm(phylo_matched)
  gc(verbose = FALSE)

  # Compute ses for one (comm, cophenetic) pair
  run_one <- function(comm, coph) {
    n_samp <- ncol(comm)
    if (verbose) message(sprintf("Computing observed %s...", obs_name))
    beta_obs <- as.matrix(dist_fn(t(comm), coph, abundance.weighted = abundance_weighted))
    if (verbose) message(sprintf("Randomizing (%d reps, %d sample pairs)...",
                                 beta_reps, n_samp * (n_samp - 1) / 2))
    rand_arr <- array(NA_real_, dim = c(n_samp, n_samp, beta_reps))
    for (rep in seq_len(beta_reps)) {
      rand_coph <- picante::taxaShuffle(coph)
      rand_arr[, , rep] <- as.matrix(dist_fn(t(comm), rand_coph,
                                             abundance.weighted = abundance_weighted))
      if (verbose && rep %% 100 == 0) {
        message(sprintf("  Rep %d / %d (%s)", rep, beta_reps, date()))
      }
    }
    if (verbose) message(sprintf("Computing %s...", ses_name))
    ses_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_samp)
    rownames(ses_mat) <- colnames(comm)
    colnames(ses_mat) <- colnames(comm)
    for (col in seq_len(n_samp - 1)) {
      for (row in (col + 1):n_samp) {
        rand_vals <- rand_arr[row, col, ]
        mu <- mean(rand_vals, na.rm = TRUE)
        sigma <- sd(rand_vals, na.rm = TRUE)
        if (sigma > 0) {
          ses_mat[row, col] <- (beta_obs[row, col] - mu) / sigma
        }
      }
    }
    obs_long <- data.frame(from = rownames(ses_mat), beta_obs,
                           stringsAsFactors = FALSE) %>%
      tidyr::pivot_longer(cols = -from, names_to = "to", values_to = obs_name)
    ses_long <- data.frame(from = rownames(ses_mat), ses_mat,
                           stringsAsFactors = FALSE) %>%
      tidyr::pivot_longer(cols = -from, names_to = "to", values_to = ses_name)
    dplyr::left_join(obs_long, ses_long, by = c("from", "to")) %>%
      dplyr::filter(!is.na(.data[[ses_name]]))
  }

  # --- Determine groups ---
  if (is.null(sample)) {
    return(run_one(otu_matched, cophenetic_full))
  }

  sample <- as.data.frame(sample)
  if (!is.null(sample_id) && sample_id %in% names(sample)) {
    rownames(sample) <- as.character(sample[[sample_id]])
  } else {
    rn <- rownames(sample)
    if (is.null(rn) || identical(rn, as.character(seq_len(nrow(sample))))) {
      for (col in names(sample)) {
        if (is.character(sample[[col]]) && all(colnames(otu_matched) %in% sample[[col]])) {
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
    grp_samples <- names(sample_groups[sample_groups == grp])
    grp_samples <- intersect(grp_samples, colnames(otu_matched))
    if (length(grp_samples) < 2) {
      if (verbose) message(sprintf("Group '%s': skipped (< 2 samples)", grp))
      next
    }
    if (verbose) message(sprintf("\n=== Group: %s (%d samples) ===", grp, length(grp_samples)))
    grp_data <- otu_matched[, grp_samples, drop = FALSE]
    keep <- rowSums(grp_data) > 0
    grp_data <- grp_data[keep, , drop = FALSE]
    grp_coph <- cophenetic_full[rownames(grp_data), rownames(grp_data)]
    res <- run_one(grp_data, grp_coph)
    res$group <- grp
    results[[grp]] <- res
    rm(grp_data, grp_coph)
    gc(verbose = FALSE)
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
