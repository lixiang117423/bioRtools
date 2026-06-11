#' Calculate Beta Nearest Taxon Index (betaNTI)
#'
#' Computes betaNTI for phylogenetic beta diversity analysis following the
#' method described in Stegen et al. (2013). betaNTI measures whether
#' observed phylogenetic turnover between sample pairs is significantly
#' different from random expectation via taxa-label null model randomization.
#' Values > +2 indicate deterministic (variable selection) processes;
#' values < -2 indicate deterministic (homogeneous selection) processes.
#'
#' When \code{sample} and \code{group_col} are provided, samples are split
#' by group and βNTI is calculated independently within each group. The
#' phylogenetic distance matrix (cophenetic) is computed once and reused
#' across groups to save memory.
#'
#' @param data Community data matrix (rows = taxa, cols = samples) or a
#'   data.frame where the first column contains taxon IDs.
#' @param tree A phylogenetic tree object (\code{ape::phylo}).
#' @param beta_reps Number of null model randomizations (default: 999).
#' @param abundance_weighted Logical; weight by species abundance (default: TRUE).
#' @param verbose Print progress messages (default: TRUE).
#' @param sample Optional data frame of sample metadata. Must contain a column
#'   with sample IDs matching \code{colnames(data)}.
#' @param sample_id Column name in \code{sample} containing sample IDs that match
#'   \code{colnames(data)}. Default is "sample".
#' @param group_col Column name in \code{sample} containing group labels
#'   (default: "group"). When provided with \code{sample}, βNTI is calculated
#'   independently within each group.
#'
#' @return A data frame with columns: \code{from}, \code{to},
#'   \code{beta_mntd}, \code{beta_nti}, \code{process}.
#'   When \code{sample} is provided, an additional \code{group} column
#'   indicates which group each result belongs to.
#'
#' @references
#' Stegen, J. C., Lin, X., Fredrickson, J. K., et al. (2013). Quantifying
#' community assembly processes and identifying features that impose them.
#' \emph{ISME J}, 7(11), 2069-2079.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("tree.nwk")
#' otu <- read.delim("asv_table.txt", row.names = 1)
#' meta <- read.delim("sample_metadata.txt", row.names = 1)
#'
#' # Calculate for all samples
#' res <- calc_beta_nti(otu, tree, beta_reps = 999)
#'
#' # Calculate within each group separately
#' res <- calc_beta_nti(otu, tree, beta_reps = 999,
#'   sample = meta, group_col = "treatment")
#' }
calc_beta_nti <- function(data, tree, beta_reps = 999,
                          abundance_weighted = TRUE, verbose = TRUE,
                          sample = NULL, sample_id = "sample",
                          group_col = "group") {

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

  # Match tree tips and data taxa (once)
  common_taxa <- intersect(tree$tip.label, rownames(data_mat))
  if (length(common_taxa) == 0) {
    stop("No matching taxa between tree tips and data row names")
  }
  phylo_matched <- ape::drop.tip(tree, setdiff(tree$tip.label, common_taxa))
  otu_matched <- data_mat[phylo_matched$tip.label, , drop = FALSE]

  # Compute cophenetic matrix once (most memory-intensive step)
  if (verbose) message(sprintf("Matched: %d taxa, %d samples", nrow(otu_matched), ncol(otu_matched)))
  if (verbose) message("Computing cophenetic distance matrix...")
  cophenetic_full <- ape::cophenetic.phylo(phylo_matched)
  rm(phylo_matched)
  gc(verbose = FALSE)

  # --- Determine groups ---
  if (!is.null(sample)) {
    sample <- as.data.frame(sample)

    # Use specified or auto-detect sample ID column
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

      # Subset samples and filter zero-count taxa
      grp_data <- otu_matched[, grp_samples, drop = FALSE]
      keep <- rowSums(grp_data) > 0
      grp_data <- grp_data[keep, , drop = FALSE]
      grp_coph <- cophenetic_full[rownames(grp_data), rownames(grp_data)]

      res <- calc_beta_nti_core(grp_data, grp_coph, beta_reps, abundance_weighted, verbose)
      res$group <- grp
      results[[grp]] <- res

      rm(grp_data, grp_coph)
      gc(verbose = FALSE)
    }

    out <- do.call(rbind, results)
    rownames(out) <- NULL
    return(out)
  }

  # No grouping: run on all samples
  calc_beta_nti_core(otu_matched, cophenetic_full, beta_reps, abundance_weighted, verbose)
}

#' @keywords internal
calc_beta_nti_core <- function(otu_matched, cophenetic_mat, beta_reps,
                               abundance_weighted, verbose) {

  n_samp <- ncol(otu_matched)

  if (verbose) message("Computing betaMNTD (observed)...")

  # Observed betaMNTD
  beta_mntd <- as.matrix(picante::comdistnt(
    t(otu_matched), cophenetic_mat,
    abundance.weighted = abundance_weighted
  ))

  if (verbose) {
    message(sprintf("Randomizing (%d reps, %d sample pairs)...",
      beta_reps, n_samp * (n_samp - 1) / 2))
  }

  # Null distribution via taxa shuffle
  rand_bmntd <- array(NA_real_, dim = c(n_samp, n_samp, beta_reps))

  for (rep in seq_len(beta_reps)) {
    rand_tree_coph <- picante::taxaShuffle(cophenetic_mat)
    rand_bmntd[, , rep] <- as.matrix(picante::comdistnt(
      t(otu_matched), rand_tree_coph,
      abundance.weighted = abundance_weighted,
      exclude.conspecifics = FALSE
    ))
    if (verbose && rep %% 100 == 0) {
      message(sprintf("  Rep %d / %d (%s)", rep, beta_reps, date()))
    }
  }

  # Calculate betaNTI
  if (verbose) message("Computing betaNTI...")
  beta_nti <- matrix(NA_real_, nrow = n_samp, ncol = n_samp)
  rownames(beta_nti) <- colnames(otu_matched)
  colnames(beta_nti) <- colnames(otu_matched)

  for (col in seq_len(n_samp - 1)) {
    for (row in (col + 1):n_samp) {
      rand_vals <- rand_bmntd[row, col, ]
      mu <- mean(rand_vals, na.rm = TRUE)
      sigma <- sd(rand_vals, na.rm = TRUE)
      if (sigma > 0) {
        beta_nti[row, col] <- (beta_mntd[row, col] - mu) / sigma
      }
    }
  }

  # Format output
  mntd_long <- data.frame(
    from = rownames(beta_nti),
    beta_mntd,
    stringsAsFactors = FALSE
  ) %>%
    tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "beta_mntd")

  nti_long <- data.frame(
    from = rownames(beta_nti),
    beta_nti,
    stringsAsFactors = FALSE
  ) %>%
    tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "beta_nti")

  result <- dplyr::left_join(mntd_long, nti_long, by = c("from", "to")) %>%
    dplyr::filter(!is.na(.data$beta_nti)) %>%
    dplyr::mutate(
      process = dplyr::case_when(
        .data$beta_nti > 2 ~ "Deterministic (variable selection)",
        .data$beta_nti < -2 ~ "Deterministic (homogeneous selection)",
        TRUE ~ "Stochastic"
      )
    )

  if (verbose) {
    n_det <- sum(abs(result$beta_nti) > 2, na.rm = TRUE)
    n_sto <- sum(abs(result$beta_nti) <= 2, na.rm = TRUE)
    message(sprintf("betaNTI: %d pairs, deterministic=%d (%.1f%%), stochastic=%d (%.1f%%)",
      nrow(result), n_det, 100 * n_det / nrow(result),
      n_sto, 100 * n_sto / nrow(result)))
  }

  result
}
