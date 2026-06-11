#' Calculate Beta Nearest Taxon Index (betaNTI)
#'
#' Computes betaNTI for phylogenetic beta diversity analysis following the
#' method described in Stegen et al. (2013). betaNTI measures whether
#' observed phylogenetic turnover between sample pairs is significantly
#' different from random expectation via taxa-label null model randomization.
#' Values > +2 indicate deterministic (variable selection) processes;
#' values < -2 indicate deterministic (homogeneous selection) processes.
#'
#' When \code{sample} and \code{group_col} are provided, each sample pair is
#' annotated with group labels and a \code{comparison} column, enabling
#' downstream summary by group pair.
#'
#' @param data Community data matrix (rows = taxa, cols = samples) or a
#'   data.frame where the first column contains taxon IDs.
#' @param tree A phylogenetic tree object (\code{ape::phylo}).
#' @param beta_reps Number of null model randomizations (default: 999).
#' @param abundance_weighted Logical; weight by species abundance (default: TRUE).
#' @param verbose Print progress messages (default: TRUE).
#' @param sample Optional data frame of sample metadata. Row names should match
#'   \code{colnames(data)}. If row names are missing, the function will try to
#'   find a column whose values match \code{colnames(data)}.
#' @param group_col Column name in \code{sample} containing group labels
#'   (default: "group"). Only used when \code{sample} is provided.
#' @param ref_group Optional character string specifying a reference group.
#'   When provided, only pairs involving this group are retained. Default is
#'   NULL (retain all pairs).
#'
#' @return A data frame with columns: \code{from}, \code{to},
#'   \code{beta_mntd}, \code{beta_nti}, \code{process}.
#'   When \code{sample} is provided, additional columns: \code{group_from},
#'   \code{group_to}, \code{comparison}.
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
#' # Without group annotation
#' res <- calc_beta_nti(otu, tree, beta_reps = 999)
#'
#' # With group annotation and reference group
#' res <- calc_beta_nti(otu, tree, beta_reps = 999,
#'   sample = meta, group_col = "treatment", ref_group = "Control")
#' }
calc_beta_nti <- function(data, tree, beta_reps = 999,
                          abundance_weighted = TRUE, verbose = TRUE,
                          sample = NULL, group_col = "group",
                          ref_group = NULL) {

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

  # Match tree tips and data taxa manually
  common_taxa <- intersect(tree$tip.label, rownames(data_mat))
  if (length(common_taxa) == 0) {
    stop("No matching taxa between tree tips and data row names")
  }
  phylo_matched <- ape::drop.tip(tree, setdiff(tree$tip.label, common_taxa))
  otu_matched <- data_mat[phylo_matched$tip.label, , drop = FALSE]
  cophenetic_mat <- ape::cophenetic.phylo(phylo_matched)

  if (verbose) {
    message(sprintf("Matched: %d taxa, %d samples", nrow(otu_matched), ncol(otu_matched)))
    message("Computing betaMNTD (observed)...")
  }

  # --- Observed betaMNTD ---
  # comdistnt expects samples as columns (community data transposed)
  beta_mntd <- as.matrix(picante::comdistnt(
    t(otu_matched), cophenetic_mat,
    abundance.weighted = abundance_weighted
  ))

  n_samp <- ncol(otu_matched)
  if (verbose) {
    message(sprintf("Randomizing (%d reps, %d sample pairs)...",
      beta_reps, n_samp * (n_samp - 1) / 2))
  }

  # --- Null distribution via taxa shuffle ---
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

  # --- Calculate betaNTI ---
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

  # --- Format output ---
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
    dplyr::filter(!is.na(.data$beta_nti))

  # Classify processes
  result <- result %>%
    dplyr::mutate(
      process = dplyr::case_when(
        .data$beta_nti > 2 ~ "Deterministic (variable selection)",
        .data$beta_nti < -2 ~ "Deterministic (homogeneous selection)",
        TRUE ~ "Stochastic"
      )
    )

  # --- Annotate with group information ---
  if (!is.null(sample)) {
    sample <- as.data.frame(sample)

    # Auto-detect sample ID column
    rn <- rownames(sample)
    if (is.null(rn) || identical(rn, as.character(seq_len(nrow(sample))))) {
      for (col in names(sample)) {
        if (is.character(sample[[col]]) && all(colnames(otu_matched) %in% sample[[col]])) {
          rownames(sample) <- sample[[col]]
          break
        }
      }
    }

    if (!group_col %in% names(sample)) {
      stop(sprintf("Column '%s' not found in sample data", group_col))
    }

    sample_groups <- stats::setNames(
      as.character(sample[[group_col]]),
      rownames(sample)
    )

    result$group_from <- sample_groups[result$from]
    result$group_to <- sample_groups[result$to]

    # Comparison label: consistent ordering (alphabetical)
    result <- result %>%
      dplyr::mutate(
        comparison = ifelse(group_from <= group_to,
          paste(group_from, "vs", group_to),
          paste(group_to, "vs", group_from))
      )

    # Filter by ref_group
    if (!is.null(ref_group)) {
      groups <- unique(c(result$group_from, result$group_to))
      if (!ref_group %in% groups) {
        stop("'ref_group' must be one of: ", paste(sort(groups), collapse = ", "))
      }
      result <- result %>%
        dplyr::filter(.data$group_from == ref_group | .data$group_to == ref_group)
    }

    # Group-level summary
    if (verbose) {
      summary_df <- result %>%
        dplyr::group_by(.data$comparison) %>%
        dplyr::summarise(
          n_pairs = dplyr::n(),
          mean_beta_nti = round(mean(.data$beta_nti, na.rm = TRUE), 3),
          median_beta_nti = round(stats::median(.data$beta_nti, na.rm = TRUE), 3),
          n_deterministic = sum(abs(.data$beta_nti) > 2, na.rm = TRUE),
          pct_deterministic = round(100 * sum(abs(.data$beta_nti) > 2, na.rm = TRUE) / dplyr::n(), 1),
          .groups = "drop"
        )

      message("\nbetaNTI summary by comparison:")
      for (i in seq_len(nrow(summary_df))) {
        message(sprintf("  %s: %d pairs, mean=%.3f, median=%.3f, deterministic=%d (%.1f%%)",
          summary_df$comparison[i],
          summary_df$n_pairs[i],
          summary_df$mean_beta_nti[i],
          summary_df$median_beta_nti[i],
          summary_df$n_deterministic[i],
          summary_df$pct_deterministic[i]))
      }
    }
  } else {
    if (verbose) {
      n_det <- sum(abs(result$beta_nti) > 2, na.rm = TRUE)
      n_sto <- sum(abs(result$beta_nti) <= 2, na.rm = TRUE)
      message(sprintf("\nbetaNTI summary: %d pairs total", nrow(result)))
      message(sprintf("  |betaNTI| > 2 (deterministic): %d (%.1f%%)",
                      n_det, 100 * n_det / nrow(result)))
      message(sprintf("  |betaNTI| <= 2 (stochastic):    %d (%.1f%%)",
                      n_sto, 100 * n_sto / nrow(result)))
    }
  }

  result
}
