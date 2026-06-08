#' Calculate Beta Nearest Taxon Index (betaNTI)
#'
#' Computes betaNTI for phylogenetic beta diversity analysis. betaNTI measures
#' whether observed phylogenetic turnover between sample pairs is significantly
#' different from random expectation. Values > +2 indicate deterministic
#' (environmental) processes; values < -2 indicate homogeneous selection.
#'
#' @param data Community data matrix (rows = taxa, cols = samples) or a
#'   data.frame where the first column contains taxon IDs.
#' @param tree A phylogenetic tree object (\code{ape::phylo}).
#' @param beta.reps Number of randomizations (default: 999).
#' @param abundance.weighted Logical; weight by species abundance (default: TRUE).
#' @param verbose Print progress messages (default: TRUE).
#'
#' @return A data frame in long format with columns: \code{from}, \code{to},
#'   \code{beta_mntd}, \code{beta_nti}.
#' @export
#' @seealso \code{\link[picante]{comdistnt}}, \code{\link[picante]{match.phylo.data}}
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("tree.nwk")
#' otu <- read.delim("asv_table.txt", row.names = 1)
#'
#' res <- calc_beta_nti(otu, tree, beta.reps = 999)
#' head(res)
#' }
calc_beta_nti <- function(data, tree, beta.reps = 999,
                          abundance.weighted = TRUE, verbose = TRUE) {

  # --- Prepare data ---
  data <- as.data.frame(data)
  if (!is.numeric(data[[1]])) {
    feat_names <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    feat_names <- rownames(data)
    if (is.null(feat_names)) feat_names <- paste0("Taxa_", seq_len(nrow(data)))
  }
  data_mat <- as.matrix(data)
  rownames(data_mat) <- feat_names

  # Match phylogenetic data
  matched <- picante::match.phylo.data(tree, data_mat)
  otu_matched <- matched$data
  phylo_matched <- matched$phylo
  cophenetic_mat <- ape::cophenetic.phylo(phylo_matched)

  if (verbose) {
    message(sprintf("Matched: %d taxa, %d samples", nrow(otu_matched), ncol(otu_matched)))
    message(sprintf("Computing betaMNTD (observed)..."))
  }

  # --- Observed betaMNTD ---
  beta_mntd <- as.matrix(picante::comdistnt(
    t(otu_matched), cophenetic_mat,
    abundance.weighted = abundance.weighted
  ))

  n_samp <- ncol(otu_matched)
  if (verbose) message(sprintf("Randomizing (%d reps, %d sample pairs)...", beta.reps, n_samp * (n_samp - 1) / 2))

  # --- Null distribution via taxa shuffle ---
  rand_bmntd <- array(NA_real_, dim = c(n_samp, n_samp, beta.reps))

  for (rep in seq_len(beta.reps)) {
    rand_tree_coph <- picante::taxaShuffle(cophenetic_mat)
    rand_bmntd[, , rep] <- as.matrix(picante::comdistnt(
      t(otu_matched), rand_tree_coph,
      abundance.weighted = abundance.weighted,
      exclude.conspecifics = FALSE
    ))
    if (verbose && rep %% 100 == 0) {
      message(sprintf("  Rep %d / %d (%s)", rep, beta.reps, date()))
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
  # Upper triangle of beta_mntd
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

  if (verbose) {
    n_det <- sum(abs(result$beta_nti) > 2, na.rm = TRUE)
    n_sto <- sum(abs(result$beta_nti) <= 2, na.rm = TRUE)
    message(sprintf("\nbetaNTI summary: %d pairs total", nrow(result)))
    message(sprintf("  |betaNTI| > 2 (deterministic): %d (%.1f%%)",
                    n_det, 100 * n_det / nrow(result)))
    message(sprintf("  |betaNTI| <= 2 (stochastic):    %d (%.1f%%)",
                    n_sto, 100 * n_sto / nrow(result)))
  }

  result
}
