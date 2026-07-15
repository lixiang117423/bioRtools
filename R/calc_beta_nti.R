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
  result <- calc_ses_beta_core(
    data = data, tree = tree, beta_reps = beta_reps,
    abundance_weighted = abundance_weighted, verbose = verbose,
    sample = sample, sample_id = sample_id, group_col = group_col,
    dist_fun = "mntd", obs_name = "beta_mntd", ses_name = "beta_nti"
  )

  # 3-bucket coarse phylogenetic partition (legacy; for the full 5-process
  # partition use calc_assembly_process()).
  result <- dplyr::mutate(
    result,
    process = dplyr::case_when(
      .data$beta_nti > 2  ~ "Deterministic (variable selection)",
      .data$beta_nti < -2 ~ "Deterministic (homogeneous selection)",
      TRUE ~ "Stochastic"
    )
  )

  # Canonical column order (preserves pre-refactor layout exactly).
  result <- dplyr::select(
    result, from, to, beta_mntd, beta_nti, process,
    dplyr::any_of("group")
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
