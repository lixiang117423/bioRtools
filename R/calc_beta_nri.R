#' Calculate Beta Net Relatedness Index (betaNRI)
#'
#' betaNRI is the standardized effect size of betaMPD (mean pairwise distance),
#' the deeper-phylogenetic complement to betaNTI (which uses nearest-taxon
#' distance, betaMNTD). Used alongside RCbray for the full Stegen (2013)
#' community-assembly process framework. Requires a phylogenetic tree.
#'
#' Mirrors [calc_beta_nti()]; the only difference is the underlying phylogenetic
#' turnover metric: betaNRI uses \code{picante::comdist} (MPD, whole-tree signal),
#' betaNTI uses \code{picante::comdistnt} (MNTD, terminal-tip signal).
#'
#' @param data Community matrix (rows = taxa, cols = samples) or a data.frame
#'   whose first column holds taxon IDs.
#' @param tree A phylogenetic tree (\code{ape::phylo}).
#' @param beta_reps Number of null-model randomizations (default 999).
#' @param abundance_weighted Weight by abundance (default TRUE).
#' @param verbose Print progress (default TRUE).
#' @param sample Optional sample metadata for per-group computation.
#' @param sample_id Column in \code{sample} matching \code{colnames(data)}.
#' @param group_col Column in \code{sample} with group labels.
#'
#' @return Long data.frame: \code{from, to, beta_mpd, beta_nri} (plus
#'   \code{group} when grouped).
#'
#' @references Stegen, J. C., et al. (2013). ISME J, 7(11):2069-2079.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("tree.nwk")
#' otu <- read.delim("asv_table.txt", row.names = 1)
#' res <- calc_beta_nri(otu, tree, beta_reps = 999)
#' }
calc_beta_nri <- function(data, tree, beta_reps = 999,
                          abundance_weighted = TRUE, verbose = TRUE,
                          sample = NULL, sample_id = "sample",
                          group_col = "group") {
  if (is.null(tree)) stop("A phylogenetic tree is required")
  calc_ses_beta_core(
    data = data, tree = tree, beta_reps = beta_reps,
    abundance_weighted = abundance_weighted, verbose = verbose,
    sample = sample, sample_id = sample_id, group_col = group_col,
    dist_fun = "mpd", obs_name = "beta_mpd", ses_name = "beta_nri"
  )
}
