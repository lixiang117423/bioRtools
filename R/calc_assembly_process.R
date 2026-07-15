#' Infer community-assembly processes (Stegen framework)
#'
#' Classifies each sample pair into one of five ecological assembly processes
#' from a phylogenetic null-model metric (betaNTI or betaNRI) and the
#' compositional RCbray, following Stegen et al. (2013):
#'
#' \itemize{
#'   \item Variable selection: phylo > +threshold
#'   \item Homogeneous selection: phylo < -threshold
#'   \item Dispersal limitation: RCbray > +threshold & |phylo| <= threshold
#'   \item Homogeneous dispersal: RCbray < -threshold & |phylo| <= threshold
#'   \item Drift (undominated): |RCbray| <= threshold & |phylo| <= threshold
#' }
#'
#' The phylo metric is auto-detected: if \code{beta_nti} has a \code{beta_nti}
#' column it is used; otherwise a \code{beta_nri} column is used (i.e. pass the
#' output of [calc_beta_nri()]). This refines the coarse 3-bucket \code{process}
#' column of [calc_beta_nti()] by splitting its "Stochastic" bucket using RCbray.
#'
#' @param beta_nti Output of [calc_beta_nti()] or [calc_beta_nri()]:
#'   long data.frame with \code{from, to} and a \code{beta_nti}/\code{beta_nri}
#'   column (plus optional \code{group}).
#' @param rc_bray Output of [calc_rc_bray()]: long data.frame with
#'   \code{from, to, rc_bray} (plus optional \code{group}).
#' @param thresholds Named numeric vector \code{c(phylo = 2, comm = 0.95)}.
#'   \code{phylo} is the |betaNTI/betaNRI| threshold; \code{comm} the |RCbray|
#'   threshold.
#' @param verbose Print a percentage summary (default TRUE).
#'
#' @return A list:
#'   \item{pair_process}{The joined long table with an added ordered
#'     \code{process} factor.}
#'   \item{summary}{Per-process percentage; with a \code{group} column when the
#'     inputs were grouped.}
#'
#' @references Stegen, J. C., et al. (2013). ISME J, 7(11):2069-2079.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' bnti <- calc_beta_nti(otu, tree, beta_reps = 999)
#' rc   <- calc_rc_bray(otu, runs = 1000)
#' res  <- calc_assembly_process(bnti, rc)
#' res$summary
#' }
calc_assembly_process <- function(beta_nti, rc_bray,
                                  thresholds = c(phylo = 2, comm = 0.95),
                                  verbose = TRUE) {

  if (!is.data.frame(beta_nti) || !is.data.frame(rc_bray)) {
    stop("'beta_nti' and 'rc_bray' must be data frames")
  }
  if (!all(c("from", "to") %in% names(beta_nti)) ||
      !all(c("from", "to") %in% names(rc_bray))) {
    stop("Both inputs need 'from' and 'to' columns")
  }

  phylo_col <- if ("beta_nti" %in% names(beta_nti)) {
    "beta_nti"
  } else if ("beta_nri" %in% names(beta_nti)) {
    "beta_nri"
  } else {
    stop("找不到 beta_nti 或 beta_nri 列;请传入 calc_beta_nti() 或 calc_beta_nri() 的输出")
  }

  has_group <- "group" %in% names(beta_nti)
  if (has_group != ("group" %in% names(rc_bray))) {
    stop("beta_nti 与 rc_bray 的 group 设置不一致;请用相同的 sample/group_col 计算")
  }

  pair_key <- function(df) {
    apply(df[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "|"))
  }

  b <- beta_nti %>%
    dplyr::mutate(.pair = pair_key(.)) %>%
    dplyr::select(from, to, .pair, dplyr::all_of(phylo_col), dplyr::any_of("group"))
  r <- rc_bray %>%
    dplyr::mutate(.pair = pair_key(.)) %>%
    dplyr::select(.pair, rc_bray, dplyr::any_of("group"))

  # Require the pair sets to match exactly — a mismatch means betaNTI/betaNRI
  # and RCbray were computed on different sample sets (user error).
  if (!setequal(b$.pair, r$.pair)) {
    stop("beta_nti 和 rc_bray 的样本对不一致;请确保两者来自同一组样本")
  }
  joined <- dplyr::inner_join(b, r, by = ".pair")
  if (has_group) {
    joined <- joined %>%
      dplyr::mutate(group = dplyr::coalesce(group.x, group.y)) %>%
      dplyr::select(-group.x, -group.y)
  }

  t_phylo <- thresholds[["phylo"]]
  t_comm <- thresholds[["comm"]]
  joined <- joined %>%
    dplyr::mutate(
      process = dplyr::case_when(
        .data[[phylo_col]] >  t_phylo ~ "Variable selection",
        .data[[phylo_col]] < -t_phylo ~ "Homogeneous selection",
        rc_bray >  t_comm             ~ "Dispersal limitation",
        rc_bray < -t_comm             ~ "Homogeneous dispersal",
        TRUE                          ~ "Drift"
      )
    )

  proc_levels <- c("Variable selection", "Homogeneous selection",
                   "Dispersal limitation", "Homogeneous dispersal", "Drift")
  joined$process <- factor(joined$process, levels = proc_levels)

  if (has_group) {
    summary_df <- joined %>%
      dplyr::group_by(group, process, .drop = FALSE) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(percentage = n / sum(n) * 100) %>%
      dplyr::select(group, process, percentage) %>%
      dplyr::arrange(group, process)
  } else {
    summary_df <- joined %>%
      dplyr::count(process, .drop = FALSE) %>%
      dplyr::mutate(percentage = n / sum(n) * 100) %>%
      dplyr::select(process, percentage) %>%
      dplyr::arrange(process)
  }

  pair_process <- joined %>%
    dplyr::select(from, to, dplyr::all_of(phylo_col), rc_bray, process,
                  dplyr::any_of("group"))

  if (verbose) {
    cat("Assembly-process percentages:\n")
    print(as.data.frame(summary_df))
  }

  list(pair_process = as.data.frame(pair_process),
       summary = as.data.frame(summary_df))
}
