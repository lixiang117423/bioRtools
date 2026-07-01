#' Pairwise Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' Performs OPLS-DA on all pairwise group combinations. For each pair of groups,
#' subsets the data, fits an OPLS-DA model via \code{ropls::opls}, and extracts
#' VIP scores, sample scores, and model performance metrics.
#'
#' \strong{Deprecated.} Use \code{\link{opls_analysis}} with
#' \code{pairwise = TRUE} instead. This function remains as a thin wrapper and
#' will be removed in the next major version.
#'
#' @param data Numerical matrix or data frame. Rows are samples and columns are
#'   features by default; the transpose (features as rows, samples as columns)
#'   is also accepted — see \code{feature_as_row}. Rows should match the samples
#'   in \code{sample}.
#' @param feature_as_row Logical or \code{NA}. \code{NA} (default) auto-detects
#'   the orientation by matching sample IDs from \code{sample} against the row
#'   and column names of \code{data}; \code{TRUE} forces features-as-rows;
#'   \code{FALSE} forces samples-as-rows. When detected or forced, the matrix is
#'   transposed internally so a manual \code{t()} is not needed.
#' @param sample Data frame of sample metadata. Row names should match
#'   \code{rownames(data)}. If row names are missing, the function will try to
#'   find a column whose values match \code{rownames(data)} and use it as row names.
#' @param group_col Column name in \code{sample} containing group labels
#'   (default: "group").
#' @param vip_threshold VIP score threshold for marking important variables
#'   (default: 1.0).
#' @param ortho_components Number of orthogonal components (default: 1).
#' @param scaling Scaling method passed to \code{ropls::opls}. One of
#'   "standard" (default), "pareto", "center", or "none".
#' @param validation Cross-validation method: "CV" (default) or "none".
#' @param cv_folds Number of cross-validation folds (default: 7).
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{vip_scores}}{Data frame of VIP scores for all comparisons,
#'       with columns \code{feature}, \code{vip}, \code{important}, and
#'       \code{comparison}. Sorted by comparison then VIP descending.}
#'     \item{\code{model_summary}}{Data frame with one row per comparison,
#'       containing R2Y, Q2Y, n_samples, n_variables, n_important, and the
#'       comparison label.}
#'     \item{\code{scores}}{Data frame of sample scores for all comparisons,
#'       with columns \code{t1}, \code{to1} (if available), \code{sample_id},
#'       \code{group}, and \code{comparison}.}
#'     \item{\code{models}}{Named list of raw \code{ropls} model objects.}
#'   }
#'
#' @seealso \code{\link{opls_analysis}} for single OPLS-DA analysis,
#'   \code{\link{find_dams_deseq2}} for pairwise DESeq2 differential analysis.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(100 * 18), nrow = 18, ncol = 100,
#'   dimnames = list(paste0("S", 1:18), paste0("Metab_", 1:100)))
#' meta <- data.frame(row.names = rownames(mat),
#'   treatment = rep(c("A", "B", "C"), each = 6))
#'
#' \dontrun{
#' res <- pairwise_oplsda(mat, meta, group_col = "treatment")
#' res$model_summary
#' head(res$vip_scores)
#' }
pairwise_oplsda <- function(data, sample, group_col = "group",
                            vip_threshold = 1.0, ortho_components = 1,
                            scaling = "standard", validation = "CV",
                            cv_folds = 7, feature_as_row = NA) {
  .Deprecated(
    "opls_analysis",
    package = "bioRtools",
    msg = paste("'pairwise_oplsda' is deprecated and will be removed in the",
                "next major version; use 'opls_analysis(..., pairwise = TRUE)'.")
  )

  res <- opls_analysis(
    data           = data,
    sample         = sample,
    group_col      = group_col,
    vip_threshold  = vip_threshold,
    ortho_components = ortho_components,
    scaling        = scaling,
    validation     = validation,
    cv_folds       = cv_folds,
    pairwise       = TRUE,
    feature_as_row = feature_as_row
  )

  # Reshape into the legacy shape.
  vip_scores <- res$vip_scores
  vip_scores <- vip_scores[order(vip_scores$comparison, -vip_scores$vip), ]
  rownames(vip_scores) <- NULL

  model_summary <- res$model_summary$pairwise_results
  model_summary <- do.call(rbind, lapply(model_summary, as.data.frame,
                                         stringsAsFactors = FALSE))
  rownames(model_summary) <- NULL

  scores <- res$scores
  rownames(scores) <- NULL

  list(
    vip_scores    = vip_scores,
    model_summary = model_summary,
    scores        = scores,
    models        = res$models
  )
}
