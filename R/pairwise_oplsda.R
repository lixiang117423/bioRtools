#' Pairwise Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' Performs OPLS-DA on all pairwise group combinations. For each pair of groups,
#' subsets the data, fits an OPLS-DA model via \code{ropls::opls}, and extracts
#' VIP scores, sample scores, and model performance metrics.
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

  # --- Prepare data ---------------------------------------------------------
  data <- as.data.frame(data)
  data <- orient_to_sample_row(data, sample, NULL, feature_as_row, FALSE)
  data_matrix <- as.matrix(data)
  sample <- as.data.frame(sample)

  # Auto-detect sample ID column
  rn <- rownames(sample)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(sample))))) {
    for (col in names(sample)) {
      if (is.character(sample[[col]]) && all(rownames(data_matrix) %in% sample[[col]])) {
        rownames(sample) <- sample[[col]]
        break
      }
    }
  }

  # Align samples
  common <- intersect(rownames(data_matrix), rownames(sample))
  if (length(common) == 0) {
    stop("No matching samples between data row names and sample row names")
  }
  data_matrix <- data_matrix[common, , drop = FALSE]
  sample      <- sample[common, , drop = FALSE]

  if (!group_col %in% names(sample)) {
    stop(sprintf("Column '%s' not found in sample data", group_col))
  }

  groups <- unique(as.character(sample[[group_col]]))
  if (length(groups) < 2) stop("At least 2 groups are required")

  variable_names <- colnames(data_matrix)
  if (is.null(variable_names)) {
    variable_names <- paste0("Var_", seq_len(ncol(data_matrix)))
  }

  # --- Pairwise OPLS-DA -----------------------------------------------------
  pairs <- utils::combn(groups, 2, simplify = FALSE)

  results <- lapply(pairs, function(pair) {
    idx <- sample[[group_col]] %in% pair
    data_sub  <- data_matrix[idx, , drop = FALSE]
    samp_sub  <- sample[idx, , drop = FALSE]
    group_sub <- factor(as.character(samp_sub[[group_col]]), levels = pair)
    comp_label <- paste(pair[1], "vs", pair[2])

    min_size <- min(table(group_sub))
    if (min_size < 3) {
      warning(sprintf("%s: group size < 3 (%d), skipping", comp_label, min_size))
      return(NULL)
    }

    # Remove zero-variance variables in this subset
    vars <- apply(data_sub, 2, var, na.rm = TRUE)
    keep <- vars > 0
    if (sum(keep) < 2) {
      warning(sprintf("%s: fewer than 2 non-constant variables, skipping", comp_label))
      return(NULL)
    }
    data_sub <- data_sub[, keep, drop = FALSE]
    var_names_sub <- variable_names[keep]

    # Fit OPLS-DA
    cv_val <- if (validation == "CV") min(cv_folds, min_size) else 0

    opls_model <- tryCatch(
      ropls::opls(
        x        = data_sub,
        y        = group_sub,
        predI    = 1,
        orthoI   = ortho_components,
        scaleC   = scaling,
        crossvalI = cv_val,
        fig.pdfC = "none",
        info.txtC = "none"
      ),
      error = function(e) {
        warning(sprintf("%s: OPLS fitting failed - %s", comp_label, e$message))
        NULL
      }
    )
    if (is.null(opls_model)) return(NULL)

    # VIP scores
    vip_values <- tryCatch(as.numeric(opls_model@vipVn), error = function(e) rep(NA_real_, ncol(data_sub)))
    vip_df <- data.frame(
      feature   = var_names_sub,
      vip       = vip_values,
      important = !is.na(vip_values) & vip_values >= vip_threshold,
      stringsAsFactors = FALSE
    )
    vip_df <- vip_df[order(vip_df$vip, decreasing = TRUE), ]
    vip_df$comparison <- comp_label

    # Scores
    scores_df <- data.frame(
      t1 = opls_model@scoreMN[, 1],
      sample_id = rownames(data_sub),
      group = group_sub,
      comparison = comp_label,
      stringsAsFactors = FALSE
    )
    if (!is.null(opls_model@orthoScoreMN) && ncol(opls_model@orthoScoreMN) > 0) {
      scores_df$to1 <- opls_model@orthoScoreMN[, 1]
    }

    # Model summary
    smry <- opls_model@summaryDF
    n_imp <- sum(vip_df$important, na.rm = TRUE)

    summary_row <- data.frame(
      comparison   = comp_label,
      R2Y          = if (!is.null(smry) && "R2Y" %in% names(smry)) smry$R2Y else NA_real_,
      Q2Y          = if (!is.null(smry) && "Q2Y" %in% names(smry)) smry$Q2Y else NA_real_,
      n_samples    = nrow(data_sub),
      n_variables  = ncol(data_sub),
      n_important  = n_imp,
      stringsAsFactors = FALSE
    )

    list(vip = vip_df, scores = scores_df, summary = summary_row, model = opls_model)
  })

  # --- Combine results ------------------------------------------------------
  # Remove NULL entries
  results <- results[!sapply(results, is.null)]

  if (length(results) == 0) {
    message("All pairwise OPLS-DA comparisons failed")
    return(list(
      vip_scores    = data.frame(),
      model_summary = data.frame(),
      scores        = data.frame(),
      models        = list()
    ))
  }

  vip_combined <- do.call(rbind, lapply(results, `[[`, "vip"))
  rownames(vip_combined) <- NULL
  vip_combined <- vip_combined[order(vip_combined$comparison, -vip_combined$vip), ]

  summary_combined <- do.call(rbind, lapply(results, `[[`, "summary"))
  rownames(summary_combined) <- NULL

  scores_combined <- do.call(rbind, lapply(results, `[[`, "scores"))
  rownames(scores_combined) <- NULL

  models_list <- setNames(lapply(results, `[[`, "model"),
                          sapply(results, function(x) x$summary$comparison))

  list(
    vip_scores    = vip_combined,
    model_summary = summary_combined,
    scores        = scores_combined,
    models        = models_list
  )
}
