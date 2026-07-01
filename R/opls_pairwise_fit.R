#' Fit binary OPLS-DA models for a list of group pairs
#'
#' Internal helper shared by the pairwise modes of [opls_analysis()]. For each
#' pair, subsets the data, fits a binary OPLS-DA via \code{ropls::opls}, and
#' collects the model, sample scores, VIP scores, and a summary row. Used by
#' both the all-pairs mode (\code{pairwise = TRUE}) and the ref-anchored mode
#' (\code{ref_group} set, >2 groups).
#'
#' @param data_matrix Numeric matrix, samples x variables.
#' @param group Factor or character vector of group membership
#'   (length = \code{nrow(data_matrix)}).
#' @param pairs List of length-2 character vectors. \code{pair[[1]]} is the
#'   reference for that comparison, \code{pair[[2]]} the treatment. Comparison
#'   label is \code{"<pair[2]> vs <pair[1]>"}.
#' @param variable_names Character vector of variable (column) names.
#' @param ortho_components Number of orthogonal components (default 1).
#' @param scaling Scaling method forwarded to \code{ropls::opls}.
#' @param validation \code{"CV"} or \code{"none"}.
#' @param cv_folds Number of CV folds (reduced to min group size if needed).
#' @param verbose Logical; print per-comparison progress.
#'
#' @return Named list: \code{models} (named list of ropls objects keyed by
#'   comparison), \code{scores} (data frame), \code{vip} (data frame with
#'   feature/vip/group/ref_group/comparison), \code{summaries} (named list).
#' @keywords internal
fit_pairwise_opls <- function(data_matrix, group, pairs, variable_names,
                              ortho_components = 1, scaling = "standard",
                              validation = "CV", cv_folds = 7, verbose = FALSE) {
  group <- as.factor(group)
  models <- list()
  scores_list <- list()
  vip_list <- list()
  summaries <- list()

  for (pair in pairs) {
    ref <- pair[1]
    trt <- pair[2]
    comp_label <- paste(trt, "vs", ref)

    idx <- which(group %in% c(ref, trt))
    sub_data <- data_matrix[idx, , drop = FALSE]
    sub_group <- droplevels(factor(as.character(group[idx]), levels = c(ref, trt)))

    min_size <- min(table(sub_group))
    if (min_size < 3) {
      warning(sprintf("%s: group size < 3 (%d), skipping", comp_label, min_size))
      next
    }

    if (verbose) {
      message("Fitting OPLS-DA: ", comp_label, " (n=", nrow(sub_data), ")")
    }

    cvI <- if (validation == "CV") min(cv_folds, min_size) else 0

    model <- tryCatch(
      ropls::opls(sub_data, sub_group, predI = 1, orthoI = ortho_components,
                  scaleC = scaling, crossvalI = cvI,
                  fig.pdfC = "none", info.txtC = "none"),
      error = function(e) {
        warning("OPLS-DA failed for ", comp_label, ": ", e$message)
        NULL
      }
    )
    if (is.null(model)) next

    models[[comp_label]] <- model

    score_df <- data.frame(
      t1 = model@scoreMN[, 1],
      sample_id = rownames(sub_data),
      group = as.character(sub_group),
      comparison = comp_label,
      stringsAsFactors = FALSE
    )
    if (!is.null(model@orthoScoreMN) && ncol(model@orthoScoreMN) > 0) {
      score_df$to1 <- model@orthoScoreMN[, 1]
    }
    scores_list[[comp_label]] <- score_df

    vip_vals <- tryCatch(as.numeric(model@vipVn),
                         error = function(e) rep(NA_real_, ncol(sub_data)))
    vip_list[[comp_label]] <- data.frame(
      feature = variable_names,
      vip = vip_vals,
      group = trt,
      ref_group = ref,
      comparison = comp_label,
      stringsAsFactors = FALSE
    )

    s <- model@summaryDF
    summaries[[comp_label]] <- list(
      comparison = comp_label,
      group = trt, ref_group = ref,
      R2Y = if (!is.null(s) && "R2Y(cum)" %in% names(s)) s$`R2Y(cum)` else NA_real_,
      Q2Y = if (!is.null(s) && "Q2(cum)" %in% names(s)) s$`Q2(cum)` else NA_real_,
      R2X = if (!is.null(s) && "R2X(cum)" %in% names(s)) s$`R2X(cum)` else NA_real_,
      n_samples = nrow(sub_data),
      n_variables = ncol(sub_data)
    )
  }

  scores_df <- if (length(scores_list)) do.call(rbind, scores_list) else data.frame()
  vip_df <- if (length(vip_list)) do.call(rbind, vip_list) else data.frame()
  rownames(scores_df) <- NULL
  rownames(vip_df) <- NULL

  list(models = models, scores = scores_df, vip = vip_df, summaries = summaries)
}
