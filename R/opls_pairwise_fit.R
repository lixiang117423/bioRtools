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

#' Compute log2FC + p-value differential stats for a list of group pairs
#'
#' Internal helper shared by the pairwise modes of [opls_analysis()]. For each
#' pair, computes per-variable log2FC (treatment vs reference) and a p-value via
#' t-test or Wilcoxon (auto-chosen by Shapiro normality test when
#' \code{test_method = "auto"}), adjusts p-values per comparison, joins VIP
#' scores, and labels regulation (Up/Down/NS).
#'
#' @param data_matrix Numeric matrix, samples x variables.
#' @param group Factor or character vector of group membership.
#' @param pairs List of length-2 character vectors (\code{pair[1]=ref},
#'   \code{pair[2]=treatment}).
#' @param variable_names Character vector of column names.
#' @param vip Data frame with columns \code{feature}, \code{comparison},
#'   \code{vip} (from [fit_pairwise_opls()]).
#' @param test_method \code{"auto"}, \code{"t-test"}, or \code{"wilcoxon"}.
#' @param p_adjust_method P-value adjustment method (see \code{stats::p.adjust}).
#' @param p_threshold Regulation significance threshold.
#'
#' @return Data frame: \code{feature, group, ref_group, group_mean, group_sd,
#'   group_n, ref_mean, ref_sd, ref_n, log2_fc, p_value, test_method,
#'   comparison, p_adjust, vip, regulation}.
#' @keywords internal
compute_pairwise_diff <- function(data_matrix, group, pairs, variable_names,
                                  vip, test_method = "auto",
                                  p_adjust_method = "BH", p_threshold = 0.05) {
  group <- as.factor(group)
  diff_list <- list()

  for (pair in pairs) {
    ref <- pair[1]
    trt <- pair[2]
    ref_idx <- which(group == ref)
    trt_idx <- which(group == trt)
    ref_data <- data_matrix[ref_idx, , drop = FALSE]
    trt_data <- data_matrix[trt_idx, , drop = FALSE]

    n_ref <- nrow(ref_data)
    n_trt <- nrow(trt_data)
    nv <- length(variable_names)

    log2fc <- numeric(nv)
    pvals <- numeric(nv)
    tmethods <- character(nv)
    ref_mean <- numeric(nv)
    trt_mean <- numeric(nv)
    ref_sd <- numeric(nv)
    trt_sd <- numeric(nv)

    for (j in seq_len(nv)) {
      rv <- ref_data[, j]
      tv <- trt_data[, j]
      ref_mean[j] <- mean(rv, na.rm = TRUE)
      trt_mean[j] <- mean(tv, na.rm = TRUE)
      ref_sd[j] <- sd(rv, na.rm = TRUE)
      trt_sd[j] <- sd(tv, na.rm = TRUE)
      log2fc[j] <- log2((trt_mean[j] + 1e-10) / (ref_mean[j] + 1e-10))

      vrv <- rv[!is.na(rv)]
      vtv <- tv[!is.na(tv)]
      if (length(vrv) >= 2 && length(vtv) >= 2 &&
          stats::var(vrv) > 0 && stats::var(vtv) > 0) {
        use_t <- if (test_method == "auto") {
          rn <- length(vrv) >= 3 &&
            tryCatch(stats::shapiro.test(vrv)$p.value > 0.05, error = function(e) FALSE)
          tn <- length(vtv) >= 3 &&
            tryCatch(stats::shapiro.test(vtv)$p.value > 0.05, error = function(e) FALSE)
          rn && tn
        } else {
          test_method == "t-test"
        }
        if (use_t) {
          pvals[j] <- stats::t.test(vtv, vrv)$p.value
          tmethods[j] <- "t-test"
        } else {
          pvals[j] <- stats::wilcox.test(vtv, vrv)$p.value
          tmethods[j] <- "wilcoxon"
        }
      } else {
        pvals[j] <- NA
        tmethods[j] <- NA_character_
      }
    }

    diff_list[[paste(trt, "vs", ref)]] <- data.frame(
      feature = variable_names,
      group = trt,
      ref_group = ref,
      group_mean = round(trt_mean, 4),
      group_sd = round(trt_sd, 4),
      group_n = n_trt,
      ref_mean = round(ref_mean, 4),
      ref_sd = round(ref_sd, 4),
      ref_n = n_ref,
      log2_fc = round(log2fc, 4),
      p_value = pvals,
      test_method = tmethods,
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, diff_list)
  rownames(res) <- NULL
  res$comparison <- paste0(res$group, " vs ", res$ref_group)

  res <- res %>%
    dplyr::group_by(comparison) %>%
    dplyr::mutate(p_adjust = stats::p.adjust(p_value, method = p_adjust_method)) %>%
    dplyr::ungroup()

  # VIP is unique per (feature, comparison); join on both to stay unambiguous
  # in all-pairs mode where a treatment group appears in multiple comparisons.
  vip_lookup <- vip[, c("feature", "comparison", "vip")]
  res <- res %>% dplyr::left_join(vip_lookup, by = c("feature", "comparison"))
  res$vip <- round(res$vip, 4)

  res$regulation <- ifelse(
    is.na(res$p_adjust), "NS",
    ifelse(res$log2_fc > 0 & res$p_adjust < p_threshold, "Up",
      ifelse(res$log2_fc < 0 & res$p_adjust < p_threshold, "Down", "NS"))
  )
  reg_order <- c("Up", "Down", "NS")
  res$regulation <- factor(res$regulation, levels = reg_order)
  res <- res[order(res$regulation, -abs(res$log2_fc)), ]
  res$regulation <- as.character(res$regulation)
  res
}
