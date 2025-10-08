#' Identify Core Microbiome
#'
#' @description
#' Identify core microbiome members using occurrence frequency, relative abundance,
#' and neutral model fitting. Combines multiple approaches to robustly define
#' which taxa constitute the core microbiome.
#'
#' @param otu_table A numeric matrix or data frame with OTUs/ASVs as rows and
#'   samples as columns. Row names should be OTU/ASV identifiers.
#' @param metadata A data frame containing sample metadata with a column
#'   matching \code{sample_col} and a grouping column specified by \code{group_col}.
#' @param sample_col Character string. Name of the column in \code{metadata}
#'   that contains sample identifiers matching column names in \code{otu_table}.
#' @param group_col Character string. Name of the column in \code{metadata}
#'   used for grouping samples (e.g., "treatment", "timepoint"). Default is "group".
#' @param n_ranked Integer. Number of top-ranked OTUs to evaluate for core
#'   membership. Default is 500.
#' @param threshold Numeric. Threshold for determining core cutoff using the
#'   "last increase" method. Represents the minimum proportional increase in
#'   Bray-Curtis dissimilarity (e.g., 0.02 = 2% increase). Default is 0.02.
#'
#' @details
#' This function implements a multi-step approach to identify core microbiome members:
#'
#' \enumerate{
#'   \item \strong{Ranking}: OTUs are ranked by an index combining occurrence
#'     frequency and consistency across groups
#'   \item \strong{Contribution}: Calculates each OTU's contribution to
#'     community dissimilarity using Bray-Curtis distances
#'   \item \strong{Cutoff determination}: Uses two methods:
#'     \itemize{
#'       \item Elbow method (first-order difference)
#'       \item Last increase method (based on threshold parameter)
#'     }
#'   \item \strong{Neutral model}: Fits Sloan's neutral model to assess
#'     whether taxa follow neutral expectations
#' }
#'
#' @return A list with class "core_microbiome" containing:
#'   \item{elbow}{Integer. Core size determined by elbow method}
#'   \item{last_call}{Integer. Core size determined by last increase method}
#'   \item{neutral_model}{Data frame with neutral model predictions}
#'   \item{ranked_otus}{Data frame of ranked OTUs with dissimilarity metrics}
#'   \item{fitting_results}{Data frame combining ranking and neutral model results}
#'   \item{plot_ranked}{ggplot object showing ranked OTUs}
#'   \item{plot_neutral}{ggplot object showing neutral model fit}
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @importFrom dplyr mutate group_by summarise ungroup left_join arrange desc
#'   transmute select filter case_when slice
#' @importFrom tidyr pivot_longer gather
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_vline annotate
#'   scale_y_continuous labs theme element_blank element_text
#' @importFrom magrittr %>%
#' @importFrom rlang .data := !!
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' otu_data <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' rownames(otu_data) <- paste0("OTU", 1:100)
#' colnames(otu_data) <- paste0("Sample", 1:10)
#'
#' metadata <- data.frame(
#'   sample_id = paste0("Sample", 1:10),
#'   group = rep(c("A", "B"), each = 5)
#' )
#'
#' # Identify core microbiome
#' core_result <- identify_core_microbiome(
#'   otu_table = otu_data,
#'   metadata = metadata,
#'   sample_col = "sample_id",
#'   group_col = "group",
#'   n_ranked = 50,
#'   threshold = 0.02
#' )
#'
#' # View results
#' print(core_result)
#' core_result$plot_ranked
#' core_result$plot_neutral
#' }
#'
#' @export
identify_core_microbiome <- function(otu_table,
                                      metadata,
                                      sample_col,
                                      group_col = "group",
                                      n_ranked = 500,
                                      threshold = 0.02) {
  # Input validation
  validate_core_inputs(otu_table, metadata, sample_col, group_col, n_ranked, threshold)
  
  # Convert to data frames if needed
  otu_df <- as.data.frame(otu_table)
  meta_df <- as.data.frame(metadata)
  
  # Calculate OTU occurrence and abundance metrics
  otu_metrics <- calculate_otu_metrics(otu_df)
  
  # Rank OTUs by occurrence and consistency
  otu_ranking <- rank_otus_by_occurrence(
    otu_df,
    meta_df,
    sample_col,
    group_col
  )
  
  # Calculate contribution to community dissimilarity
  bc_ranked <- calculate_bray_curtis_contribution(
    otu_df,
    otu_ranking,
    n_ranked
  )
  
  # Determine core cutoffs
  elbow <- determine_elbow_cutoff(bc_ranked)
  last_call <- determine_last_increase_cutoff(bc_ranked, threshold)
  
  # Fit neutral model
  neutral_results <- fit_neutral_model(otu_df)
  
  # Combine results
  fitting_results <- combine_ranking_and_neutral(
    otu_metrics,
    otu_ranking,
    last_call
  )
  
  # Create plots
  plot_ranked <- plot_ranked_otus(bc_ranked, elbow, last_call, threshold)
  plot_neutral <- plot_neutral_model(fitting_results, neutral_results)
  
  # Create output object
  result <- structure(
    list(
      elbow = elbow,
      last_call = last_call,
      neutral_model = neutral_results,
      ranked_otus = bc_ranked,
      fitting_results = fitting_results,
      plot_ranked = plot_ranked,
      plot_neutral = plot_neutral
    ),
    class = "core_microbiome"
  )
  
  return(result)
}


#' Validate Core Microbiome Inputs
#' @keywords internal
#' @noRd
validate_core_inputs <- function(otu_table, metadata, sample_col,
                                  group_col, n_ranked, threshold) {
  # Check OTU table
  if (!is.matrix(otu_table) && !is.data.frame(otu_table)) {
    stop("otu_table must be a matrix or data frame", call. = FALSE)
  }
  
  if (is.null(rownames(otu_table)) || any(rownames(otu_table) == "")) {
    stop("otu_table must have row names (OTU identifiers)", call. = FALSE)
  }
  
  # Check metadata
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame", call. = FALSE)
  }
  
  if (!sample_col %in% names(metadata)) {
    stop(
      "sample_col '", sample_col, "' not found in metadata",
      call. = FALSE
    )
  }
  
  if (!group_col %in% names(metadata)) {
    stop(
      "group_col '", group_col, "' not found in metadata",
      call. = FALSE
    )
  }
  
  # Check sample matching
  otu_samples <- colnames(otu_table)
  meta_samples <- metadata[[sample_col]]
  
  if (!all(otu_samples %in% meta_samples)) {
    stop(
      "Not all OTU table samples found in metadata ",
      sample_col,
      call. = FALSE
    )
  }
  
  # Check parameters
  if (!is.numeric(n_ranked) || n_ranked < 1) {
    stop("n_ranked must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be between 0 and 1", call. = FALSE)
  }
  
  if (n_ranked > nrow(otu_table)) {
    warning(
      "n_ranked (", n_ranked, ") exceeds number of OTUs (",
      nrow(otu_table), "). Using all OTUs.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Calculate OTU Occurrence and Abundance Metrics
#' @keywords internal
#' @noRd
calculate_otu_metrics <- function(otu_df) {
  n_samples <- ncol(otu_df)
  
  otu_df %>%
    tibble::rownames_to_column(var = "otu") %>%
    tidyr::pivot_longer(
      cols = -"otu",
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    dplyr::mutate(
      present = dplyr::if_else(.data$abundance > 0, 1, 0)
    ) %>%
    dplyr::group_by(.data$otu) %>%
    dplyr::mutate(
      otu_occurrences = sum(.data$present)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::mutate(
      sample_total = sum(.data$abundance),
      rel_abundance = .data$abundance / .data$sample_total
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$otu) %>%
    dplyr::mutate(
      otu_occurrence_freq = .data$otu_occurrences / n_samples,
      otu_mean_rel_abun = mean(.data$rel_abundance)
    ) %>%
    dplyr::ungroup()
}


#' Rank OTUs by Occurrence Pattern
#' @keywords internal
#' @noRd
rank_otus_by_occurrence <- function(otu_df, meta_df, sample_col, group_col) {
  group_sym <- rlang::sym(group_col)
  
  otu_df %>%
    tibble::rownames_to_column(var = "otu") %>%
    dplyr::mutate(otu = factor(.data$otu)) %>%
    tidyr::pivot_longer(
      cols = -"otu",
      names_to = sample_col,
      values_to = "abundance"
    ) %>%
    dplyr::left_join(meta_df, by = sample_col) %>%
    dplyr::group_by(.data$otu, !!group_sym) %>%
    dplyr::summarise(
      freq = sum(.data$abundance > 0) / dplyr::n(),
      is_core = dplyr::if_else(.data$freq == 1, 1, 0),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$otu) %>%
    dplyr::summarise(
      sum_freq = sum(.data$freq),
      sum_core = sum(.data$is_core),
      n_groups = dplyr::n(),
      index = (.data$sum_freq + .data$sum_core) / .data$n_groups,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$index))
}


#' Calculate Bray-Curtis Contribution of Ranked OTUs
#' @keywords internal
#' @noRd
calculate_bray_curtis_contribution <- function(otu_df, otu_ranking, n_ranked) {
  n_otus_to_use <- min(n_ranked, nrow(otu_ranking))
  n_samples <- ncol(otu_df)
  sample_names <- colnames(otu_df)
  
  # Initialize results
  bc_accumulation <- list()
  
  # Calculate for each rank
  for (i in seq_len(n_otus_to_use)) {
    selected_otus <- as.character(otu_ranking$otu[1:i])
    subset_matrix <- as.matrix(otu_df[selected_otus, , drop = FALSE])
    
    # Calculate pairwise Bray-Curtis
    bc_values <- calculate_pairwise_bc(subset_matrix)
    bc_accumulation[[i]] <- bc_values
    
    if (i %% 50 == 0) {
      message("Processed ", i, "/", n_otus_to_use, " OTUs")
    }
  }
  
  # Calculate for full community
  full_bc <- calculate_pairwise_bc(as.matrix(otu_df))
  
  # Summarize results
  bc_summary <- tibble::tibble(
    rank = factor(seq_len(n_otus_to_use)),
    mean_bc = sapply(bc_accumulation, mean)
  ) %>%
    dplyr::mutate(
      proportion_bc = .data$mean_bc / max(.data$mean_bc)
    )
  
  # Calculate increase
  increase <- c(0, bc_summary$mean_bc[-1] / bc_summary$mean_bc[-nrow(bc_summary)])
  bc_summary$increase_bc <- increase
  
  # Calculate first-order differences for elbow method
  bc_summary$fo_diff <- sapply(
    seq_len(nrow(bc_summary)),
    function(pos) {
      left <- (bc_summary$mean_bc[pos] - bc_summary$mean_bc[1]) / pos
      right <- (bc_summary$mean_bc[nrow(bc_summary)] - bc_summary$mean_bc[pos]) /
        (nrow(bc_summary) - pos)
      return(left - right)
    }
  )
  
  return(bc_summary)
}


#' Calculate Pairwise Bray-Curtis Dissimilarity
#' @keywords internal
#' @noRd
calculate_pairwise_bc <- function(matrix) {
  n_samples <- ncol(matrix)
  sample_totals <- colSums(matrix)
  
  bc_values <- numeric()
  
  for (i in seq_len(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      bc <- sum(abs(matrix[, i] - matrix[, j])) / (2 * sample_totals[1])
      bc_values <- c(bc_values, bc)
    }
  }
  
  return(bc_values)
}


#' Determine Core Size Using Elbow Method
#' @keywords internal
#' @noRd
determine_elbow_cutoff <- function(bc_ranked) {
  which.max(bc_ranked$fo_diff)
}


#' Determine Core Size Using Last Increase Method
#' @keywords internal
#' @noRd
determine_last_increase_cutoff <- function(bc_ranked, threshold) {
  above_threshold <- bc_ranked %>%
    dplyr::filter(.data$increase_bc >= 1 + threshold)
  
  if (nrow(above_threshold) == 0) {
    warning(
      "No OTUs met the threshold increase of ",
      threshold * 100,
      "%. Using elbow method instead.",
      call. = FALSE
    )
    return(determine_elbow_cutoff(bc_ranked))
  }
  
  as.numeric(as.character(dplyr::last(above_threshold$rank)))
}


#' Combine Ranking Results with Neutral Model
#' @keywords internal
#' @noRd
combine_ranking_and_neutral <- function(otu_metrics, otu_ranking, last_call) {
  core_otus <- as.character(otu_ranking$otu[seq_len(last_call)])
  
  otu_metrics %>%
    dplyr::select(
      "otu",
      "otu_mean_rel_abun",
      "otu_occurrence_freq"
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      is_core = dplyr::if_else(.data$otu %in% core_otus, "Core", "Non-core")
    )
}


#' Plot Ranked OTUs
#' @keywords internal
#' @noRd
plot_ranked_otus <- function(bc_ranked, elbow, last_call, threshold) {
  thread_value <- bc_ranked[bc_ranked$rank == last_call, ]$proportion_bc
  
  plot_data <- bc_ranked %>%
    dplyr::mutate(
      rank = as.numeric(as.character(.data$rank)),
      status = dplyr::if_else(
        .data$proportion_bc < thread_value,
        "Core",
        "Non-core"
      )
    ) %>%
    dplyr::slice(1:min(100, dplyr::n()))
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$rank, y = .data$proportion_bc)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$status), size = 2) +
    ggplot2::geom_vline(
      xintercept = elbow,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = last_call,
      linetype = "dashed",
      color = "blue",
      linewidth = 1
    ) +
    ggplot2::annotate(
      "text",
      x = elbow + 10,
      y = 0.1,
      label = paste0("Elbow method (", elbow, ")"),
      color = "red"
    ) +
    ggplot2::annotate(
      "text",
      x = last_call + 10,
      y = 0.5,
      label = paste0("Last ", threshold * 100, "% increase (", last_call, ")"),
      color = "blue"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    ggplot2::labs(
      x = "Ranked OTUs",
      y = "Proportion of Bray-Curtis Dissimilarity",
      color = "Status"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = c(0.8, 0.2),
      legend.text = ggplot2::element_text(size = 12)
    )
}


#' Plot Neutral Model Fit
#' @keywords internal
#' @noRd
plot_neutral_model <- function(fitting_results, neutral_results) {
  ggplot2::ggplot(fitting_results) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = log10(.data$otu_mean_rel_abun),
        y = .data$otu_occurrence_freq,
        color = .data$is_core
      ),
      size = 3,
      alpha = 0.7
    ) +
    ggplot2::geom_line(
      data = neutral_results,
      ggplot2::aes(
        x = log10(.data$mean_rel_abun),
        y = .data$freq_pred
      ),
      color = "black",
      linewidth = 1
    ) +
    ggplot2::geom_line(
      data = neutral_results,
      ggplot2::aes(
        x = log10(.data$mean_rel_abun),
        y = .data$pred_upper
      ),
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = neutral_results,
      ggplot2::aes(
        x = log10(.data$mean_rel_abun),
        y = .data$pred_lower
      ),
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    ggplot2::labs(
      x = "Log10(Mean Relative Abundance)",
      y = "Occurrence Frequency",
      color = "Status"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = c(0.2, 0.8),
      legend.text = ggplot2::element_text(size = 12)
    )
}


#' Print Method for Core Microbiome Results
#' @export
print.core_microbiome <- function(x, ...) {
  cat("Core Microbiome Analysis Results\n")
  cat("=================================\n\n")
  cat("Core size (elbow method):", x$elbow, "OTUs\n")
  cat("Core size (last increase):", x$last_call, "OTUs\n\n")
  cat("Use plot(x$plot_ranked) to visualize ranked OTUs\n")
  cat("Use plot(x$plot_neutral) to visualize neutral model fit\n")
  invisible(x)
}


#' Fit Neutral Model to Community Data
#'
#' @description
#' Fits Sloan's neutral community model to OTU occurrence data. The neutral
#' model predicts species occurrence frequencies based on their abundance in
#' the metacommunity and dispersal limitation.
#'
#' @param otu_table A numeric matrix with OTUs as rows and samples as columns.
#'   Should ideally be rarefied to equal sequencing depth.
#' @param pool Optional matrix defining the source community pool. If NULL
#'   (default), uses the mean across all samples.
#' @param return_stats Logical. If TRUE, returns model fitting statistics.
#'   If FALSE, returns observed and predicted values for each OTU. Default is FALSE.
#'
#' @details
#' The neutral model from Sloan et al. (2006) predicts occurrence frequencies
#' based on:
#' \itemize{
#'   \item Mean relative abundance in the metacommunity (p)
#'   \item Migration rate parameter (m)
#'   \item Community size (N)
#' }
#'
#' The model assumes that species dynamics are driven by random birth, death,
#' and dispersal, with no selective differences between species.
#'
#' @return If \code{return_stats = TRUE}, returns a data frame with one row
#'   containing model fit statistics (R², AIC, BIC, migration parameter m).
#'   If \code{return_stats = FALSE}, returns a data frame with columns:
#'   \item{otu}{OTU identifier}
#'   \item{mean_rel_abun}{Mean relative abundance}
#'   \item{freq_observed}{Observed occurrence frequency}
#'   \item{freq_pred}{Predicted occurrence frequency}
#'   \item{pred_lower}{Lower 95% prediction interval}
#'   \item{pred_upper}{Upper 95% prediction interval}
#'
#' @references
#' Sloan, W. T. et al. (2006). Quantifying the roles of immigration and chance
#' in shaping prokaryote community structure. Environmental Microbiology, 8(4), 732-740.
#'
#' Burns, A. R. et al. (2016). Contribution of neutral processes to the assembly
#' of gut microbial communities in the zebrafish over host development.
#' The ISME Journal, 10(3), 655-664.
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @importFrom minpack.lm nlsLM
#' @importFrom Hmisc binconf
#' @importFrom stats4 mle AIC BIC
#' @importFrom stats pbeta pbinom ppois dnorm
#'
#' @examples
#' \dontrun{
#' # Create example data
#' otu_data <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' rownames(otu_data) <- paste0("OTU", 1:100)
#'
#' # Fit model and get predictions
#' predictions <- fit_sloan_neutral_model(otu_data, return_stats = FALSE)
#'
#' # Get model statistics
#' stats <- fit_sloan_neutral_model(otu_data, return_stats = TRUE)
#' print(stats)
#' }
#'
#' @export
fit_sloan_neutral_model <- function(otu_table,
                                     pool = NULL,
                                     return_stats = FALSE) {
  # Validate inputs
  if (!is.matrix(otu_table) && !is.data.frame(otu_table)) {
    stop("otu_table must be a matrix or data frame", call. = FALSE)
  }
  
  # Transpose if needed (samples should be rows)
  if (is.null(rownames(otu_table))) {
    stop("otu_table must have row names (OTU identifiers)", call. = FALSE)
  }
  
  spp <- t(as.matrix(otu_table))
  
  # Calculate mean community size
  N <- mean(rowSums(spp))
  
  # Calculate mean relative abundance
  if (is.null(pool)) {
    p_mean <- colMeans(spp)
    p_mean <- p_mean[p_mean != 0]
    p <- p_mean / N
  } else {
    p_mean <- colMeans(pool)
    p_mean <- p_mean[p_mean != 0]
    p <- p_mean / N
  }
  
  # Calculate occurrence frequencies
  spp_binary <- (spp > 0) * 1
  freq <- colMeans(spp_binary)
  freq <- freq[freq != 0]
  
  # Combine and remove zeros
  combined <- merge(
    data.frame(otu = names(p), p = p, row.names = NULL),
    data.frame(otu = names(freq), freq = freq, row.names = NULL),
    by = "otu"
  )
  combined <- combined[order(combined$p), ]
  
  p <- combined$p
  freq <- combined$freq
  names(p) <- combined$otu
  names(freq) <- combined$otu
  
  # Detection limit
  d <- 1 / N
  
  # Fit migration parameter using non-linear least squares
  tryCatch(
    {
      m_fit <- minpack.lm::nlsLM(
        freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE),
        start = list(m = 0.1)
      )
      
      m_ci <- confint(m_fit, "m", level = 0.95)
      
      # Maximum likelihood estimation
      sncm_ll <- function(m, sigma) {
        R <- freq - pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
        R <- dnorm(R, 0, sigma)
        -sum(log(R))
      }
      
      m_mle <- stats4::mle(
        sncm_ll,
        start = list(m = 0.1, sigma = 0.1),
        nobs = length(p)
      )
      
      # Model fit statistics
      aic_fit <- AIC(m_mle, k = 2)
      bic_fit <- BIC(m_mle)
      
      # Predictions
      freq_pred <- pbeta(
        d,
        N * coef(m_fit) * p,
        N * coef(m_fit) * (1 - p),
        lower.tail = FALSE
      )
      
      # R-squared
      rsqr <- 1 - sum((freq - freq_pred)^2) / sum((freq - mean(freq))^2)
      
      # RMSE
      rmse <- sqrt(sum((freq - freq_pred)^2) / (length(freq) - 1))
      
      # Prediction intervals
      pred_ci <- Hmisc::binconf(
        freq_pred * nrow(spp),
        nrow(spp),
        alpha = 0.05,
        method = "wilson",
        return.df = TRUE
      )
      
      if (return_stats) {
        # Return statistics
        return(
          data.frame(
            m = coef(m_fit),
            m_ci = coef(m_fit) - m_ci[1],
            m_mle = m_mle@coef["m"],
            rsquared = rsqr,
            rmse = rmse,
            aic = aic_fit,
            bic = bic_fit,
            N = N,
            n_samples = nrow(spp),
            richness = length(p),
            detection_limit = d
          )
        )
      } else {
        # Return predictions
        result <- data.frame(
          otu = names(p),
          mean_rel_abun = p,
          freq_observed = freq,
          freq_pred = freq_pred,
          pred_lower = pred_ci[, 2],
          pred_upper = pred_ci[, 3],
          row.names = NULL
        )
        return(result)
      }
    },
    error = function(e) {
      warning(
        "Neutral model fitting failed: ",
        e$message,
        "\nReturning NULL",
        call. = FALSE
      )
      return(NULL)
    }
  )
}


#' Internal wrapper for neutral model fitting
#' @keywords internal
#' @noRd
fit_neutral_model <- function(otu_df) {
  result <- fit_sloan_neutral_model(
    otu_table = otu_df,
    pool = NULL,
    return_stats = FALSE
  )
  
  if (is.null(result)) {
    warning("Neutral model fitting failed, returning empty data frame", call. = FALSE)
    return(data.frame(
      otu = character(),
      mean_rel_abun = numeric(),
      freq_observed = numeric(),
      freq_pred = numeric(),
      pred_lower = numeric(),
      pred_upper = numeric()
    ))
  }
  
  return(result)
}


# Global variables
utils::globalVariables(c(".", ".data"))