#' Perform t-Distributed Stochastic Neighbor Embedding (t-SNE) Analysis
#'
#' Performs t-SNE dimensionality reduction on high-dimensional data and provides
#' comprehensive results including embeddings, parameter diagnostics, and
#' visualization. t-SNE is a non-linear technique particularly effective at
#' revealing local structure and clusters in complex datasets.
#'
#' @param data Numeric data frame or matrix of abundance values. Rows are
#'   samples and columns are variables (features) by default; the transpose
#'   (features as rows, samples as columns) is also accepted — see
#'   \code{feature_as_row}. All variables should be numeric.
#' @param feature_as_row Logical or \code{NA}. \code{NA} (default) auto-detects
#'   the orientation by matching sample IDs from \code{sample} against the row
#'   and column names of \code{data}; \code{TRUE} forces features-as-rows;
#'   \code{FALSE} forces samples-as-rows. When detected or forced, the matrix is
#'   transposed internally so a manual \code{t()} is not needed.
#' @param sample Data frame containing sample metadata. Must include a
#'   \code{sample_id} column (or similar) matching the row names of \code{data}.
#' @param dims Output dimensionality (default: 2). Use 2 for standard 2D
#'   visualization, or 3 for 3D analysis.
#' @param perplexity Perplexity parameter controlling the balance between local
#'   and global aspects of the data (default: 30). Typical values range from 5
#'   to 50. Should be smaller when sample size is small (e.g., \code{floor((n - 1) / 3)}).
#' @param theta Barnes-Hut approximation parameter (default: 0.5). Speeds up
#'   computation at the cost of accuracy. Set to 0.0 for exact t-SNE (slower).
#' @param max_iter Maximum number of iterations (default: 1000). More iterations
#'   may improve convergence for complex datasets.
#' @param pca Logical; whether to perform initial PCA reduction before t-SNE
#'   (default: TRUE). Recommended when number of variables is large.
#' @param pca_dims Number of PCs to retain as input to t-SNE when
#'   \code{pca = TRUE} (default: 50). Ignored if \code{pca = FALSE}.
#' @param scale_data Logical; whether to scale variables to unit variance before
#'   running t-SNE (default: TRUE).
#' @param seed Random seed for reproducibility (default: 42). t-SNE results are
#'   stochastic; set a seed for reproducible outputs.
#' @param color_by Column name in \code{sample} for point colors (default: "group").
#' @param shape_by Column name in \code{sample} for point shapes (default: same as
#'   \code{color_by}).
#' @param plot_type Type of plots to generate (default: "all"). Options:
#'   \itemize{
#'     \item "all": Generate all available plots
#'     \item "scores": Only the t-SNE scatter plot
#'     \item "none": No plots (results only)
#'   }
#' @param conf_ellipses Logical; whether to add confidence ellipses around groups
#'   (default: FALSE).
#' @param ellipse_level Confidence level for ellipses (default: 0.95).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{tsne_model}}{Complete Rtsne model object}
#'     \item{\code{sample_scores}}{Data frame with t-SNE coordinates and sample metadata:
#'       \itemize{
#'         \item \code{sample_id}: Sample identifiers
#'         \item \code{tSNE1, tSNE2, ...}: t-SNE coordinates
#'         \item All columns from sample metadata
#'       }}
#'     \item{\code{plots}}{List of ggplot2 objects (when \code{plot_type != "none"}):
#'       \itemize{
#'         \item \code{score_plot}: t-SNE scatter plot with grouping
#'       }}
#'     \item{\code{summary_stats}}{List with key summary statistics:
#'       \itemize{
#'         \item \code{n_samples}: Number of samples
#'         \item \code{n_variables}: Number of input variables
#'         \item \code{perplexity}: Perplexity used
#'         \item \code{iterations}: Number of iterations run
#'         \item \code{final_cost}: Final KL divergence cost
#'       }}
#'   }
#'
#' @details
#' \strong{t-SNE Workflow:}
#' \enumerate{
#'   \item Optional PCA pre-processing to reduce dimensionality
#'   \item Compute pairwise affinities in high-dimensional space using Gaussian kernels
#'   \item Initialize low-dimensional embeddings
#'   \item Optimize embeddings via gradient descent to match affinity distributions
#' }
#'
#' \strong{Perplexity Guidelines:}
#' \itemize{
#'   \item Perplexity roughly corresponds to the number of nearest neighbors considered
#'   \item Typical range: 5--50; default is 30
#'   \item Must satisfy \code{perplexity < (n_samples - 1) / 3}
#'   \item Lower values emphasize fine local structure; higher values capture broader patterns
#' }
#'
#' @note
#' \itemize{
#'   \item t-SNE is stochastic; results vary between runs. Use \code{seed} for reproducibility.
#'   \item Cluster sizes and distances between clusters in t-SNE plots are not directly meaningful.
#'   \item For very large datasets, Barnes-Hut t-SNE (\code{theta > 0}) is recommended.
#'   \item t-SNE does not produce an explicit mapping; new data cannot be projected without re-running.
#' }
#'
#' @references
#' van der Maaten, L. and Hinton, G. (2008). Visualizing Data using t-SNE.
#' Journal of Machine Learning Research, 9(Nov), 2579--2605.
#'
#' van der Maaten, L. (2014). Accelerating t-SNE using Tree-Based Algorithms.
#' Journal of Machine Learning Research, 15(1), 3221--3245.
#'
#' @seealso
#' \code{\link{pca_analysis}} for linear dimensionality reduction,
#' \code{\link[Rtsne]{Rtsne}} for the underlying implementation
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # Prepare data and sample metadata
#' iris_data <- iris[, 1:4]
#' iris_samples <- data.frame(
#'   sample_id = paste0("sample_", 1:nrow(iris)),
#'   species = iris$Species
#' )
#' rownames(iris_data) <- iris_samples$sample_id
#'
#' # Basic t-SNE analysis
#' tsne_res <- tsne_analysis(
#'   data = iris_data,
#'   sample = iris_samples,
#'   color_by = "species",
#'   perplexity = 30
#' )
#'
#' # View results
#' head(tsne_res$sample_scores)
#' print(tsne_res$plots$score_plot)
#'
#' # Small sample size: lower perplexity
#' tsne_small <- tsne_analysis(
#'   data = iris_data[1:30, ],
#'   sample = iris_samples[1:30, ],
#'   perplexity = 5,
#'   max_iter = 500
#' )
#' }
tsne_analysis <- function(data, sample, dims = 2, perplexity = 30,
                          theta = 0.5, max_iter = 1000,
                          pca = TRUE, pca_dims = 50, scale_data = TRUE,
                          seed = 42,
                          color_by = "group", shape_by = NULL,
                          plot_type = "all",
                          conf_ellipses = FALSE, ellipse_level = 0.95,
                          feature_as_row = NA) {

  # --- Input validation ---
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix")
  }

  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }

  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("'data' cannot be empty")
  }

  # Resolve orientation before as.matrix() so data and data_matrix stay aligned
  data <- orient_to_sample_row(data, sample, NULL, feature_as_row, FALSE)

  data_matrix <- as.matrix(data)
  if (!is.numeric(data_matrix)) {
    stop("'data' must contain only numeric values")
  }

  # Sample matching
  data_samples <- rownames(data)
  if (is.null(data_samples)) {
    data_samples <- paste0("Sample_", 1:nrow(data))
    rownames(data) <- data_samples
  }

  if (!"sample_id" %in% names(sample)) {
    sample_id_cols <- c("sample", "Sample", "ID", "id", "sample_name")
    found_col <- intersect(sample_id_cols, names(sample))
    if (length(found_col) > 0) {
      sample$sample_id <- sample[[found_col[1]]]
    } else {
      stop("Sample metadata must contain a 'sample_id' column or similar identifier")
    }
  }

  if (!all(data_samples %in% sample$sample_id)) {
    missing_samples <- setdiff(data_samples, sample$sample_id)
    stop("Samples missing from metadata: ", paste(head(missing_samples, 5), collapse = ", "))
  }

  sample_aligned <- sample[match(data_samples, sample$sample_id), ]

  # Parameter validation
  if (!is.numeric(dims) || length(dims) != 1 || dims < 1 || dims != round(dims)) {
    stop("'dims' must be a positive integer")
  }

  if (!is.numeric(perplexity) || length(perplexity) != 1 || perplexity <= 0) {
    stop("'perplexity' must be a positive number")
  }

  n_samples <- nrow(data_matrix)
  max_perp <- floor((n_samples - 1) / 3)
  if (perplexity >= max_perp) {
    perplexity <- max_perp
    warning("Perplexity too large for ", n_samples, " samples, reduced to ", perplexity)
  }

  if (is.null(shape_by)) {
    shape_by <- color_by
  } else if (!shape_by %in% names(sample_aligned)) {
    stop("'shape_by' column '", shape_by, "' not found in sample metadata")
  }

  if (!color_by %in% names(sample_aligned)) {
    stop("'color_by' column '", color_by, "' not found in sample metadata")
  }

  valid_plot_types <- c("all", "scores", "none")
  if (!plot_type %in% valid_plot_types) {
    stop("'plot_type' must be one of: ", paste(valid_plot_types, collapse = ", "))
  }

  # --- Data preprocessing ---
  # Remove zero-variance variables
  var_check <- apply(data_matrix, 2, var, na.rm = TRUE)
  keep_vars <- var_check > 0
  if (sum(!keep_vars) > 0) {
    warning("Removed ", sum(!keep_vars), " zero-variance variable(s)")
    data_matrix <- data_matrix[, keep_vars, drop = FALSE]
  }

  # Scale if requested
  if (scale_data) {
    data_matrix <- scale(data_matrix)
  }

  # Handle NA
  if (any(is.na(data_matrix))) {
    na_count <- sum(is.na(data_matrix))
    warning("Data contains ", na_count, " missing values, replacing with column means")
    for (j in seq_len(ncol(data_matrix))) {
      data_matrix[is.na(data_matrix[, j]), j] <- mean(data_matrix[, j], na.rm = TRUE)
    }
  }

  # --- Run t-SNE ---
  if (!is.null(seed)) {
    set.seed(seed)
  }

  tsne_model <- tryCatch({
    Rtsne::Rtsne(
      X = data_matrix,
      dims = dims,
      perplexity = perplexity,
      theta = theta,
      max_iter = max_iter,
      pca = pca,
      initial_dims = if (pca) min(pca_dims, ncol(data_matrix)) else ncol(data_matrix),
      check_duplicates = FALSE,
      verbose = FALSE
    )
  }, error = function(e) {
    stop("t-SNE analysis failed: ", e$message)
  })

  # --- Extract results ---
  tsne_coords <- as.data.frame(tsne_model$Y)
  coord_names <- paste0("tSNE", 1:dims)
  names(tsne_coords) <- coord_names
  tsne_coords$sample_id <- data_samples

  sample_scores <- tsne_coords %>%
    dplyr::left_join(sample_aligned, by = "sample_id")

  # --- Plots ---
  plots <- list()

  if (plot_type %in% c("all", "scores") && dims >= 2) {
    score_plot <- sample_scores %>%
      ggplot2::ggplot(ggplot2::aes(
        x = .data$tSNE1,
        y = .data$tSNE2,
        color = !!rlang::sym(color_by),
        shape = !!rlang::sym(shape_by)
      )) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        x = "tSNE1",
        y = "tSNE2",
        title = "t-SNE Plot",
        subtitle = paste0("Perplexity: ", perplexity,
                          " | Iterations: ", max_iter,
                          " | Samples: ", n_samples),
        color = tools::toTitleCase(gsub("[._]", " ", color_by)),
        shape = tools::toTitleCase(gsub("[._]", " ", shape_by))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11),
        axis.title = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 11)
      )

    if (conf_ellipses) {
      score_plot <- score_plot +
        ggplot2::stat_ellipse(
          ggplot2::aes(color = !!rlang::sym(color_by)),
          level = ellipse_level, type = "norm", alpha = 0.3
        )
    }

    plots$score_plot <- score_plot
  }

  # --- Summary statistics ---
  summary_stats <- list(
    n_samples = n_samples,
    n_variables = ncol(data),
    n_variables_used = ncol(data_matrix),
    dims = dims,
    perplexity = perplexity,
    theta = theta,
    iterations = max_iter,
    pca_preprocessing = pca,
    pca_dims_used = if (pca) min(pca_dims, ncol(data_matrix)) else NA,
    scaled = scale_data,
    seed = seed,
    final_cost = if (!is.null(tsne_model$itercosts)) {
      tail(tsne_model$itercosts, 1)
    } else NA
  )

  # --- Return ---
  results <- list(
    tsne_model = tsne_model,
    sample_scores = sample_scores,
    plots = if (plot_type != "none") plots else NULL,
    summary_stats = summary_stats
  )

  attr(results, "analysis_type") <- "t-SNE"
  attr(results, "n_samples") <- n_samples
  attr(results, "n_variables") <- ncol(data)
  attr(results, "perplexity") <- perplexity

  if (interactive()) {
    cat("t-SNE Analysis Summary:\n")
    cat("=======================\n")
    cat("Samples:", n_samples, "\n")
    cat("Variables:", ncol(data), "\n")
    cat("Output dimensions:", dims, "\n")
    cat("Perplexity:", perplexity, "\n")
    cat("Iterations:", max_iter, "\n")
    cat("PCA preprocessing:", if (pca) "Yes" else "No", "\n")
    cat("Data scaling:", if (scale_data) "Yes" else "No", "\n")
    if (!is.na(summary_stats$final_cost)) {
      cat("Final cost (KL divergence):", round(summary_stats$final_cost, 4), "\n")
    }
    cat("\n")
  }

  results
}
