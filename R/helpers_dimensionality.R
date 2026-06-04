#' Dimensionality Reduction Analysis Helper Functions
#'
#' Shared utility functions for PCA, PCoA, and RDA analyses.
#' These functions extract common patterns in dimensionality reduction methods,
#' reducing code duplication and improving maintainability.
#'
#' @details
#' Functions provided:
#' \itemize{
#'   \item \code{extract_eigenvalues()}: Extract and format eigenvalues/explained variance
#'   \item \code{prepare_ordination_data()}: Extract and format ordination coordinates
#'   \item \code{calculate_variance_explained()}: Compute variance explained percentages
#'   \item \code{create_scree_plot()}: Generate variance explained visualization
#'   \item \code{create_biplot_arrows()}: Format biplot arrows for environmental variables
#'   \item \code{validate_ordination_inputs()}: Check ordination analysis prerequisites
#' }
#'
#' @keywords internal

#' Extract Eigenvalues from Ordination Results
#'
#' Standardizes eigenvalue extraction across different ordination methods.
#' Handles PCA (FactoMineR), PCoA (ape), and RDA (vegan) formats.
#'
#' @param eigen_object The eigenvalue object from an ordination analysis
#' @param analysis_type Character: "pca", "pcoa", or "rda"
#' @param n_components Maximum number of components to return
#'
#' @return Data frame with columns: component, eigenvalue, variance_percent, cumulative_percent
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # PCA eigenvalues
#' pca_eig <- extract_eigenvalues(pca_result$eig, "pca")
#'
#' # PCoA eigenvalues
#' pcoa_eig <- extract_eigenvalues(pcoa_result$values, "pcoa")
#' }
extract_eigenvalues <- function(eigen_object, analysis_type = "pca", n_components = NULL) {
  if (is.null(eigen_object)) {
    stop("Eigenvalue object cannot be NULL")
  }

  # Handle different eigenvalue formats
  if (analysis_type == "pca") {
    # FactoMineR returns $eig with eigenvalues in first column
    if (is.matrix(eigen_object)) {
      eigs <- as.numeric(eigen_object[, 1])
    } else {
      eigs <- as.numeric(eigen_object)
    }
  } else if (analysis_type == "pcoa") {
    # ape::pcoa returns data frame with Relative_eig in second column
    if (is.data.frame(eigen_object)) {
      eigs <- as.numeric(eigen_object[, 2])
    } else {
      eigs <- as.numeric(eigen_object)
    }
  } else if (analysis_type == "rda") {
    # vegan::rda eigenvalues are numeric vector
    eigs <- as.numeric(eigen_object)
  } else {
    stop("Unknown analysis_type: ", analysis_type)
  }

  # Remove negative eigenvalues and NAs
  eigs <- eigs[eigs > 0 & !is.na(eigs)]

  if (length(eigs) == 0) {
    stop("No valid eigenvalues found")
  }

  # Limit to n_components if specified
  if (!is.null(n_components)) {
    eigs <- eigs[seq_len(min(n_components, length(eigs)))]
  }

  # Calculate variance explained percentages
  total_variance <- sum(eigs)
  variance_pct <- (eigs / total_variance) * 100
  cumulative_pct <- cumsum(variance_pct)

  # Create component names based on type
  if (analysis_type == "pca") {
    comp_names <- paste0("PC", seq_len(length(eigs)))
  } else if (analysis_type == "pcoa") {
    comp_names <- paste0("PCo", seq_len(length(eigs)))
  } else {
    comp_names <- paste0("RDA", seq_len(length(eigs)))
  }

  # Return formatted data frame
  data.frame(
    component = comp_names,
    eigenvalue = round(eigs, 6),
    variance_percent = round(variance_pct, 2),
    cumulative_percent = round(cumulative_pct, 2),
    stringsAsFactors = FALSE
  )
}

#' Prepare Ordination Coordinate Data
#'
#' Standardizes extraction and formatting of sample coordinates from ordination results.
#' Merges coordinates with sample metadata.
#'
#' @param coord_matrix Matrix or data frame of sample coordinates (samples × components)
#' @param sample_data Data frame with sample metadata (rows = samples)
#' @param analysis_type Character: "pca", "pcoa", or "rda"
#' @param n_components Maximum number of components to include
#'
#' @return Data frame with sample IDs, coordinates, and metadata merged
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # PCA scores
#' pca_coords <- extract_from_pca_result$scores
#' prepared <- prepare_ordination_data(pca_coords, sample_metadata, "pca")
#' }
prepare_ordination_data <- function(coord_matrix, sample_data, analysis_type = "pca",
                                     n_components = NULL) {
  if (!is.data.frame(sample_data) && !is.matrix(sample_data)) {
    stop("sample_data must be a data frame or matrix")
  }

  # Convert to data frame if needed
  if (is.matrix(coord_matrix)) {
    coord_df <- as.data.frame(coord_matrix)
  } else {
    coord_df <- coord_matrix
  }

  # Limit to n_components if specified
  if (!is.null(n_components)) {
    coord_df <- coord_df[, seq_len(min(n_components, ncol(coord_df))), drop = FALSE]
  }

  # Create component names
  n_cols <- ncol(coord_df)
  if (analysis_type == "pca") {
    comp_names <- paste0("PC", seq_len(n_cols))
  } else if (analysis_type == "pcoa") {
    comp_names <- paste0("PCo", seq_len(n_cols))
  } else {
    comp_names <- paste0("RDA", seq_len(n_cols))
  }
  colnames(coord_df) <- comp_names

  # Add sample IDs from row names
  sample_id <- rownames(coord_df)
  if (is.null(sample_id)) {
    sample_id <- seq_len(nrow(coord_df))
  }

  # Create result data frame
  result <- cbind(sample_id = sample_id, coord_df, stringsAsFactors = FALSE)

  # Merge with sample metadata if first column is sample ID
  if (nrow(sample_data) == nrow(result)) {
    # Check if first column of sample_data contains sample IDs
    if ("sample" %in% colnames(sample_data) || colnames(sample_data)[1] == "sample") {
      sample_data_renamed <- sample_data
      if (colnames(sample_data)[1] != "sample_id") {
        colnames(sample_data_renamed)[1] <- "sample_id"
      }
      result <- merge(result, sample_data_renamed, by = "sample_id", all.x = TRUE)
    } else {
      # Direct merge if row order matches
      result <- cbind(result, sample_data)
    }
  }

  result
}

#' Calculate Variance Explained Statistics
#'
#' Computes summary statistics for variance explained by ordination components.
#'
#' @param eigenvalues Data frame from \code{extract_eigenvalues()}
#' @param n_components_to_sum Number of components to sum for total variance
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{total_variance}: Total variance explained by first n components
#'     \item \code{kaiser_criterion}: Number of components with eigenvalue > 1
#'     \item \code{broken_stick}: Components above broken stick model
#'   }
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' eigs <- extract_eigenvalues(pca_result$eig, "pca")
#' stats <- calculate_variance_explained(eigs, n_components_to_sum = 3)
#' cat("Total variance in PC1-PC3:", stats$total_variance, "%\n")
#' }
calculate_variance_explained <- function(eigenvalues, n_components_to_sum = NULL) {
  if (is.null(n_components_to_sum)) {
    n_components_to_sum <- min(2, nrow(eigenvalues))
  }

  # Total variance explained
  total_var <- sum(eigenvalues$variance_percent[seq_len(
    min(n_components_to_sum, nrow(eigenvalues))
  )])

  # Kaiser criterion (eigenvalue > 1, only for PCA)
  # Estimate from variance percent: eigenvalue > 1 means > 1/n_features percent variance
  kaiser_count <- sum(eigenvalues$eigenvalue > 1)

  # Broken stick model (theoretical expectation)
  n <- nrow(eigenvalues)
  j <- seq_len(n)
  broken_stick <- (100 / n) * sum(1 / j)

  list(
    total_variance = round(total_var, 2),
    kaiser_criterion = kaiser_count,
    broken_stick = round(broken_stick, 2)
  )
}

#' Create Scree Plot
#'
#' Generates a publication-ready scree plot showing variance explained by components.
#'
#' @param eigenvalues Data frame from \code{extract_eigenvalues()}
#' @param title Plot title (default: "Scree Plot")
#' @param theme_type ggplot2 theme to apply (default: "bw")
#'
#' @return ggplot2 object
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' eigs <- extract_eigenvalues(pca_result$eig, "pca")
#' scree <- create_scree_plot(eigs)
#' print(scree)
#' }
create_scree_plot <- function(eigenvalues, title = "Scree Plot", theme_type = "bw") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  # Get theme function
  theme_fn <- switch(theme_type,
    "bw" = ggplot2::theme_bw(),
    "minimal" = ggplot2::theme_minimal(),
    "classic" = ggplot2::theme_classic(),
    ggplot2::theme_bw()
  )

  # Create scree plot
  p <- ggplot2::ggplot(eigenvalues, ggplot2::aes(x = .data$component,
                                                   y = .data$variance_percent)) +
    ggplot2::geom_point(size = 3, color = "#0072B2") +
    ggplot2::geom_line(color = "#0072B2", group = 1, linewidth = 0.8) +
    ggplot2::labs(
      title = title,
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

  p
}

#' Format Environmental Variable Arrows for Biplot
#'
#' Prepares environmental variable scores for visualization as arrows in ordination biplots.
#'
#' @param env_scores Matrix or data frame of environmental variable coordinates
#' @param arrow_scale Scaling factor for arrow lengths (default: 1)
#' @param min_length Minimum arrow length to display (filters short arrows)
#'
#' @return Data frame with arrow coordinates ready for ggplot2::geom_segment()
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' env_arrows <- create_biplot_arrows(rda_result$env_scores, arrow_scale = 0.8)
#' }
create_biplot_arrows <- function(env_scores, arrow_scale = 1, min_length = 0) {
  if (is.matrix(env_scores)) {
    env_df <- as.data.frame(env_scores)
  } else {
    env_df <- env_scores
  }

  # Add variable names
  if (is.null(rownames(env_df))) {
    env_df$variable <- paste0("Var", seq_len(nrow(env_df)))
  } else {
    env_df$variable <- rownames(env_df)
  }

  # Scale arrows
  env_df[, seq_len(ncol(env_df) - 1)] <- env_df[, seq_len(ncol(env_df) - 1)] * arrow_scale

  # Calculate arrow length
  env_df$arrow_length <- sqrt(rowSums(env_df[, seq_len(ncol(env_df) - 1)]^2))

  # Filter by minimum length
  env_df <- env_df[env_df$arrow_length >= min_length, ]

  env_df
}

#' Validate Ordination Analysis Inputs
#'
#' Checks that input data and parameters are valid for ordination analysis.
#'
#' @param data Numeric matrix or data frame (samples × features)
#' @param sample_data Data frame with sample metadata
#' @param sample_id_col Column name in sample_data matching data row names (default: "sample")
#'
#' @return Invisibly returns TRUE if validation passes, stops with error otherwise
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' validate_ordination_inputs(expr_data, sample_metadata)
#' }
validate_ordination_inputs <- function(data, sample_data, sample_id_col = "sample") {
  # Check data
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(err_invalid_input("data", "matrix or data frame"))
  }

  if (nrow(data) == 0 || ncol(data) == 0) {
    stop(err_invalid_input("data", "non-empty matrix with samples and features"))
  }

  # Check for numeric data
  numeric_cols <- sapply(data, is.numeric)
  if (!all(numeric_cols)) {
    stop(err_non_numeric("data columns"))
  }

  # Check sample data
  if (!is.data.frame(sample_data)) {
    stop(err_invalid_input("sample_data", "data frame"))
  }

  if (nrow(sample_data) != nrow(data)) {
    stop(err_dimension_mismatch("sample_data", "data", nrow(sample_data), nrow(data)))
  }

  # Check sample ID column
  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(sample_data)) {
      stop("Sample ID column '", sample_id_col, "' not found in sample_data")
    }
  }

  invisible(TRUE)
}

#' Extract Component Names for Axis Labels
#'
#' Creates formatted axis labels with variance explained percentages.
#'
#' @param eigenvalues Data frame from \code{extract_eigenvalues()}
#' @param component Character: component name (e.g., "PC1", "PCo1")
#'
#' @return Formatted string like "PC1 (45.23%)"
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' eigs <- extract_eigenvalues(pca_result$eig, "pca")
#' x_label <- extract_component_label(eigs, "PC1")
#' }
extract_component_label <- function(eigenvalues, component) {
  row_idx <- which(eigenvalues$component == component)

  if (length(row_idx) == 0) {
    return(component)
  }

  var_pct <- eigenvalues$variance_percent[row_idx]
  paste0(component, " (", var_pct, "%)")
}

#' Prepare Eigenvalue Labels for Multi-Plot Display
#'
#' Formats eigenvalues for use in faceted or multi-component visualizations.
#'
#' @param eigenvalues Data frame from \code{extract_eigenvalues()}
#' @param components Character vector of components to label (e.g., c("PC1", "PC2"))
#'
#' @return Character vector of formatted labels
#'
#' @keywords internal
prepare_component_labels <- function(eigenvalues, components) {
  labels <- sapply(components, function(comp) {
    extract_component_label(eigenvalues, comp)
  }, USE.NAMES = FALSE)

  labels
}
