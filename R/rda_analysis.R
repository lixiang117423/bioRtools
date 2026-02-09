#' Perform Redundancy Analysis (RDA) with environmental variables
#'
#' @description
#' Perform Redundancy Analysis (RDA) to explore the relationship between
#' community composition and environmental variables. RDA is a constrained
#' ordination method that identifies linear combinations of environmental
#' variables that best explain the variation in species composition.
#'
#' @param data A data frame where rows represent samples and columns represent
#'   species/OTUs/taxa. Should contain abundance or count data.
#' @param physicochemical A data frame containing environmental/physicochemical
#'   variables where rows represent samples and columns represent environmental
#'   variables. Sample names should match those in \code{data}.
#' @param sample A data frame containing sample information where the first
#'   column is 'sample' and subsequent columns contain grouping variables.
#' @param method_decostand Character string specifying the standardization
#'   method for \code{vegan::decostand()}. Default is "hellinger". Other
#'   options include "normalize", "standardize", "pa", "chi.square", etc.
#' @param scale Logical indicating whether to scale environmental variables
#'   in \code{vegan::rda()}. Default is FALSE.
#' @param color Character string specifying the column name in \code{sample}
#'   used for coloring points in the plot. Default is "group".
#' @param point_size Numeric value specifying the size of sample points.
#'   Default is 3.
#' @param arrow_size Numeric value specifying the size/length of environmental
#'   variable arrows. Default is 0.8.
#' @param arrow_color Character string specifying the color of environmental
#'   variable arrows. Default is "#222222".
#' @param arrow_alpha Numeric value between 0 and 1 for arrow transparency.
#'   Default is 0.8.
#' @param label_repel Logical indicating whether to use repelling labels for
#'   environmental variables to avoid overlap. Default is TRUE.
#' @param permutations Integer specifying the number of permutations for
#'   environmental fitting significance test. Default is 999.
#'
#' @return A list containing five components:
#' \describe{
#'   \item{result.rda}{The complete RDA result object from \code{vegan::rda()},
#'     which can be used for further analysis.}
#'   \item{plot.rda}{A ggplot object showing the RDA biplot with sample points
#'     and environmental variable arrows.}
#'   \item{sample.scores}{A data frame containing sample coordinates on RDA axes
#'     merged with sample information.}
#'   \item{environmental.scores}{A data frame containing environmental variable
#'     coordinates for plotting arrows.}
#'   \item{environmental.fit}{A data frame containing R² and p-values for
#'     environmental variables from permutation tests.}
#' }
#'
#' @details
#' RDA is a multivariate technique that:
#' \itemize{
#'   \item Combines regression and PCA to explore species-environment relationships
#'   \item Constrains ordination axes to be linear combinations of environmental variables
#'   \item Shows both constrained (explained by environment) and unconstrained variation
#'   \item Provides significance testing for environmental variables
#' }
#'
#' The function automatically:
#' \itemize{
#'   \item Standardizes species data using the specified method
#'   \item Performs RDA with all environmental variables
#'   \item Tests environmental variable significance using permutation tests
#'   \item Creates a publication-ready biplot
#' }
#'
#' @note
#' \itemize{
#'   \item Data should be pre-filtered to remove rare species if needed
#'   \item Environmental variables should be checked for collinearity
#'   \item Hellinger transformation is recommended for abundance data
#'   \item Consider the number of environmental variables relative to sample size
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example data
#' data(df.rda.otu)
#' data(df.rda.chem)
#'
#' # Prepare sample information
#' df.rda.chem %>%
#'   rownames() %>%
#'   as.data.frame() %>%
#'   magrittr::set_names("sample") %>%
#'   dplyr::mutate(
#'     group = stringr::str_split(sample, "_") %>% sapply("[", 1)
#'   ) -> df.sample
#'
#' # Basic RDA analysis
#' rda_result <- rda_analysis(
#'   data = df.rda.otu,
#'   physicochemical = df.rda.chem,
#'   sample = df.sample
#' )
#'
#' # View the plot
#' rda_result$plot.rda
#'
#' # Check environmental variable significance
#' rda_result$environmental.fit
#'
#' # Customized RDA with different parameters
#' rda_custom <- rda_analysis(
#'   data = df.rda.otu,
#'   physicochemical = df.rda.chem,
#'   sample = df.sample,
#'   method_decostand = "standardize",
#'   scale = TRUE,
#'   color = "group",
#'   point_size = 4
#' )
#'
#' # Access RDA summary statistics
#' summary(rda_result$result.rda)
#'
#' # Check variance explained
#' rda_result$result.rda$CCA$eig / rda_result$result.rda$tot.chi
#'
rda_analysis <- function(data,
                         physicochemical,
                         sample,
                         method_decostand = "hellinger",
                         scale = FALSE,
                         color = "group",
                         point_size = 3,
                         arrow_size = 0.8,
                         arrow_color = "#222222",
                         arrow_alpha = 0.8,
                         label_repel = TRUE,
                         permutations = 999) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.data.frame(physicochemical)) {
    stop("'physicochemical' must be a data frame")
  }

  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }

  # Check if sample names match between data and physicochemical
  if (!identical(rownames(data), rownames(physicochemical))) {
    stop("Sample names (row names) must match between 'data' and 'physicochemical'")
  }

  # Check if color column exists in sample data
  if (!color %in% colnames(sample)) {
    stop("Color column '", color, "' not found in sample data")
  }

  # Validate method_decostand
  valid_methods <- c("total", "max", "freq", "normalize", "range", "standardize",
    "pa", "chi.square", "hellinger", "log", "rrank")
  if (!method_decostand %in% valid_methods) {
    stop("'method_decostand' must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Check for missing values
  if (any(is.na(data))) {
    warning("Missing values found in 'data'. Consider imputation or removal.")
  }

  if (any(is.na(physicochemical))) {
    warning("Missing values found in 'physicochemical'. Consider imputation or removal.")
  }

  # Step 1: Standardize species/OTU data
  tryCatch(
    {
      data_standardized <- vegan::decostand(data, method = method_decostand)
    },
    error = function(e) {
      stop("Data standardization failed: ", e$message)
    })

  # Step 2: Perform RDA
  tryCatch(
    {
      rda_result <- vegan::rda(data_standardized ~ ., physicochemical, scale = scale)
    },
    error = function(e) {
      stop("RDA analysis failed: ", e$message)
    })

  # Step 3: Extract environmental variable coordinates (biplot scores)
  tryCatch(
    {
      environmental_scores <- vegan::scores(rda_result, display = "bp")

      if (is.null(environmental_scores)) {
        stop("Could not extract environmental variable scores")
      }

      environmental_coords <- environmental_scores %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "variable") %>%
        dplyr::mutate(
          # Scale arrows for better visualization
          RDA1_scaled = RDA1 * arrow_size,
          RDA2_scaled = RDA2 * arrow_size
        )
    },
    error = function(e) {
      stop("Failed to extract environmental scores: ", e$message)
    })

  # Step 4: Extract sample coordinates
  tryCatch(
    {
      sample_scores <- vegan::scores(rda_result, display = "sites") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample")

      # Merge with sample information
      sample_coords <- sample_scores %>%
        dplyr::left_join(sample, by = "sample")

      # Check if merge was successful
      if (nrow(sample_coords) != nrow(sample_scores)) {
        warning("Not all samples could be matched with sample information")
      }
    },
    error = function(e) {
      stop("Failed to extract sample scores: ", e$message)
    })

  # Step 5: Calculate eigenvalues and variance explained
  eigenvalues <- rda_result$CCA$eig
  total_variance <- rda_result$tot.chi
  variance_explained <- eigenvalues / total_variance * 100

  # Create axis labels with variance explained
  x_label <- paste0("RDA1 (", round(variance_explained[1], 1), "%)")
  y_label <- paste0("RDA2 (", round(variance_explained[2], 1), "%)")

  # Step 6: Create the RDA biplot
  rda_plot <- sample_coords %>%
    ggplot2::ggplot(ggplot2::aes(x = RDA1, y = RDA2)) +
    # Add reference lines
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    # Add sample points
    ggplot2::geom_point(
      ggplot2::aes(color = !!rlang::sym(color)),
      size = point_size,
      alpha = 0.8
    ) +
    # Add environmental variable arrows
    ggplot2::geom_segment(
      data = environmental_coords,
      ggplot2::aes(
        x = 0, y = 0,
        xend = RDA1_scaled, yend = RDA2_scaled
      ),
      arrow = ggplot2::arrow(angle = 20, length = ggplot2::unit(0.3, "cm")),
      color = arrow_color,
      alpha = arrow_alpha,
      linewidth = 0.8
    )

  # Add labels with or without repelling
  if (label_repel) {
    rda_plot <- rda_plot +
      ggrepel::geom_label_repel(
        data = environmental_coords,
        ggplot2::aes(x = RDA1_scaled, y = RDA2_scaled, label = variable),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50"
      )
  } else {
    rda_plot <- rda_plot +
      ggplot2::geom_text(
        data = environmental_coords,
        ggplot2::aes(x = RDA1_scaled, y = RDA2_scaled, label = variable),
        size = 3,
        vjust = -0.5
      )
  }

  # Apply styling and labels
  rda_plot <- rda_plot +
    ggsci::scale_color_d3() +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = "Redundancy Analysis (RDA)",
      subtitle = paste0("Standardization: ", method_decostand,
        " | Environmental variables: ", ncol(physicochemical))
    ) +
    bioRtools::theme_bio() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11)
    )

  # Step 7: Perform environmental fitting with permutation tests
  tryCatch(
    {
      env_fit <- vegan::envfit(rda_result, physicochemical, permutations = permutations)

      environmental_significance <- data.frame(
        variable = rownames(env_fit$vectors$arrows),
        r_squared = env_fit$vectors$r,
        p_value = env_fit$vectors$pvals,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::arrange(p_value) %>%
        dplyr::mutate(
          significance = dplyr::case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            p_value < 0.1 ~ ".",
            TRUE ~ ""
          ),
          r_squared_percent = round(r_squared * 100, 2)
        )
    },
    error = function(e) {
      warning("Environmental fitting failed: ", e$message)
      environmental_significance <- data.frame()
    })

  # Return comprehensive results following project conventions
  return(list(
    result.rda = rda_result,
    plot.rda = rda_plot,
    sample.scores = sample_coords,
    environmental.scores = environmental_coords,
    environmental.fit = environmental_significance
  ))
}
