#' Perform Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)
#'
#' @description
#' Perform sparse Partial Least Squares Discriminant Analysis (sPLS-DA) for
#' supervised dimension reduction and classification of high-dimensional data.
#' sPLS-DA is particularly useful for metabolomics, genomics, and other omics
#' data where variable selection and class discrimination are important.
#'
#' @param data Numeric matrix or data frame with observations as rows and
#'   variables (features) as columns. Missing values (NAs) are allowed.
#' @param group A factor or character vector indicating the class membership
#'   for each observation. Must have the same length as the number of rows in data.
#' @param ncomp Positive integer specifying the number of components to include
#'   in the model. Default is 3.
#' @param keepX Numeric vector indicating the number of variables to keep in
#'   each component. If NULL (default), all variables are kept (non-sparse).
#' @param scale Logical indicating whether to scale the data. Default is TRUE.
#' @param center Logical indicating whether to center the data. Default is TRUE.
#' @param max_iter Maximum number of iterations for the algorithm. Default is 500.
#' @param tol Convergence tolerance for the algorithm. Default is 1e-06.
#' @param validation Character string specifying the validation method for
#'   parameter tuning. Options: "Mfold", "loo". Default is NULL (no validation).
#' @param folds Number of folds for cross-validation if validation is used.
#'   Default is 10.
#' @param nrepeat Number of repeats for cross-validation. Default is 5.
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A list containing six components:
#' \describe{
#'   \item{result.splsda}{The complete sPLS-DA result object from
#'     \code{mixOmics::splsda()}, which can be used for further analysis.}
#'   \item{sample.scores}{A data frame containing sample coordinates on
#'     sPLS-DA components with sample names.}
#'   \item{variance.explained}{A data frame containing the variance explained
#'     by each component with formatted percentages and axis labels.}
#'   \item{variable.loadings}{A data frame containing variable loadings for
#'     each component (if available).}
#'   \item{classification.performance}{Cross-validation performance metrics
#'     (if validation was performed).}
#'   \item{model.parameters}{A summary of model parameters and settings used.}
#' }
#'
#' @details
#' sPLS-DA combines the benefits of PLS-DA and variable selection:
#' \itemize{
#'   \item Supervised dimension reduction for classification
#'   \item Automatic variable selection through sparsity
#'   \item Handles high-dimensional data with small sample sizes
#'   \item Provides interpretable results with selected biomarkers
#' }
#'
#' The function automatically:
#' \itemize{
#'   \item Validates input data and group matching
#'   \item Performs sPLS-DA with specified parameters
#'   \item Extracts and formats key results
#'   \item Calculates variance explained by each component
#'   \item Provides ready-to-use plotting coordinates
#' }
#'
#' @note
#' \itemize{
#'   \item Groups should be balanced when possible for optimal performance
#'   \item Consider parameter tuning using \code{tune.splsda()} for optimal results
#'   \item Missing values are handled by the mixOmics package internally
#'   \item For large datasets, consider reducing ncomp to improve computational efficiency
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example metabolomics data
#' data("df.splsda.meta")
#' data("df.splsda.sample")
#'
#' # Basic sPLS-DA analysis
#' splsda_result <- spls_analysis(
#'   data = df.splsda.meta,
#'   group = as.factor(df.splsda.sample$day)
#' )
#'
#' # View sample scores
#' head(splsda_result$sample.scores)
#'
#' # Check variance explained
#' splsda_result$variance.explained
#'
#' # View model parameters
#' splsda_result$model.parameters
#'
#' # Customized sPLS-DA with variable selection
#' splsda_sparse <- spls_analysis(
#'   data = df.splsda.meta,
#'   group = as.factor(df.splsda.sample$day),
#'   ncomp = 2,
#'   keepX = c(50, 30),  # Select 50 variables for comp1, 30 for comp2
#'   scale = TRUE
#' )
#'
#' # sPLS-DA with cross-validation
#' splsda_cv <- spls_analysis(
#'   data = df.splsda.meta,
#'   group = as.factor(df.splsda.sample$day),
#'   validation = "Mfold",
#'   folds = 5,
#'   nrepeat = 3
#' )
#'
#' # Access cross-validation results
#' splsda_cv$classification.performance
#'
#' # Create a scatter plot
#' library(ggplot2)
#' splsda_result$sample.scores %>%
#'   dplyr::mutate(group = as.factor(df.splsda.sample$day)) %>%
#'   ggplot(aes(x = comp1, y = comp2, color = group)) +
#'   geom_point(size = 3, alpha = 0.8) +
#'   labs(
#'     x = splsda_result$variance.explained$label[1],
#'     y = splsda_result$variance.explained$label[2],
#'     title = "sPLS-DA Score Plot"
#'   ) +
#'   theme_minimal()
#'
spls_analysis <- function(data,
                          group,
                          ncomp = 3,
                          keepX = NULL,
                          scale = TRUE,
                          center = TRUE,
                          max_iter = 500,
                          tol = 1e-06,
                          validation = NULL,
                          folds = 10,
                          nrepeat = 5,
                          verbose = TRUE) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame")
  }

  # Convert to matrix if data frame
  if (is.data.frame(data)) {
    if (!all(sapply(data, is.numeric))) {
      stop("All columns in 'data' must be numeric")
    }
    data <- as.matrix(data)
  }

  # Check dimensions
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("'data' cannot be empty")
  }

  # Validate group
  if (length(group) != nrow(data)) {
    stop("Length of 'group' must equal the number of rows in 'data'")
  }

  # Convert group to factor if not already
  if (!is.factor(group)) {
    group <- as.factor(group)
    if (verbose) {
      message("Converting 'group' to factor with levels: ", paste(levels(group), collapse = ", "))
    }
  }

  # Check for minimum samples per group
  group_counts <- table(group)
  if (any(group_counts < 2)) {
    warning("Some groups have fewer than 2 samples. This may affect model stability.")
  }

  # Validate parameters
  if (!is.numeric(ncomp) || ncomp < 1 || ncomp != round(ncomp)) {
    stop("'ncomp' must be a positive integer")
  }

  if (ncomp > min(nrow(data) - 1, ncol(data))) {
    warning("'ncomp' is larger than the maximum possible components. Adjusting automatically.")
    ncomp <- min(nrow(data) - 1, ncol(data))
  }

  # Validate keepX if provided
  if (!is.null(keepX)) {
    if (length(keepX) != ncomp) {
      stop("Length of 'keepX' must equal 'ncomp'")
    }
    if (any(keepX > ncol(data))) {
      stop("Values in 'keepX' cannot exceed the number of variables")
    }
  }

  # Validate validation method
  if (!is.null(validation)) {
    valid_methods <- c("Mfold", "loo")
    if (!validation %in% valid_methods) {
      stop("'validation' must be one of: ", paste(valid_methods, collapse = ", "))
    }
  }

  if (verbose) {
    message("Starting sPLS-DA analysis...")
    message("Data dimensions: ", nrow(data), " samples × ", ncol(data), " variables")
    message("Groups: ", paste(names(group_counts), " (n=", group_counts, ")", collapse = ", "))
    message("Components: ", ncomp)
    if (!is.null(keepX)) {
      message("Variable selection: ", paste(keepX, collapse = ", "), " variables per component")
    }
  }

  # Perform sPLS-DA
  tryCatch(
    {
      splsda_result <- mixOmics::splsda(
        X = data,
        Y = group,
        ncomp = ncomp,
        keepX = keepX,
        scale = scale,
        center = center,
        max.iter = max_iter,
        tol = tol
      )
    },
    error = function(e) {
      stop("sPLS-DA analysis failed: ", e$message)
    })

  # Extract sample scores
  sample_scores <- splsda_result$variates$X %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::rename_with(~ paste0("comp", seq_along(.x)), -sample)

  # Calculate variance explained and create labels
  variance_explained <- splsda_result$prop_expl_var$X %>%
    as.data.frame() %>%
    magrittr::set_names("proportion") %>%
    dplyr::mutate(
      component = paste0("comp", seq_len(nrow(.))),
      percentage = round(proportion * 100, 2),
      percentage_label = paste0(percentage, "%"),
      label = paste0("Component ", seq_len(nrow(.)), " (", percentage_label, ")")
    ) %>%
    dplyr::select(component, proportion, percentage, percentage_label, label)

  # Extract variable loadings if available
  variable_loadings <- NULL
  if (!is.null(splsda_result$loadings$X)) {
    variable_loadings <- splsda_result$loadings$X %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "variable") %>%
      dplyr::rename_with(~ paste0("comp", seq_along(.x)), -variable)
  }

  # Perform cross-validation if requested
  classification_performance <- NULL
  if (!is.null(validation)) {
    if (verbose) {
      message("Performing cross-validation...")
    }

    tryCatch(
      {
        cv_result <- mixOmics::perf(
          splsda_result,
          validation = validation,
          folds = folds,
          nrepeat = nrepeat,
          progressBar = FALSE
        )

        # Extract performance metrics
        classification_performance <- data.frame(
          component = paste0("comp", seq_len(ncomp)),
          overall_error = cv_result$error.rate$overall[, "max.dist"],
          balanced_error = cv_result$error.rate$BER[, "max.dist"],
          stringsAsFactors = FALSE
        )

        if (verbose) {
          message("Cross-validation completed. Overall error rates:")
          for (i in seq_len(nrow(classification_performance))) {
            message("  ", classification_performance$component[i], ": ",
              round(classification_performance$overall_error[i] * 100, 2), "%")
          }
        }
      },
      error = function(e) {
        warning("Cross-validation failed: ", e$message)
      })
  }

  # Create model parameters summary
  model_parameters <- data.frame(
    parameter = c("ncomp", "scale", "center", "max_iter", "tol", "n_samples",
      "n_variables", "n_groups", "sparse"),
    value = c(ncomp, scale, center, max_iter, tol, nrow(data), ncol(data),
      nlevels(group), !is.null(keepX)),
    stringsAsFactors = FALSE
  )

  if (!is.null(keepX)) {
    keepX_summary <- data.frame(
      parameter = paste0("keepX_comp", seq_along(keepX)),
      value = keepX,
      stringsAsFactors = FALSE
    )
    model_parameters <- rbind(model_parameters, keepX_summary)
  }

  if (verbose) {
    message("sPLS-DA analysis completed successfully!")
    message("Total variance explained by ", ncomp, " components: ",
      round(sum(variance_explained$percentage), 2), "%")
  }

  # Return comprehensive results following project conventions
  return(list(
    result.splsda = splsda_result,
    sample.scores = sample_scores,
    variance.explained = variance_explained,
    variable.loadings = variable_loadings,
    classification.performance = classification_performance,
    model.parameters = model_parameters
  ))
}
