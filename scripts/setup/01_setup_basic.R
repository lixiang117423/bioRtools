#' Setup Basic Package Configuration
#'
#' This script configures the basic package settings including license,
#' GitHub links, and build/git ignore files.
#'
#' @param github_user Character. GitHub username. Default is NULL.
#' @param github_repo Character. GitHub repository name. Default is NULL.
#' @param license Character. License type. Options: "mit", "apache", "gpl".
#'   Default is "mit".
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' setup_basic_config(
#'   github_user = "username",
#'   github_repo = "bioRtools",
#'   license = "mit"
#' )
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

setup_basic_config <- function(github_user = NULL,
                                github_repo = NULL,
                                license = "mit") {

  # Check if usethis is available
  if (!requireNamespace("usethis", quietly = TRUE)) {
    stop("Package 'usethis' is required. Please install it with: install.packages('usethis')",
         call. = FALSE)
  }

  message("Setting up basic package configuration...")

  # Add license
  message("✓ Adding ", toupper(license), " license")
  tryCatch({
    switch(license,
           "mit" = usethis::use_mit_license(),
           "apache" = usethis::use_apache_license(),
           "gpl" = usethis::use_gpl3_license(),
           stop("Unknown license type: ", license)
    )
  }, error = function(e) {
    warning("Failed to add license: ", conditionMessage(e))
  })

  # Add GitHub links if user and repo provided
  if (!is.null(github_user) && !is.null(github_repo)) {
    message("✓ Adding GitHub links")
    tryCatch({
      usethis::use_github_links()
    }, error = function(e) {
      warning("Failed to add GitHub links: ", conditionMessage(e))
    })
  }

  # Add build ignores
  message("✓ Configuring .Rbuildignore")
  build_ignores <- c(
    "README.md",
    "deve.log.R",
    "test_function.R",
    "*.Rproj",
    "raw.data/",
    ".Rprofile",
    ".Rhistory"
  )

  for (item in build_ignores) {
    tryCatch({
      usethis::use_build_ignore(item)
    }, error = function(e) {
      # Ignore if already exists
    })
  }

  # Add git ignores
  message("✓ Configuring .gitignore")
  git_ignores <- c(
    "raw.data/",
    "*.Rproj.user",
    ".Rproj.user/",
    ".Rhistory",
    "test_*.png",
    "test_*.pdf",
    "*.tar.gz",
    "*.Rcheck/"
  )

  for (item in git_ignores) {
    tryCatch({
      usethis::use_git_ignore(item)
    }, error = function(e) {
      # Ignore if already exists
    })
  }

  message("\n✓ Basic configuration complete!\n")
  invisible(TRUE)
}


#' Add Package Imports
#'
#' Add @importFrom declarations to NAMESPACE for commonly used functions.
#' This function maintains a record of external dependencies.
#'
#' @param imports Character vector. Import declarations in the format
#'   "package::function". If NULL, adds a default set of imports.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Add default imports
#' add_package_imports()
#'
#' # Add custom imports
#' add_package_imports(c(
#'   "dplyr::mutate",
#'   "ggplot2::aes"
#' ))
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

add_package_imports <- function(imports = NULL) {

  if (!requireNamespace("usethis", quietly = TRUE)) {
    stop("Package 'usethis' is required.", call. = FALSE)
  }

  # Default imports commonly used in bioRtools
  default_imports <- c(
    "dplyr" = c("%>%", "arrange", "mutate", "select", "filter",
               "left_join", "rename", "count", "n", "slice_max",
               "lag", "inner_join"),
    "ggplot2" = c("ggplot", "aes", "geom_point", "geom_hline",
                   "geom_vline", "geom_smooth", "labs", "theme",
                   "element_line", "element_rect", "element_text",
                   "element_blank", "unit", "rel", "guides",
                   "geom_segment", "position_fill"),
    "tidyr" = c("pivot_longer", "nest"),
    "tibble" = c("rownames_to_column"),
    "purrr" = c("map", "map_dbl"),
    "stringr" = c("str_replace", "str_split"),
    "rlang" = c("sym", "!!"),
    "scales" = c("percent"),
    "stats" = c("anova", "quantile", "IQR", "reorder"),
    "ggrepel" = c("geom_text_repel", "geom_label_repel"),
    "ggthemes" = c("theme_foundation")
  )

  if (is.null(imports)) {
    message("Adding default package imports...")
    imports_list <- default_imports
  } else {
    message("Adding custom package imports...")
    # Parse custom imports
    imports_list <- list()
    for (imp in imports) {
      parts <- strsplit(imp, "::")[[1]]
      if (length(parts) == 2) {
        pkg <- parts[1]
        func <- parts[2]
        if (is.null(imports_list[[pkg]])) {
          imports_list[[pkg]] <- character(0)
        }
        imports_list[[pkg]] <- c(imports_list[[pkg]], func)
      }
    }
  }

  # Add imports using usethis
  for (pkg in names(imports_list)) {
    funcs <- imports_list[[pkg]]
    for (func in funcs) {
      tryCatch({
        usethis::use_import_from(pkg, func)
      }, error = function(e) {
        message("  Note: ", pkg, "::", func, " - ", conditionMessage(e))
      })
    }
  }

  message("✓ Package imports complete!\n")
  invisible(TRUE)
}


#' Format All R Code
#'
#' Format all R files in the package using styler.
#'
#' @param strict Logical. If TRUE, uses strict formatting. Default is FALSE.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' format_package_code(strict = FALSE)
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

format_package_code <- function(strict = FALSE) {

  if (!requireNamespace("styler", quietly = TRUE)) {
    stop("Package 'styler' is required. Please install it with: install.packages('styler')",
         call. = FALSE)
  }

  message("Formatting R code with styler...")

  # Get all R files
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)

  if (length(r_files) == 0) {
    message("No R files found in R/ directory")
    return(invisible(TRUE))
  }

  message("Found ", length(r_files), " R files")

  # Format each file
  for (i in seq_along(r_files)) {
    file <- r_files[i]
    message(sprintf("  [%d/%d] %s", i, length(r_files), basename(file)))

    tryCatch({
      styler::style_file(file, strict = strict)
    }, error = function(e) {
      warning("Failed to format ", basename(file), ": ", conditionMessage(e))
    })
  }

  message("\n✓ Code formatting complete!\n")
  invisible(TRUE)
}
