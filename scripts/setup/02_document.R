#' Generate Package Documentation
#'
#' Generate roxygen2 documentation for all functions in the package.
#'
#' @param check Logical. If TRUE, runs devtools::check() after documentation.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' generate_documentation(check = FALSE)
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

generate_documentation <- function(check = FALSE) {

  message("Generating package documentation...")

  # Check for required packages
  if (!requireNamespace("roxygen2", quietly = TRUE)) {
    stop("Package 'roxygen2' is required.", call. = FALSE)
  }

  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required.", call. = FALSE)
  }

  # Remove old NAMESPACE to ensure clean regeneration
  if (file.exists("NAMESPACE")) {
    message("  Removing old NAMESPACE")
    file.remove("NAMESPACE")
  }

  # Generate documentation
  message("  Running roxygen2::roxygenise()")
  tryCatch({
    roxygen2::roxygenise()
  }, error = function(e) {
    stop("Failed to generate documentation: ", conditionMessage(e), call. = FALSE)
  })

  message("✓ Documentation generated")
  message("  - NAMESPACE updated")
  message("  - man/ directory updated")

  # Optionally run check
  if (check) {
    message("\nRunning devtools::check()...")
    tryCatch({
      devtools::check()
    }, error = function(e) {
      warning("Check failed: ", conditionMessage(e))
    })
  }

  message("\n✓ Documentation generation complete!\n")
  invisible(TRUE)
}


#' Build and Install Package
#'
#' Build the package source tarball and install it locally.
#'
#' @param build Logical. If TRUE, builds the tarball. Default is TRUE.
#' @param install Logical. If TRUE, installs the package. Default is TRUE.
#' @param force Logical. If TRUE, forces installation. Default is TRUE.
#' @param check Logical. If TRUE, runs devtools::check() on built tarball.
#'   Default is FALSE.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' build_install_package(
#'   build = TRUE,
#'   install = TRUE,
#'   check = FALSE
#' )
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

build_install_package <- function(build = TRUE,
                                    install = TRUE,
                                    force = TRUE,
                                    check = FALSE) {

  message("Building and installing package...")

  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required.", call. = FALSE)
  }

  # Regenerate documentation first
  message("\n1. Generating documentation...")
  generate_documentation(check = FALSE)

  # Build package
  if (build) {
    message("\n2. Building package...")
    tryCatch({
      devtools::build()
      message("  ✓ Package built successfully")
    }, error = function(e) {
      stop("Failed to build package: ", conditionMessage(e), call. = FALSE)
    })
  }

  # Install package
  if (install) {
    message("\n3. Installing package...")
    tryCatch({
      devtools::install_local(force = force)
      message("  ✓ Package installed successfully")
    }, error = function(e) {
      stop("Failed to install package: ", conditionMessage(e), call. = FALSE)
    })
  }

  # Check package
  if (check) {
    message("\n4. Checking package...")

    # Find built tarball
    tarballs <- list.files(pattern = "\\.tar\\.gz$")
    if (length(tarballs) > 0) {
      tarball <- tarballs[1]
      tryCatch({
        devtools::check_built(tarball)
        message("  ✓ Package checked successfully")
      }, error = function(e) {
        warning("Check failed: ", conditionMessage(e))
      })
    } else {
      warning("No tarball found to check")
    }
  }

  message("\n✓ Build and install complete!\n")
  invisible(TRUE)
}


#' Update Package Version
#'
#' Update the version number in DESCRIPTION file.
#'
#' @param type Character. Type of version increment. Options: "major",
#'   "minor", "patch", or "dev". Default is "dev".
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Increment patch version (bug fix)
#' update_version("patch")
#'
#' # Increment minor version (new feature)
#' update_version("minor")
#'
#' # Increment major version (breaking change)
#' update_version("major")
#'
#' # Set to development version
#' update_version("dev")
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

update_version <- function(type = "dev") {

  if (!requireNamespace("usethis", quietly = TRUE)) {
    stop("Package 'usethis' is required.", call. = FALSE)
  }

  valid_types <- c("major", "minor", "patch", "dev")

  if (!type %in% valid_types) {
    stop("Invalid version type. Must be one of: ",
         paste(valid_types, collapse = ", "),
         call. = FALSE)
  }

  message("Updating package version...")

  tryCatch({
    switch(type,
           "major" = usethis::use_version("major"),
           "minor" = usethis::use_version("minor"),
           "patch" = usethis::use_version("patch"),
           "dev" = usethis::use_version("dev")
    )
    message("✓ Version updated to ", type, " version")
  }, error = function(e) {
    stop("Failed to update version: ", conditionMessage(e), call. = FALSE)
  })

  invisible(TRUE)
}


#' Run Package Tests
#'
#' Run all package tests using testthat.
#'
#' @param filter Character. If provided, only runs tests matching this pattern.
#' @param reporter Character. Reporter type. Default is "summary".
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Run all tests
#' run_package_tests()
#'
#' # Run specific tests
#' run_package_tests(filter = "plot_multi_volcano")
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

run_package_tests <- function(filter = NULL, reporter = "summary") {

  message("Running package tests...")

  if (!requireNamespace("testthat", quietly = TRUE)) {
    stop("Package 'testthat' is required.", call. = FALSE)
  }

  tryCatch({
    if (is.null(filter)) {
      testthat::test(reporter = reporter)
    } else {
      testthat::test_file(filter, reporter = reporter)
    }
    message("\n✓ Tests complete!")
  }, error = function(e) {
    warning("Tests failed: ", conditionMessage(e))
  })

  invisible(TRUE)
}
