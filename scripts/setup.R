#' BioRtools Package Development Script
#'
#' This script provides a complete toolkit for R package development.
#' It automates common development tasks such as:
#' - Package configuration (license, GitHub links)
#' - Documentation generation
#' - Code formatting
#' - Version updates
#' - Building and installing
#' - Running tests
#'
#' @section Usage:
#'
#' Source this script in an R session:
#'
#' \pre{`source("scripts/setup.R")`}
#'
#' Then call the individual functions as needed:
#'
#' \pre{`# Setup basic configuration`}
#' \pre{`setup_basic_config(github_user = "yourname", github_repo = "bioRtools")`}
#'
#' \pre{`# Add package imports`}
#' \pre{`add_package_imports()}`}
#'
#' \pre{`# Format all code`}
#' \pre{`format_package_code()`}
#'
#' \pre{`# Generate documentation`}
#' \pre{`generate_documentation()`}
#'
#' \pre{`# Build and install}
#' \pre{`build_install_package()}`}
#'
#' \pre{`# Run tests`}
#' \pre{`run_package_tests()}`}
#'
#' Or use the master setup function:
#'
#' \pre{`# Complete setup (does everything)}
#' \pre{`setup_package()}`
#'
#' @section Quick Start:
#'
#' For first-time setup:
#'
#' \pre{`# 1. Load this script`}
#' \pre{`source("scripts/setup.R")`}
#'
#' \pre{`# 2. Run complete setup}
#' \pre{`setup_package(github_user = "yourname", github_repo = "bioRtools")`}
#'
#' For routine development:
#'
#' \pre{`# Just format code and update docs}
#' \pre{`format_and_document()}`}
#'
#' \pre{`# Just build and install}
#' \pre{`build_install_package(check = FALSE)}`}
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}

# Required packages
required_packages <- c("usethis", "styler", "roxygen2", "devtools")

# Check for required packages
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  message("Installing required packages...")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
  message("✓ Required packages installed")
}

# Load helper scripts
# Get the directory where this script is located
get_script_dir <- function() {
  # Try multiple methods to get script directory
  script_path <- NULL

  # Method 1: thisfile() (R 4.0+)
  if (exists("thisfile", mode = "function")) {
    script_path <- tryCatch(thisfile(), error = function(e) NULL)
  }

  # Method 2: sys.frame(1)$filename
  if (is.null(script_path)) {
    script_path <- tryCatch(sys.frame(1)$filename, error = function(e) NULL)
  }

  # Method 3: Use current working directory with relative path
  if (is.null(script_path) || !file.exists(script_path)) {
    # Assume script is in scripts/setup.R relative to package root
    script_path <- file.path(getwd(), "scripts", "setup.R")
  }

  if (!is.null(script_path) && file.exists(script_path)) {
    return(dirname(script_path))
  }

  # Fallback: return NULL and handle it later
  return(NULL)
}

script_dir <- get_script_dir()

if (is.null(script_dir)) {
  # If we can't determine script directory, use relative path
  script_dir <- file.path(getwd(), "scripts")
}

helper_dir <- file.path(script_dir, "setup")

# Check if helper directory exists
if (!dir.exists(helper_dir)) {
  # Try alternative paths
  possible_paths <- c(
    file.path(getwd(), "scripts", "setup"),
    file.path(dirname(getwd()), "scripts", "setup"),
    "scripts/setup",
    "../scripts/setup"
  )

  for (path in possible_paths) {
    if (dir.exists(path)) {
      helper_dir <- path
      break
    }
  }
}

# Load helper scripts if directory exists
if (dir.exists(helper_dir)) {
  source_paths <- list.files(helper_dir, pattern = "\\.R$", full.names = TRUE)

  for (script in source_paths) {
    tryCatch(source(script), error = function(e) {
      warning(paste("Failed to load", script, ":", e$message))
    })
  }
  message("✓ Loaded ", length(source_paths), " helper script(s)\n")
} else {
  warning("Could not find helper scripts directory. Some functions may not be available.")
}


#' Complete Package Setup
#'
#' This function performs all package setup steps in sequence:
#' 1. Basic configuration
#' 2. Add imports
#' 3. Format code
#' 4. Generate documentation
#' 5. Build and install
#'
#' @param github_user Character. GitHub username.
#' @param github_repo Character. GitHub repository name.
#' @param license Character. License type. Default is "mit".
#' @param build Logical. Whether to build the package. Default is TRUE.
#' @param install Logical. Whether to install the package. Default is TRUE.
#' @param format_code Logical. Whether to format code. Default is TRUE.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Complete setup
#' setup_package(
#'   github_user = "yourname",
#'   github_repo = "bioRtools",
#'   license = "mit"
#' )
#'
#' # Setup without building/installing
#' setup_package(
#'   github_user = "yourname",
#'   github_repo = "bioRtools",
#'   build = FALSE,
#'   install = FALSE
#' )
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

setup_package <- function(github_user = NULL,
                          github_repo = NULL,
                          license = "mit",
                          build = TRUE,
                          install = TRUE,
                          format_code = TRUE) {

  message("========================================")
  message("   BioRtools Package Setup")
  message("========================================\n")

  # 1. Basic configuration
  message("Step 1/4: Basic configuration")
  message("----------------------------------------")
  setup_basic_config(github_user, github_repo, license)

  # 2. Add imports
  message("\nStep 2/4: Adding package imports")
  message("----------------------------------------")
  add_package_imports()

  # 3. Format code
  if (format_code) {
    message("\nStep 3/4: Formatting code")
    message("----------------------------------------")
    format_package_code()
  }

  # 4. Documentation and build
  message("\nStep 4/4: Documentation, build, and install")
  message("----------------------------------------")
  generate_documentation()

  if (build || install) {
    build_install_package(build = build, install = install)
  }

  message("\n========================================")
  message("   Setup Complete!")
  message("========================================\n")

  message("Next steps:")
  message("  1. Review the changes")
  message("  2. Test the package: run_package_tests()")
  message("  3. Commit to git")

  invisible(TRUE)
}


#' Format and Document
#'
#' Quick function to format code and regenerate documentation.
#' Useful for routine development workflow.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' format_and_document()
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

format_and_document <- function() {

  message("Formatting code and updating documentation...")

  format_package_code()
  generate_documentation(check = FALSE)

  message("✓ Format and document complete!\n")
  invisible(TRUE)
}


#' Quick Install
#'
#' Quick function to regenerate documentation and install the package.
#' Useful for testing changes during development.
#'
#' @param force Logical. Force reinstallation. Default is TRUE.
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' quick_install(force = TRUE)
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export

quick_install <- function(force = TRUE) {

  message("Quick install...")

  generate_documentation(check = FALSE)
  pak::pak("local::.", upgrade = FALSE)

  message("✓ Package installed!\n")
  invisible(TRUE)
}


# Print usage information when script is sourced
message("\n")
message("BioRtools Development Toolkit")
message("============================")
message("\nAvailable functions:")
message("  - setup_package()          : Complete package setup")
message("  - setup_basic_config()      : Configure license, GitHub")
message("  - add_package_imports()     : Add @importFrom declarations")
message("  - format_package_code()     : Format all R code with styler")
message("  - generate_documentation()  : Generate roxygen2 docs")
message("  - build_install_package()    : Build and install package")
message("  - update_version()           : Update package version")
message("  - run_package_tests()        : Run package tests")
message("  - format_and_document()      : Quick format + document")
message("  - quick_install()            : Quick install after changes")
message("\nExample usage:")
message("  setup_package(github_user = 'yourname', github_repo = 'bioRtools')")
message("  format_and_document()")
message("  quick_install()")
message("\nType ?function_name for help\n")
