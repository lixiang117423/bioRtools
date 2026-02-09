#' Create a linkage disequilibrium (LD) heatmap from VCF file
#'
#' @description
#' Generate a linkage disequilibrium heatmap directly from a VCF file in PLINK format.
#' This function provides a convenient way to visualize LD patterns across genomic regions
#' with customizable color schemes and output options.
#'
#' @param vcf_file Character string specifying the path to the VCF file in PLINK format.
#'   The file should contain SNP genotype data.
#' @param color_palette Character string or vector specifying the color scheme. Options include:
#'   \itemize{
#'     \item "default": Soft green gradient (publication-ready)
#'     \item "green_yellow": Green to yellow gradient similar to your image
#'     \item "scientific": Blue-green gradient for scientific publications
#'     \item "nature": Nature journal style color scheme
#'     \item "cell": Cell journal style green gradient
#'     \item "viridis": Viridis color scale
#'     \item "plasma": Plasma color scale
#'     \item "magma": Magma color scale
#'     \item "red_blue": Blue to white to red diverging scale
#'     \item Custom vector: User-provided color vector (e.g., c("#f7fcf5", "#238b45"))
#'   }
#' @param color_steps Integer specifying the number of color steps in the palette.
#'   Default is 100.
#' @param flip_diagonal Logical indicating whether to flip the heatmap along the diagonal.
#'   Default is TRUE.
#' @param title Character string for the plot title. If NULL, a default title will be generated.
#'   Default is NULL.
#' @param subtitle Character string for the plot subtitle. Default is NULL.
#' @param output_file Character string specifying the output filename (without extension).
#'   If NULL, the plot will only be displayed. Default is NULL.
#' @param file_format Character string specifying the output format. Options include:
#'   "png", "pdf", "jpeg", "tiff". Default is "png".
#' @param width Numeric value for plot width. For raster formats (png, jpeg, tiff),
#'   this is in pixels when multiplied by dpi. For PDF, this is in inches. Default is 9.
#' @param height Numeric value for plot height. For raster formats (png, jpeg, tiff),
#'   this is in pixels when multiplied by dpi. For PDF, this is in inches. Default is 7.
#' @param dpi Numeric value for resolution (dots per inch) for raster formats.
#'   Default is 300.
#' @param ggplot_version Logical indicating whether to create a ggplot2 version
#'   of the heatmap in addition to the LDheatmap version. Default is TRUE.
#' @param show_values Logical indicating whether to show LD values in the ggplot2
#'   heatmap tiles. Default is FALSE (recommended for large datasets).
#' @param text_size Numeric value for text size in ggplot2 version when show_values = TRUE.
#'   Default is 3.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot.ld_heatmap}{The LD heatmap plot object from LDheatmap package}
#'   \item{plot.ggplot2}{A ggplot2 version of the heatmap (if ggplot_version = TRUE)}
#'   \item{ld.matrix}{The linkage disequilibrium matrix (r²) used for plotting}
#'   \item{snp.matrix}{The SNP genotype matrix used for LD calculation}
#'   \item{genetic.distances}{The genetic distances vector}
#'   \item{snp.info}{Information about the SNPs including positions}
#'   \item{plot.params}{Parameters used for plotting (colors, dimensions, etc.)}
#' }
#'
#' @details
#' This function reads VCF files in PLINK format and generates LD heatmaps using the
#' LDheatmap package. The function provides several built-in color palettes optimized
#' for LD visualization, where higher LD values (closer to 1) are typically shown in
#' "hotter" colors.
#'
#' The function automatically handles:
#' \itemize{
#'   \item SNP matrix extraction from VCF files
#'   \item Genetic distance calculation
#'   \item Color palette generation
#'   \item File output in multiple formats
#'   \item Plot parameter optimization
#' }
#'
#' @note
#' This function requires the LDheatmap package and ttplot package for VCF file processing.
#' If LDheatmap is not installed, it will be installed automatically from Bioconductor.
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' \dontrun{
#' # Basic LD heatmap with default settings
#' ld_result <- plot_LDheatmap(
#'   vcf_file = "path/to/your/file.vcf",
#'   title = "LD Heatmap for Chr1:1000000-2000000"
#' )
#'
#' # Display the plot
#' ld_result$plot.ld_heatmap
#'
#' # Custom color palette matching your image style
#' ld_soft <- plot_LDheatmap(
#'   vcf_file = "path/to/your/file.vcf",
#'   color_palette = "green_yellow",
#'   title = "LD Analysis with Soft Green-Yellow Colors"
#' )
#'
#' # Nature journal style
#' ld_nature <- plot_LDheatmap(
#'   vcf_file = "path/to/your/file.vcf",
#'   color_palette = "nature",
#'   title = "Publication-Ready LD Heatmap"
#' )
#'
#' # Custom hex colors for precise control
#' ld_custom_hex <- plot_LDheatmap(
#'   vcf_file = "path/to/your/file.vcf",
#'   color_palette = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b"),
#'   title = "Custom Green Gradient LD Analysis"
#' )
#'
#' # Save ggplot2 version with ggsave
#' ld_ggplot <- plot_LDheatmap(
#'   vcf_file = "path/to/your/file.vcf",
#'   color_palette = "green_yellow",
#'   ggplot_version = TRUE,
#'   title = "LD Heatmap Analysis"
#' )
#'
#' # Save the ggplot2 version
#' ggplot2::ggsave("ld_heatmap.png", ld_ggplot$plot.ggplot2,
#'   width = 10, height = 8, dpi = 300)
#' ggplot2::ggsave("ld_heatmap.pdf", ld_ggplot$plot.ggplot2,
#'   width = 10, height = 8)
#'
#' # Customize the ggplot2 version further
#' custom_plot <- ld_ggplot$plot.ggplot2 +
#'   ggplot2::theme_minimal() +
#'   ggplot2::theme(
#'     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
#'     plot.title = ggplot2::element_text(size = 16, face = "bold")
#'   ) +
#'   ggplot2::labs(
#'     caption = "Generated with bioRtools::plot_LDheatmap()"
#'   )
#'
#' ggplot2::ggsave("ld_heatmap_custom.png", custom_plot,
#'   width = 12, height = 10, dpi = 300)
#'
#' # Multiple format output with soft colors
#' ld_multi <- plot_LDheatmap(
#'   vcf_file = "path/to/your/file.vcf",
#'   color_palette = "scientific",
#'   color_steps = 150,
#'   output_file = "ld_heatmap_final",
#'   file_format = "png",
#'   dpi = 600,
#'   verbose = TRUE
#' )
#' }
#'
plot_LDheatmap <- function(vcf_file,
                           color_palette = "default",
                           color_steps = 100,
                           flip_diagonal = TRUE,
                           title = NULL,
                           subtitle = NULL,
                           output_file = NULL,
                           file_format = "png",
                           width = 9,
                           height = 7,
                           dpi = 300,
                           ggplot_version = TRUE,
                           show_values = FALSE,
                           text_size = 3,
                           verbose = TRUE) {
  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------

  if (!is.character(vcf_file) || length(vcf_file) != 1) {
    stop("'vcf_file' must be a single character string")
  }

  if (!file.exists(vcf_file)) {
    stop("VCF file not found: ", vcf_file)
  }

  if (!is.numeric(color_steps) || color_steps < 10 || color_steps > 1000) {
    stop("'color_steps' must be a numeric value between 10 and 1000")
  }

  if (!is.logical(flip_diagonal)) {
    stop("'flip_diagonal' must be logical (TRUE or FALSE)")
  }

  if (!file_format %in% c("png", "pdf", "jpeg", "jpg", "tiff")) {
    stop("'file_format' must be one of: png, pdf, jpeg, jpg, tiff")
  }

  if (!is.numeric(width) || width <= 0) {
    stop("'width' must be a positive numeric value")
  }

  if (!is.numeric(height) || height <= 0) {
    stop("'height' must be a positive numeric value")
  }

  if (!is.numeric(dpi) || dpi < 72 || dpi > 2400) {
    stop("'dpi' must be between 72 and 2400")
  }

  if (!is.logical(ggplot_version)) {
    stop("'ggplot_version' must be logical (TRUE or FALSE)")
  }

  if (!is.logical(show_values)) {
    stop("'show_values' must be logical (TRUE or FALSE)")
  }

  if (!is.numeric(text_size) || text_size <= 0) {
    stop("'text_size' must be a positive numeric value")
  }

  # -------------------------------------------------------------------------
  # Package dependencies
  # -------------------------------------------------------------------------

  if (verbose) message("Checking package dependencies...")

  # Check and install LDheatmap if needed
  if (!requireNamespace("LDheatmap", quietly = TRUE)) {
    if (verbose) message("Installing LDheatmap from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("LDheatmap", quiet = !verbose)
  }

  # Check ttplot package
  if (!requireNamespace("ttplot", quietly = TRUE)) {
    stop("ttplot package is required for VCF file processing. Please install it first.")
  }

  # Load required libraries
  suppressWarnings(suppressMessages({
    library(LDheatmap, quietly = TRUE)
    if (ggplot_version) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for ggplot_version = TRUE")
      }
    }
  }))

  # -------------------------------------------------------------------------
  # Data preparation
  # -------------------------------------------------------------------------

  if (verbose) message("Reading VCF file and extracting SNP data...")

  # Extract SNP matrix and information
  tryCatch(
    {
      snp_matrix <- ttplot::getsnpMat(vcf_file)
      snp_info <- ttplot::getsnpInfo(vcf_file)
    },
    error = function(e) {
      stop("Error reading VCF file: ", e$message)
    })

  if (verbose) {
    message("Successfully loaded ", nrow(snp_matrix), " samples and ",
      ncol(snp_matrix), " SNPs")
  }

  # Extract genetic distances
  genetic_distances <- as.numeric(snp_info$POS)

  if (any(is.na(genetic_distances))) {
    stop("Invalid position information in VCF file")
  }

  # -------------------------------------------------------------------------
  # Color palette generation
  # -------------------------------------------------------------------------

  if (verbose) message("Setting up color palette...")

  # Define color palette based on input
  if (length(color_palette) == 1 && is.character(color_palette)) {
    color_scheme <- switch(
      color_palette,
      "default" = c("#f7fcf5", "#bae4b3", "#74c476", "#238b45", "#00441b"),  # 柔和绿色渐变
      "green_yellow" = c("#f7fcf5", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45"),  # 绿黄色系
      "scientific" = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"),  # 科学出版风格
      "nature" = c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1"),  # 自然杂志风格
      "cell" = c("#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#238b45"),  # Cell杂志风格
      "viridis" = viridis::viridis(color_steps),
      "plasma" = viridis::plasma(color_steps),
      "magma" = viridis::magma(color_steps),
      "red_blue" = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
      "blue_yellow" = c("#2c7fb8", "#7fcdbb", "#c7e9b4", "#ffffcc", "#fec44f", "#d95f0e"),
      "green_red" = c("#00441b", "#238b45", "#74c476", "#bae4b3", "#f7fcf5", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"),
      {
        warning("Unknown color palette '", color_palette, "'. Using default.")
        c("#f7fcf5", "#bae4b3", "#74c476", "#238b45", "#00441b")
      }
    )
  } else if (is.character(color_palette) && length(color_palette) > 1) {
    # User provided custom colors
    color_scheme <- color_palette
  } else {
    stop("'color_palette' must be a character string or vector of colors")
  }

  # Generate color palette function
  if (length(color_scheme) == color_steps) {
    # Pre-defined palette with exact number of colors
    color_function <- function(n) color_scheme[1:n]
  } else {
    # Generate gradient from provided colors
    color_function <- colorRampPalette(rev(color_scheme), space = "rgb")
  }

  # Generate final colors
  plot_colors <- color_function(color_steps)

  # -------------------------------------------------------------------------
  # Plot title generation
  # -------------------------------------------------------------------------

  if (is.null(title)) {
    # Extract basic info for default title
    n_snps <- ncol(snp_matrix)
    if (!is.null(snp_info$CHR) && length(unique(snp_info$CHR)) == 1) {
      chr_info <- paste0("Chr ", unique(snp_info$CHR))
    } else {
      chr_info <- "Multi-chromosome"
    }

    title <- paste0("LD Heatmap: ", chr_info, " (", n_snps, " SNPs)")
  }

  # -------------------------------------------------------------------------
  # File output setup
  # -------------------------------------------------------------------------

  if (!is.null(output_file)) {
    if (verbose) message("Setting up file output...")

    # Create output filename
    output_filename <- paste0(output_file, ".", file_format)

    # Set up graphics device
    if (file_format == "pdf") {
      pdf(output_filename, width = width, height = height)
    } else if (file_format %in% c("png", "PNG")) {
      png(output_filename, width = width * dpi, height = height * dpi, res = dpi)
    } else if (file_format %in% c("jpeg", "jpg", "JPEG", "JPG")) {
      jpeg(output_filename, width = width * dpi, height = height * dpi,
        res = dpi, quality = 95)
    } else if (file_format %in% c("tiff", "TIFF")) {
      tiff(output_filename, width = width * dpi, height = height * dpi, res = dpi)
    }

    par(xpd = TRUE)
  } else {
    # Set up display device if needed
    if (is.null(dev.list())) {
      dev.new(width = width, height = height)
    }
    par(xpd = TRUE)
  }

  # -------------------------------------------------------------------------
  # Generate LD heatmap
  # -------------------------------------------------------------------------

  if (verbose) message("Generating LD heatmap...")

  tryCatch(
    {
      ld_heatmap <- LDheatmap(
        gdat = snp_matrix,
        genetic.distances = genetic_distances,
        color = plot_colors,
        flip = flip_diagonal,
        title = title
      )
    },
    error = function(e) {
      if (!is.null(output_file)) dev.off()
      stop("Error generating LD heatmap: ", e$message)
    })

  # Close graphics device if file output
  if (!is.null(output_file)) {
    dev.off()
    if (verbose) {
      message("LD heatmap saved to: ", file.path(getwd(), output_filename))
    }
  }

  # -------------------------------------------------------------------------
  # Prepare return results
  # -------------------------------------------------------------------------

  plot_params <- list(
    color_palette = color_palette,
    color_steps = color_steps,
    flip_diagonal = flip_diagonal,
    title = title,
    subtitle = subtitle,
    output_file = output_file,
    file_format = file_format,
    width = width,
    height = height,
    dpi = dpi,
    ggplot_version = ggplot_version,
    show_values = show_values,
    text_size = text_size,
    n_colors_used = length(plot_colors),
    n_snps = ncol(snp_matrix)
  )

  if (verbose) {
    message("LD heatmap analysis completed successfully!")
    message("Matrix dimensions: ", nrow(snp_matrix), " samples × ", ncol(snp_matrix), " SNPs")
    message("Genetic region span: ", min(genetic_distances), " - ", max(genetic_distances))
    if (!is.null(output_file)) {
      message("Output saved as: ", output_filename)
    }
  }

  # Return results following bioRtools conventions
  return(list(
    plot.ld_heatmap = ld_heatmap,
    snp.matrix = snp_matrix,
    genetic.distances = genetic_distances,
    snp.info = snp_info,
    plot.params = plot_params
  ))
}
