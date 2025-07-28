#' Create Manhattan Plot for Genome-Wide Association Study (GWAS) Results
#'
#' This function creates a publication-ready Manhattan plot for visualizing
#' genome-wide association study (GWAS) results. The plot displays -log10(p-values)
#' across chromosomes, with each chromosome represented by alternating colors.
#' Significant SNPs are highlighted and can be labeled. This visualization is
#' essential for identifying genomic regions associated with traits or diseases.
#'
#' @param df Data frame containing GWAS results with the following required columns:
#'   \itemize{
#'     \item \code{chr}: Chromosome number (numeric or character)
#'     \item \code{posi} or \code{bp} or \code{pos}: Genomic position (base pairs)
#'     \item \code{p} or \code{pvalue} or \code{P}: P-value from association test
#'     \item \code{snp} or \code{SNP} or \code{rsid}: SNP identifier (optional, for labeling)
#'   }
#'   Additional optional columns:
#'   \itemize{
#'     \item \code{label}: Custom labels for significant SNPs
#'     \item \code{gene}: Nearest gene name for annotation
#'   }
#' @param p.threshold Significance threshold for drawing horizontal reference line
#'   (default: Bonferroni correction = 0.05 / number of SNPs). Can be set to
#'   custom values like 5e-8 (genome-wide significance) or 1e-5 (suggestive)
#' @param highlight.threshold P-value threshold for highlighting significant SNPs
#'   (default: same as p.threshold). SNPs below this threshold will be emphasized
#' @param label.threshold P-value threshold for labeling SNPs (default: p.threshold).
#'   Only SNPs more significant than this threshold will be labeled with text
#' @param chr.colors Vector of colors for alternating chromosomes (default: blue and orange).
#'   Should have at least 2 colors that will be recycled across chromosomes
#' @param point.size.range Range of point sizes for SNPs (default: c(0.5, 3)).
#'   Points are sized by -log10(p-value) within this range
#' @param max.overlaps Maximum number of overlapping labels allowed (default: 20).
#'   Helps prevent overcrowded labels in regions of high significance
#' @param title Plot title (default: "Manhattan Plot"). Set to NULL to remove title
#' @param subtitle Plot subtitle (default: "Genome-wide Association Study Results")
#'
#' @return A ggplot2 object representing the Manhattan plot. Can be further
#'   customized with additional ggplot2 layers, themes, or saved using ggsave()
#'
#' @details
#' \strong{Manhattan Plot Components:}
#' \itemize{
#'   \item \strong{X-axis}: Genomic position across all chromosomes
#'   \item \strong{Y-axis}: -log10(p-value) from association tests
#'   \item \strong{Colors}: Alternating colors distinguish chromosomes
#'   \item \strong{Point size}: Proportional to -log10(p-value)
#'   \item \strong{Reference line}: Significance threshold (typically genome-wide significance)
#'   \item \strong{Labels}: Significant SNPs can be annotated with identifiers
#' }
#' 
#' \strong{Standard Significance Thresholds:}
#' \itemize{
#'   \item \strong{Genome-wide significance}: p < 5e-8 (-log10(p) > 7.3)
#'   \item \strong{Suggestive significance}: p < 1e-5 (-log10(p) > 5)
#'   \item \strong{Bonferroni correction}: p < 0.05/N where N = number of SNPs tested
#' }
#' 
#' \strong{Data Preparation Tips:}
#' \itemize{
#'   \item Ensure chromosome names are consistent (1-22, X, Y or chr1-chr22, chrX, chrY)
#'   \item Remove SNPs with p-value = 0 (set to minimum detectable p-value)
#'   \item Consider filtering SNPs with MAF < 0.01 or poor imputation quality
#'   \item Large datasets (>1M SNPs) may require thinning for visualization
#' }
#'
#' @note
#' \itemize{
#'   \item For large datasets, consider filtering or thinning SNPs for better performance
#'   \item The plot automatically scales Y-axis based on the most significant SNP
#'   \item Labels are placed using ggrepel to avoid overlapping
#'   \item Default colors follow scientific publication standards
#' }
#'
#' @references
#' Gibson, G. (2010). Hints of hidden heritability in GWAS. Nature Genetics, 42(7), 558-560.
#' 
#' Visscher, P. M., et al. (2017). 10 Years of GWAS Discovery: Biology, Function, and Translation. 
#' American Journal of Human Genetics, 101(1), 5-22.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{find_degs_deseq2}} for gene-level association analysis
#'   \item \code{\link[ggplot2]{ggplot}} for plot customization
#'   \item \code{\link[ggrepel]{geom_label_repel}} for label positioning
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' library(ggplot2)
#' 
#' # Example 1: Basic Manhattan plot with simulated GWAS data
#' \dontrun{
#' # Generate simulated GWAS data
#' set.seed(20250312)
#' 
#' simulated_gwas <- normentR::simulateGWAS(nSNPs = 1e6, nSigCols = 3) %>%
#'   janitor::clean_names() %>%
#'   dplyr::rename(posi = bp) %>%
#'   dplyr::mutate(
#'     # Add labels for genome-wide significant SNPs
#'     label = case_when(
#'       p < 5e-8 ~ snp,  # Genome-wide significance
#'       TRUE ~ NA
#'     )
#'   )
#' 
#' # Create basic Manhattan plot
#' manhattan_basic <- manhattan_plot(simulated_gwas)
#' print(manhattan_basic)
#' 
#' # Summary of results
#' print(paste("Total SNPs:", nrow(simulated_gwas)))
#' print(paste("Genome-wide significant SNPs (p < 5e-8):", 
#'             sum(simulated_gwas$p < 5e-8, na.rm = TRUE)))
#' print(paste("Suggestive SNPs (p < 1e-5):", 
#'             sum(simulated_gwas$p < 1e-5, na.rm = TRUE)))
#' }
#' 
#' # Example 2: Custom significance thresholds
#' \dontrun{
#' # Manhattan plot with genome-wide significance threshold
#' manhattan_gwas <- manhattan_plot(
#'   df = simulated_gwas,
#'   p.threshold = 5e-8,           # Genome-wide significance
#'   highlight.threshold = 5e-8,    # Highlight genome-wide significant
#'   label.threshold = 1e-10,       # Label only very significant SNPs
#'   title = "GWAS: Height Association",
#'   subtitle = "Genome-wide significance threshold: p < 5×10⁻⁸"
#' )
#' print(manhattan_gwas)
#' 
#' # Add custom annotations
#' manhattan_annotated <- manhattan_gwas +
#'   annotate("text", x = Inf, y = Inf, 
#'            label = paste("N =", format(nrow(simulated_gwas), big.mark = ",")),
#'            hjust = 1.1, vjust = 1.5, size = 3) +
#'   labs(caption = "Data: Simulated GWAS results")
#' 
#' print(manhattan_annotated)
#' }
#' 
#' # Example 3: Multiple significance levels
#' \dontrun{
#' # Create plot with both genome-wide and suggestive thresholds
#' manhattan_multi <- manhattan_plot(
#'   df = simulated_gwas,
#'   p.threshold = 5e-8,
#'   chr.colors = c("#E31A1C", "#1F78B4"),  # Red and blue
#'   point.size.range = c(0.3, 2.5)
#' ) +
#'   # Add suggestive significance line
#'   geom_hline(yintercept = -log10(1e-5), 
#'              color = "grey60", linetype = "dotted", alpha = 0.8) +
#'   # Add annotations for thresholds
#'   annotate("text", x = 0, y = -log10(5e-8) + 0.5, 
#'            label = "Genome-wide significance", 
#'            hjust = 0, size = 3, color = "grey40") +
#'   annotate("text", x = 0, y = -log10(1e-5) + 0.5, 
#'            label = "Suggestive significance", 
#'            hjust = 0, size = 3, color = "grey60")
#' 
#' print(manhattan_multi)
#' }
#' 
#' # Example 4: Real-world data format example
#' # This shows how to prepare your own GWAS data
#' sample_gwas_data <- data.frame(
#'   chr = rep(1:22, each = 1000),
#'   posi = rep(seq(1e6, 100e6, length.out = 1000), 22),
#'   snp = paste0("rs", 1:22000),
#'   p = runif(22000, min = 1e-12, max = 0.5),
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Add some significant hits
#' sample_gwas_data$p[c(1000, 5000, 15000)] <- c(1e-9, 3e-8, 2e-11)
#' sample_gwas_data$label <- ifelse(sample_gwas_data$p < 5e-8, 
#'                                  sample_gwas_data$snp, NA)
#' 
#' \dontrun{
#' # Create Manhattan plot
#' real_world_plot <- manhattan_plot(
#'   df = sample_gwas_data,
#'   p.threshold = 5e-8,
#'   title = "Real-world GWAS Example",
#'   subtitle = "22 chromosomes, 22,000 SNPs"
#' )
#' print(real_world_plot)
#' }
#' 
#' # Example 5: Customization and export
#' \dontrun{
#' # Highly customized Manhattan plot
#' custom_manhattan <- manhattan_plot(
#'   df = simulated_gwas,
#'   p.threshold = 5e-8,
#'   chr.colors = c("#2E8B57", "#CD853F", "#4682B4", "#D2691E"),  # 4 colors
#'   point.size.range = c(0.2, 4),
#'   max.overlaps = 50,
#'   title = NULL  # Remove title for publication
#' ) +
#'   theme_minimal() +
#'   theme(
#'     panel.grid.major.y = element_line(color = "grey90", size = 0.5),
#'     panel.grid.minor.y = element_blank(),
#'     panel.grid.major.x = element_blank(),
#'     axis.title.y = element_text(size = 14),
#'     axis.text.x = element_text(size = 10),
#'     axis.text.y = element_text(size = 10),
#'     plot.margin = margin(20, 20, 20, 20)
#'   ) +
#'   labs(y = expression(-log[10](italic(P))))
#' 
#' print(custom_manhattan)
#' 
#' # Save high-resolution plot
#' ggsave("manhattan_plot.png", custom_manhattan, 
#'        width = 12, height = 6, dpi = 300, bg = "white")
#' 
#' # Save as PDF for publications
#' ggsave("manhattan_plot.pdf", custom_manhattan, 
#'        width = 12, height = 6, device = "pdf")
#' }
#' 
#' # Example 6: Quality control and diagnostics
#' \dontrun{
#' # Check data quality before plotting
#' qc_summary <- simulated_gwas %>%
#'   summarise(
#'     total_snps = n(),
#'     chromosomes = n_distinct(chr),
#'     min_p = min(p, na.rm = TRUE),
#'     max_p = max(p, na.rm = TRUE),
#'     zero_p_values = sum(p == 0, na.rm = TRUE),
#'     na_p_values = sum(is.na(p)),
#'     genome_wide_sig = sum(p < 5e-8, na.rm = TRUE),
#'     suggestive_sig = sum(p < 1e-5, na.rm = TRUE)
#'   )
#' 
#' print("GWAS Data Quality Summary:")
#' print(qc_summary)
#' 
#' # Check genomic inflation
#' observed_chi2 <- qchisq(1 - simulated_gwas$p, df = 1)
#' lambda_gc <- median(observed_chi2, na.rm = TRUE) / qchisq(0.5, df = 1)
#' print(paste("Genomic inflation factor (λ):", round(lambda_gc, 3)))
#' 
#' if (lambda_gc > 1.1) {
#'   warning("Genomic inflation detected (λ > 1.1). Consider population stratification correction.")
#' }
#' }
#'
manhattan_plot <- function(df, p.threshold = NA, highlight.threshold = NULL, 
                          label.threshold = NULL, chr.colors = c("#1F77B4", "#FF7F0E"),
                          point.size.range = c(0.5, 3), max.overlaps = 20,
                          title = "Manhattan Plot", 
                          subtitle = "Genome-wide Association Study Results") {
  
  # Input validation
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  
  if (nrow(df) == 0) {
    stop("'df' cannot be empty")
  }
  
  # Standardize column names - check for common variants
  col_mappings <- list(
    chr = c("chr", "CHR", "chromosome", "Chromosome", "CHROM"),
    posi = c("posi", "bp", "pos", "BP", "POS", "position", "Position"),
    p = c("p", "pvalue", "P", "PVALUE", "P_VALUE", "p_value"),
    snp = c("snp", "SNP", "rsid", "RSID", "rs", "RS", "marker", "MARKER"),
    label = c("label", "LABEL", "annotation", "ANNOTATION")
  )
  
  # Find existing columns and standardize names
  standardized_df <- df
  for (standard_name in names(col_mappings)) {
    possible_names <- col_mappings[[standard_name]]
    found_col <- intersect(possible_names, names(df))
    
    if (length(found_col) > 0) {
      # Use the first matching column and rename it
      standardized_df[[standard_name]] <- df[[found_col[1]]]
      if (found_col[1] != standard_name) {
        # Remove original column if it was renamed
        standardized_df[[found_col[1]]] <- NULL
      }
    }
  }
  
  # Check for required columns
  required_cols <- c("chr", "posi", "p")
  missing_cols <- setdiff(required_cols, names(standardized_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nExpected chromosome column: chr/CHR/chromosome",
               "\nExpected position column: posi/bp/pos/BP/position", 
               "\nExpected p-value column: p/pvalue/P/PVALUE"))
  }
  
  # Validate column types and content
  if (!is.numeric(standardized_df$posi)) {
    stop("Position column must be numeric")
  }
  
  if (!is.numeric(standardized_df$p)) {
    stop("P-value column must be numeric")
  }
  
  # Handle chromosome column (numeric or character)
  if (is.character(standardized_df$chr)) {
    # Remove "chr" prefix if present and convert to numeric where possible
    standardized_df$chr <- gsub("^chr", "", standardized_df$chr, ignore.case = TRUE)
    # Handle X, Y chromosomes
    standardized_df$chr <- case_when(
      toupper(standardized_df$chr) == "X" ~ "23",
      toupper(standardized_df$chr) == "Y" ~ "24", 
      toupper(standardized_df$chr) == "MT" ~ "25",
      toupper(standardized_df$chr) == "M" ~ "25",
      TRUE ~ standardized_df$chr
    )
  }
  
  # Convert chromosome to numeric
  standardized_df$chr <- as.numeric(standardized_df$chr)
  
  # Data quality checks
  if (any(is.na(standardized_df$chr))) {
    stop("Invalid chromosome values detected. Chromosomes should be 1-22, X, Y, or MT")
  }
  
  if (any(standardized_df$p <= 0, na.rm = TRUE)) {
    warning("P-values <= 0 detected. Setting to minimum detectable p-value (1e-300)")
    standardized_df$p[standardized_df$p <= 0] <- 1e-300
  }
  
  if (any(standardized_df$p > 1, na.rm = TRUE)) {
    warning("P-values > 1 detected. These will be set to 1")
    standardized_df$p[standardized_df$p > 1] <- 1
  }
  
  # Remove rows with missing essential data
  complete_rows <- complete.cases(standardized_df[, c("chr", "posi", "p")])
  if (sum(complete_rows) < nrow(standardized_df)) {
    warning(paste("Removing", nrow(standardized_df) - sum(complete_rows), 
                  "rows with missing chromosome, position, or p-value data"))
    standardized_df <- standardized_df[complete_rows, ]
  }
  
  if (nrow(standardized_df) == 0) {
    stop("No complete data rows remain after quality filtering")
  }
  
  # Set default thresholds
  if (is.na(p.threshold)) {
    p.threshold <- 0.05 / nrow(standardized_df)  # Bonferroni correction
  }
  
  if (is.null(highlight.threshold)) {
    highlight.threshold <- p.threshold
  }
  
  if (is.null(label.threshold)) {
    label.threshold <- p.threshold
  }
  
  # Validate thresholds
  if (!is.numeric(p.threshold) || p.threshold <= 0 || p.threshold > 1) {
    stop("'p.threshold' must be a number between 0 and 1")
  }
  
  # Calculate cumulative positions for continuous x-axis
  tryCatch({
    chr_data <- standardized_df %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(max_posi = max(posi, na.rm = TRUE), .groups = 'drop') %>%
      dplyr::arrange(chr) %>%
      dplyr::mutate(
        posi_add = dplyr::lag(cumsum(max_posi), default = 0),
        chr_center = posi_add + max_posi / 2
      )
  }, error = function(e) {
    stop(paste("Error calculating cumulative positions:", e$message))
  })
  
  # Add cumulative positions to main data
  gwas_data <- standardized_df %>%
    dplyr::left_join(chr_data[, c("chr", "posi_add")], by = "chr") %>%
    dplyr::mutate(posi_cum = posi + posi_add)
  
  # Calculate y-axis limits
  min_p <- min(gwas_data$p, na.rm = TRUE)
  max_log10_p <- -log10(min_p)
  y_max <- ceiling(max_log10_p) + 1
  
  # Create labels for significant SNPs
  if ("label" %in% names(gwas_data)) {
    # Use existing labels
    gwas_data$plot_label <- ifelse(gwas_data$p <= label.threshold & !is.na(gwas_data$label),
                                  gwas_data$label, NA)
  } else if ("snp" %in% names(gwas_data)) {
    # Create labels from SNP IDs for significant hits
    gwas_data$plot_label <- ifelse(gwas_data$p <= label.threshold, gwas_data$snp, NA)
  } else {
    # No labels available
    gwas_data$plot_label <- NA
  }
  
  # Prepare colors (ensure sufficient colors for all chromosomes)
  n_chrs <- length(unique(gwas_data$chr))
  if (length(chr.colors) < n_chrs) {
    chr.colors <- rep(chr.colors, ceiling(n_chrs / length(chr.colors)))
  }
  chr.colors <- chr.colors[1:n_chrs]
  
  # Create the Manhattan plot
  tryCatch({
    manhattan_plot <- gwas_data %>%
      ggplot2::ggplot(ggplot2::aes(x = posi_cum, y = -log10(p))) +
      
      # Add significance threshold line
      ggplot2::geom_hline(yintercept = -log10(p.threshold), 
                         color = "red", linetype = "dashed", alpha = 0.7) +
      
      # Add points colored by chromosome
      ggplot2::geom_point(ggplot2::aes(color = factor(chr), size = -log10(p)), 
                         alpha = 0.6) +
      
      # Add labels for significant SNPs (if any exist)
      {if (any(!is.na(gwas_data$plot_label))) {
        ggrepel::geom_label_repel(
          ggplot2::aes(label = plot_label),
          data = gwas_data[!is.na(gwas_data$plot_label), ],
          size = 3,
          box.padding = 0.3,
          point.padding = 0.2,
          segment.color = "grey50",
          segment.size = 0.3,
          max.overlaps = max.overlaps,
          min.segment.length = 0
        )
      }} +
      
      # Customize scales
      ggplot2::scale_x_continuous(
        breaks = chr_data$chr_center,
        labels = chr_data$chr,
        expand = ggplot2::expansion(mult = 0.01)
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.05)),
        limits = c(0, y_max),
        breaks = scales::pretty_breaks(n = 6)
      ) +
      ggplot2::scale_color_manual(values = chr.colors) +
      ggplot2::scale_size_continuous(range = point.size.range, guide = "none") +
      
      # Labels and theme
      ggplot2::labs(
        x = "Chromosome",
        y = expression(-log[10](italic(P))),
        title = title,
        subtitle = subtitle
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 9),
        axis.text.y = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10),
        panel.border = ggplot2::element_rect(color = "grey80", fill = NA, size = 0.5)
      )
  }, error = function(e) {
    stop(paste("Error creating Manhattan plot:", e$message))
  })
  
  # Add summary statistics as attributes
  n_total <- nrow(gwas_data)
  n_significant <- sum(gwas_data$p < highlight.threshold, na.rm = TRUE)
  n_labeled <- sum(!is.na(gwas_data$plot_label))
  
  attr(manhattan_plot, "n_snps") <- n_total
  attr(manhattan_plot, "n_significant") <- n_significant
  attr(manhattan_plot, "n_labeled") <- n_labeled
  attr(manhattan_plot, "p_threshold") <- p.threshold
  attr(manhattan_plot, "min_p") <- min_p
  attr(manhattan_plot, "chromosomes") <- sort(unique(gwas_data$chr))
  
  # Print summary if interactive
  if (interactive()) {
    cat("Manhattan Plot Summary:\n")
    cat("======================\n")
    cat("Total SNPs plotted:", format(n_total, big.mark = ","), "\n")
    cat("Chromosomes:", paste(sort(unique(gwas_data$chr)), collapse = ", "), "\n")
    cat("Significance threshold:", format(p.threshold, scientific = TRUE), "\n")
    cat("Significant SNPs:", format(n_significant, big.mark = ","), 
        paste0("(", round(100 * n_significant / n_total, 2), "%)"), "\n")
    cat("Labeled SNPs:", n_labeled, "\n")
    cat("Most significant p-value:", format(min_p, scientific = TRUE), "\n")
    
    if (n_significant > 0) {
      # Show top hits
      top_hits <- gwas_data %>%
        dplyr::filter(p < highlight.threshold) %>%
        dplyr::arrange(p) %>%
        dplyr::select(chr, posi, p, dplyr::any_of(c("snp", "plot_label"))) %>%
        head(5)
      
      cat("\nTop 5 significant hits:\n")
      print(top_hits)
    }
    cat("\n")
  }
  
  return(manhattan_plot)
}