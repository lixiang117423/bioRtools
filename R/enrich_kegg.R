#' Perform KEGG Pathway Enrichment Analysis
#'
#' This function performs KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway
#' enrichment analysis on a set of differentially expressed genes using a custom
#' KEGG annotation database. It identifies significantly over-represented metabolic
#' and signaling pathways and calculates enrichment statistics with multiple
#' testing correction.
#'
#' @param gene Character vector containing gene identifiers (e.g., gene symbols,
#'   Entrez IDs, Ensembl IDs, or other identifiers) that correspond to
#'   differentially expressed genes for pathway enrichment analysis
#' @param kegg.db Data frame containing KEGG pathway annotation information.
#'   Must include the following columns:
#'   \itemize{
#'     \item \code{gene}: Gene identifiers matching those in the gene parameter
#'     \item \code{kegg.id}: KEGG pathway identifiers (e.g., "hsa04110", "map00010")
#'     \item \code{kegg.term}: KEGG pathway descriptions (e.g., "Cell cycle", "Glycolysis")
#'     \item \code{kegg.category}: Optional. KEGG pathway categories (e.g., "Metabolism", "Signaling")
#'   }
#' @param pAdjustMethod Method for multiple testing correction. Options include:
#'   \itemize{
#'     \item \code{"BH"} (default): Benjamini-Hochberg false discovery rate
#'     \item \code{"bonferroni"}: Bonferroni correction (most conservative)
#'     \item \code{"holm"}: Holm's step-down procedure
#'     \item \code{"hochberg"}: Hochberg's procedure
#'     \item \code{"hommel"}: Hommel's procedure
#'     \item \code{"BY"}: Benjamini-Yekutieli procedure
#'     \item \code{"fdr"}: Synonym for "BH"
#'     \item \code{"none"}: No correction (not recommended)
#'   }
#' @param p.adjust Significance threshold for adjusted p-values (default: 0.05).
#'   Only KEGG pathways with adjusted p-values below this threshold will be returned
#' @param min.pathway.size Minimum number of genes required in a pathway for testing
#'   (default: 3). Pathways with fewer genes will be excluded from analysis
#' @param max.pathway.size Maximum number of genes allowed in a pathway for testing
#'   (default: 500). Very large pathways can be excluded to focus on specific processes
#'
#' @return A data frame containing pathway enrichment results with columns:
#'   \itemize{
#'     \item \code{ID}: KEGG pathway identifier
#'     \item \code{Description}: KEGG pathway description/name
#'     \item \code{GeneRatio}: Proportion of input genes in this pathway (numeric)
#'     \item \code{BgRatio}: Proportion of background genes in this pathway
#'     \item \code{pvalue}: Raw p-value from hypergeometric test
#'     \item \code{p.adjust}: Adjusted p-value using specified method
#'     \item \code{qvalue}: Q-value (if calculated)
#'     \item \code{geneID}: Gene identifiers contributing to this pathway (separated by "/")
#'     \item \code{Count}: Number of input genes in this pathway
#'     \item \code{gene.count}: Number of input genes in this pathway
#'     \item \code{total.genes}: Total number of input genes tested
#'     \item \code{enrichment.score}: -log10(p.adjust) for visualization
#'   }
#'
#' @details
#' This function uses clusterProfiler to perform over-representation analysis (ORA)
#' based on the hypergeometric distribution. The analysis tests whether genes in
#' your input list are significantly over-represented in specific KEGG pathways
#' compared to what would be expected by chance.
#'
#' KEGG pathways provide insight into:
#' \itemize{
#'   \item Metabolic processes and biochemical reactions
#'   \item Signal transduction cascades
#'   \item Cellular processes and regulation
#'   \item Disease-associated pathways
#'   \item Drug metabolism and interactions
#' }
#'
#' The function automatically converts character-based GeneRatio to numeric format
#' and adds an enrichment score for easier interpretation and visualization.
#'
#' @note
#' \itemize{
#'   \item Ensure gene identifiers in your list match those in the KEGG database
#'   \item KEGG pathway annotations may vary between organisms
#'   \item Consider pathway size when interpreting results (very large pathways may be less specific)
#'   \item For human data, common KEGG ID formats are "hsa" prefixed (e.g., "hsa04110")
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{enrich_go}} for Gene Ontology enrichment analysis
#'   \item \code{\link{find_degs_deseq2}} for differential expression analysis
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#'
#' # Example 1: Basic KEGG pathway enrichment
#' \dontrun{
#' # Load required data (assumes these datasets exist)
#' data(df.rnaseq.kegg)   # KEGG annotation database
#' data(df.rnaseq.degs)   # Differentially expressed genes
#'
#' # Basic pathway enrichment analysis
#' kegg_results <- enrich_kegg(
#'   gene = df.rnaseq.degs$gene,
#'   kegg.db = df.rnaseq.kegg
#' )
#'
#' print(head(kegg_results))
#' print(paste("Found", nrow(kegg_results), "enriched pathways"))
#'
#' # View top enriched pathways
#' top_pathways <- head(kegg_results, 10)
#' print(top_pathways[, c("ID", "Description", "pvalue", "p.adjust", "Count")])
#' }
#'
#' # Example 2: Stringent analysis parameters
#' \dontrun{
#' # More conservative analysis
#' kegg_strict <- enrich_kegg(
#'   gene = df.rnaseq.degs$gene,
#'   kegg.db = df.rnaseq.kegg,
#'   pAdjustMethod = "bonferroni",
#'   p.adjust = 0.01,
#'   min.pathway.size = 5,
#'   max.pathway.size = 200
#' )
#'
#' print(paste("Strict analysis found", nrow(kegg_strict), "pathways"))
#' }
#'
#' # Example 3: Creating a KEGG database format
#' # This demonstrates the required structure for kegg.db parameter
#' sample_kegg_db <- data.frame(
#'   gene = c("GENE1", "GENE2", "GENE3", "GENE1", "GENE4", "GENE5", "GENE6"),
#'   kegg.id = c("hsa04110", "hsa04110", "hsa04110", "hsa00010", "hsa00010", "hsa04210", "hsa04210"),
#'   kegg.term = c("Cell cycle", "Cell cycle", "Cell cycle",
#'     "Glycolysis / Gluconeogenesis", "Glycolysis / Gluconeogenesis",
#'     "Apoptosis", "Apoptosis"),
#'   kegg.category = c("Cell Growth and Death", "Cell Growth and Death", "Cell Growth and Death",
#'     "Metabolism", "Metabolism", "Cell Growth and Death", "Cell Growth and Death")
#' )
#'
#' sample_genes <- c("GENE1", "GENE2", "GENE4", "GENE6")
#'
#' \dontrun{
#' # Run pathway analysis with sample data
#' sample_results <- enrich_kegg(
#'   gene = sample_genes,
#'   kegg.db = sample_kegg_db,
#'   p.adjust = 1.0  # Accept all pathways for demonstration
#' )
#' print(sample_results)
#' }
#'
#' # Example 4: Filter by pathway categories
#' \dontrun{
#' # Focus on metabolic pathways only
#' metabolism_kegg <- df.rnaseq.kegg %>%
#'   filter(grepl("Metabolism", kegg.category, ignore.case = TRUE))
#'
#' metabolism_results <- enrich_kegg(
#'   gene = df.rnaseq.degs$gene,
#'   kegg.db = metabolism_kegg,
#'   pAdjustMethod = "BH"
#' )
#'
#' print(paste("Metabolic pathways enriched:", nrow(metabolism_results)))
#'
#' # Focus on signaling pathways
#' signaling_kegg <- df.rnaseq.kegg %>%
#'   filter(grepl("Signal|Signaling", kegg.category, ignore.case = TRUE))
#'
#' signaling_results <- enrich_kegg(
#'   gene = df.rnaseq.degs$gene,
#'   kegg.db = signaling_kegg
#' )
#'
#' print(paste("Signaling pathways enriched:", nrow(signaling_results)))
#' }
#'
#' # Example 5: Visualization preparation
#' \dontrun{
#' # Prepare data for plotting
#' kegg_for_plot <- kegg_results %>%
#'   arrange(desc(enrichment.score)) %>%
#'   head(15) %>%
#'   mutate(
#'     pathway_short = ifelse(nchar(Description) > 40,
#'       paste0(substr(Description, 1, 37), "..."),
#'       Description)
#'   )
#'
#' # The results are ready for ggplot2 visualization
#' print("Data prepared for pathway visualization")
#' print(kegg_for_plot[, c("ID", "pathway_short", "enrichment.score", "Count")])
#' }
#'
enrich_kegg <- function(gene, kegg.db, pAdjustMethod = "BH", p.adjust = 0.05,
                        min.pathway.size = 3, max.pathway.size = 500) {
  # Input validation
  if (!is.character(gene) || length(gene) == 0) {
    stop("'gene' must be a non-empty character vector")
  }

  if (!is.data.frame(kegg.db) || nrow(kegg.db) == 0) {
    stop("'kegg.db' must be a non-empty data frame")
  }

  # Check required columns in kegg.db
  required_cols <- c("gene", "kegg.id", "kegg.term")
  missing_cols <- setdiff(required_cols, names(kegg.db))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in kegg.db:", paste(missing_cols, collapse = ", ")))
  }

  # Validate adjustment method
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!pAdjustMethod %in% valid_methods) {
    stop(paste("Invalid pAdjustMethod. Choose from:", paste(valid_methods, collapse = ", ")))
  }

  # Validate parameters
  if (!is.numeric(p.adjust) || length(p.adjust) != 1 || p.adjust <= 0 || p.adjust > 1) {
    stop("'p.adjust' must be a single numeric value between 0 and 1")
  }

  if (!is.numeric(min.pathway.size) || min.pathway.size < 1) {
    stop("'min.pathway.size' must be a positive integer")
  }

  if (!is.numeric(max.pathway.size) || max.pathway.size < min.pathway.size) {
    stop("'max.pathway.size' must be numeric and >= min.pathway.size")
  }

  # Remove duplicates and missing values from input genes
  original_gene_count <- length(gene)
  gene <- unique(gene[!is.na(gene) & gene != ""])

  if (length(gene) != original_gene_count) {
    warning(paste("Removed", original_gene_count - length(gene),
      "duplicate, missing, or empty gene identifiers"))
  }

  if (length(gene) == 0) {
    stop("No valid gene identifiers remaining after cleanup")
  }

  # Clean KEGG database
  kegg.db_clean <- kegg.db %>%
    dplyr::filter(
      !is.na(gene), !is.na(kegg.id), !is.na(kegg.term),
      gene != "", kegg.id != "", kegg.term != ""
    ) %>%
    dplyr::distinct()

  if (nrow(kegg.db_clean) == 0) {
    stop("No valid annotations remain in kegg.db after removing missing values")
  }

  # Check overlap between input genes and KEGG database
  db_genes <- unique(kegg.db_clean$gene)
  overlapping_genes <- intersect(gene, db_genes)

  if (length(overlapping_genes) == 0) {
    stop("No overlap between input genes and genes in KEGG database")
  }

  if (length(overlapping_genes) < length(gene)) {
    missing_count <- length(gene) - length(overlapping_genes)
    warning(paste("Only", length(overlapping_genes), "of", length(gene),
      "input genes found in KEGG database;", missing_count, "genes missing"))
  }

  # Prepare TERM2GENE and TERM2NAME mappings
  term2gene <- kegg.db_clean %>%
    dplyr::select(kegg.id, gene) %>%
    dplyr::distinct()

  term2name <- kegg.db_clean %>%
    dplyr::select(kegg.id, kegg.term) %>%
    dplyr::distinct() %>%
    # Handle cases where one KEGG ID has multiple descriptions
    dplyr::group_by(kegg.id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  # Check pathway sizes and filter
  pathway_sizes <- term2gene %>%
    dplyr::count(kegg.id, name = "pathway_size") %>%
    dplyr::filter(pathway_size >= min.pathway.size, pathway_size <= max.pathway.size)

  if (nrow(pathway_sizes) == 0) {
    stop(paste("No KEGG pathways meet the size criteria (",
      min.pathway.size, "-", max.pathway.size, " genes)"))
  }

  # Filter TERM2GENE to include only pathways meeting size criteria
  term2gene_filtered <- term2gene %>%
    dplyr::filter(kegg.id %in% pathway_sizes$kegg.id)

  # Check if we have enough pathways to test
  n_pathways_to_test <- length(unique(term2gene_filtered$kegg.id))
  if (n_pathways_to_test == 0) {
    stop("No pathways available for testing after size filtering")
  }

  # Perform enrichment analysis
  tryCatch(
    {
      kegg_enrich <- clusterProfiler::enricher(
        gene = gene,
        TERM2GENE = term2gene_filtered,
        TERM2NAME = term2name,
        qvalueCutoff = p.adjust,
        pAdjustMethod = pAdjustMethod,
        minGSSize = min.pathway.size,
        maxGSSize = max.pathway.size
      )
    },
    error = function(e) {
      stop(paste("Error in clusterProfiler::enricher:", e$message))
    })

  # Extract and process results
  if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
    warning("No significantly enriched KEGG pathways found with current parameters")

    # Return empty data frame with expected structure
    empty_result <- data.frame(
      ID = character(0),
      Description = character(0),
      GeneRatio = numeric(0),
      BgRatio = character(0),
      pvalue = numeric(0),
      p.adjust = numeric(0),
      qvalue = numeric(0),
      geneID = character(0),
      Count = integer(0),
      gene.count = integer(0),
      total.genes = integer(0),
      enrichment.score = numeric(0)
    )

    # Add metadata
    attr(empty_result, "input_genes") <- length(gene)
    attr(empty_result, "genes_in_database") <- length(overlapping_genes)
    attr(empty_result, "pathways_tested") <- n_pathways_to_test

    return(empty_result)
  }

  # Process results
  kegg_results <- kegg_enrich@result

  # Convert GeneRatio from character to numeric and add derived columns
  kegg_results$gene.count <- sapply(kegg_results$GeneRatio, function(x) {
    ratio_parts <- strsplit(as.character(x), "/")[[1]]
    as.numeric(ratio_parts[1])
  })

  kegg_results$total.genes <- sapply(kegg_results$GeneRatio, function(x) {
    ratio_parts <- strsplit(as.character(x), "/")[[1]]
    as.numeric(ratio_parts[2])
  })

  # Calculate numeric GeneRatio
  kegg_results$GeneRatio <- kegg_results$gene.count / kegg_results$total.genes

  # Add enrichment score for visualization (-log10 of adjusted p-value)
  kegg_results$enrichment.score <- -log10(kegg_results$p.adjust)

  # Reorder columns for better readability
  kegg_results <- kegg_results %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust,
      qvalue, geneID, Count, gene.count, total.genes, enrichment.score,
      dplyr::everything()) %>%
    dplyr::arrange(p.adjust, pvalue)

  # Add comprehensive metadata as attributes
  attr(kegg_results, "input_genes") <- length(gene)
  attr(kegg_results, "genes_in_database") <- length(overlapping_genes)
  attr(kegg_results, "pathways_tested") <- n_pathways_to_test
  attr(kegg_results, "adjustment_method") <- pAdjustMethod
  attr(kegg_results, "significance_threshold") <- p.adjust
  attr(kegg_results, "min_pathway_size") <- min.pathway.size
  attr(kegg_results, "max_pathway_size") <- max.pathway.size

  # Print comprehensive summary if interactive
  if (interactive() && nrow(kegg_results) > 0) {
    cat("KEGG Pathway Enrichment Analysis Summary:\n")
    cat("=========================================\n")
    cat("Input genes:", length(gene), "\n")
    cat("Genes found in database:", length(overlapping_genes),
      paste0("(", round(100 * length(overlapping_genes) / length(gene), 1), "%)"), "\n")
    cat("Pathways tested:", n_pathways_to_test, "\n")
    cat("Pathway size range:", min.pathway.size, "-", max.pathway.size, "genes\n")
    cat("Significantly enriched pathways:", nrow(kegg_results), "\n")
    cat("Adjustment method:", pAdjustMethod, "\n")
    cat("Significance threshold:", p.adjust, "\n\n")

    # Show top pathways
    cat("Top enriched pathways:\n")
    top_pathways <- head(kegg_results, 5)
    for (i in 1:nrow(top_pathways)) {
      pathway_name <- top_pathways$Description[i]
      # Truncate long pathway names for display
      if (nchar(pathway_name) > 50) {
        pathway_name <- paste0(substr(pathway_name, 1, 47), "...")
      }

      cat(sprintf("%d. %s: %s\n", i, top_pathways$ID[i], pathway_name))
      cat(sprintf("   p.adj = %.2e, genes = %d/%d, enrichment = %.1f\n",
        top_pathways$p.adjust[i],
        top_pathways$gene.count[i],
        top_pathways$total.genes[i],
        top_pathways$enrichment.score[i]))
    }
    cat("\n")
  }

  return(kegg_results)
}
