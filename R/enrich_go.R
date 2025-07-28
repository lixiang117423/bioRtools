#' Perform Gene Ontology (GO) Enrichment Analysis
#'
#' This function performs Gene Ontology enrichment analysis on a set of 
#' differentially expressed genes using a custom GO annotation database.
#' It identifies significantly over-represented GO terms and calculates
#' enrichment statistics with multiple testing correction.
#'
#' @param gene Character vector containing gene identifiers (e.g., gene symbols,
#'   Ensembl IDs, or other identifiers) that correspond to differentially 
#'   expressed genes for enrichment analysis
#' @param go.db Data frame containing GO annotation information. Must include
#'   the following columns:
#'   \itemize{
#'     \item \code{gene}: Gene identifiers matching those in the gene parameter
#'     \item \code{go.id}: GO term identifiers (e.g., "GO:0008150")
#'     \item \code{go.term}: GO term descriptions (e.g., "biological_process")
#'     \item \code{go.ontology}: Optional. GO ontology category (BP, MF, CC)
#'   }
#' @param pAdjustMethod Method for multiple testing correction. Options include:
#'   \itemize{
#'     \item \code{"BH"} (default): Benjamini-Hochberg false discovery rate
#'     \item \code{"bonferroni"}: Bonferroni correction
#'     \item \code{"holm"}: Holm's step-down procedure
#'     \item \code{"hochberg"}: Hochberg's procedure
#'     \item \code{"hommel"}: Hommel's procedure  
#'     \item \code{"BY"}: Benjamini-Yekutieli procedure
#'     \item \code{"fdr"}: Synonym for "BH"
#'     \item \code{"none"}: No correction
#'   }
#' @param p.adjust Significance threshold for adjusted p-values (default: 0.05).
#'   Only GO terms with adjusted p-values below this threshold will be returned
#' @param min.gene.set Minimum number of genes required in a GO term for testing
#'   (default: 3). GO terms with fewer genes will be excluded
#' @param max.gene.set Maximum number of genes allowed in a GO term for testing
#'   (default: length of input gene list). Very large GO terms can be excluded
#'
#' @return A data frame containing enrichment results with columns:
#'   \itemize{
#'     \item \code{ID}: GO term identifier
#'     \item \code{Description}: GO term description
#'     \item \code{GeneRatio}: Proportion of input genes in this GO term (numeric)
#'     \item \code{BgRatio}: Proportion of background genes in this GO term
#'     \item \code{pvalue}: Raw p-value from hypergeometric test
#'     \item \code{p.adjust}: Adjusted p-value using specified method
#'     \item \code{qvalue}: Q-value (if calculated)
#'     \item \code{geneID}: Gene identifiers contributing to this term (separated by "/")
#'     \item \code{Count}: Number of input genes in this GO term
#'     \item \code{gene.count}: Number of input genes in this GO term
#'     \item \code{total.genes}: Total number of input genes tested
#'   }
#'
#' @details
#' This function uses the clusterProfiler package to perform over-representation
#' analysis (ORA) based on the hypergeometric distribution. The analysis tests
#' whether genes in your input list are significantly over-represented in specific
#' GO terms compared to what would be expected by chance.
#' 
#' The function converts the character-based GeneRatio (e.g., "5/100") to numeric
#' format for easier downstream analysis and plotting.
#'
#' @note
#' \itemize{
#'   \item Ensure gene identifiers in the input list match those in the GO database
#'   \item Remove duplicate genes from the input list for accurate statistics
#'   \item Consider the completeness of your GO annotation database
#'   \item Very small gene sets may not yield meaningful enrichment results
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' 
#' # Example 1: Basic GO enrichment analysis
#' # Note: These examples assume you have the required data objects
#' \dontrun{
#' data(df.rnaseq.go)     # GO annotation database
#' data(df.rnaseq.degs)   # Differentially expressed genes
#' 
#' # Basic enrichment analysis
#' go_results <- enrich_go(
#'   gene = df.rnaseq.degs$gene,
#'   go.db = df.rnaseq.go
#' )
#' 
#' print(head(go_results))
#' print(paste("Found", nrow(go_results), "enriched GO terms"))
#' }
#' 
#' # Example 2: Custom parameters
#' \dontrun{
#' # More stringent analysis with Bonferroni correction
#' go_results_strict <- enrich_go(
#'   gene = df.rnaseq.degs$gene,
#'   go.db = df.rnaseq.go,
#'   pAdjustMethod = "bonferroni",
#'   p.adjust = 0.01,
#'   min.gene.set = 5
#' )
#' 
#' print(go_results_strict)
#' }
#' 
#' # Example 3: Creating a simple GO database format
#' # This shows the required structure for go.db parameter
#' sample_go_db <- data.frame(
#'   gene = c("GENE1", "GENE2", "GENE3", "GENE1", "GENE4", "GENE5"),
#'   go.id = c("GO:0008150", "GO:0008150", "GO:0008150", 
#'             "GO:0003674", "GO:0003674", "GO:0005575"),
#'   go.term = c("biological_process", "biological_process", "biological_process",
#'               "molecular_function", "molecular_function", "cellular_component"),
#'   go.ontology = c("BP", "BP", "BP", "MF", "MF", "CC")
#' )
#' 
#' sample_genes <- c("GENE1", "GENE2", "GENE4")
#' 
#' # Run enrichment (will likely find no significant terms due to small size)
#' \dontrun{
#' sample_results <- enrich_go(
#'   gene = sample_genes,
#'   go.db = sample_go_db,
#'   p.adjust = 1.0  # Accept all terms for demonstration
#' )
#' print(sample_results)
#' }
#' 
#' # Example 4: Filtering by GO ontology
#' \dontrun{
#' # Filter GO database to only biological processes
#' bp_go_db <- df.rnaseq.go %>%
#'   filter(go.ontology == "BP" | is.na(go.ontology))
#' 
#' bp_results <- enrich_go(
#'   gene = df.rnaseq.degs$gene,
#'   go.db = bp_go_db,
#'   pAdjustMethod = "BH",
#'   p.adjust = 0.05
#' )
#' 
#' print(paste("Biological process terms found:", nrow(bp_results)))
#' }
#'
enrich_go <- function(gene, go.db, pAdjustMethod = "BH", p.adjust = 0.05, 
                      min.gene.set = 3, max.gene.set = NULL) {
  
  # Input validation
  if (!is.character(gene) || length(gene) == 0) {
    stop("'gene' must be a non-empty character vector")
  }
  
  if (!is.data.frame(go.db) || nrow(go.db) == 0) {
    stop("'go.db' must be a non-empty data frame")
  }
  
  # Check required columns in go.db
  required_cols <- c("gene", "go.id", "go.term")
  missing_cols <- setdiff(required_cols, names(go.db))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in go.db:", paste(missing_cols, collapse = ", ")))
  }
  
  # Validate adjustment method
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!pAdjustMethod %in% valid_methods) {
    stop(paste("Invalid pAdjustMethod. Choose from:", paste(valid_methods, collapse = ", ")))
  }
  
  # Validate p.adjust threshold
  if (!is.numeric(p.adjust) || length(p.adjust) != 1 || p.adjust <= 0 || p.adjust > 1) {
    stop("'p.adjust' must be a single numeric value between 0 and 1")
  }
  
  # Validate gene set size parameters
  if (!is.numeric(min.gene.set) || min.gene.set < 1) {
    stop("'min.gene.set' must be a positive integer")
  }
  
  if (is.null(max.gene.set)) {
    max.gene.set <- length(gene)
  } else if (!is.numeric(max.gene.set) || max.gene.set < min.gene.set) {
    stop("'max.gene.set' must be numeric and >= min.gene.set")
  }
  
  # Remove duplicates from input genes and warn if found
  original_gene_count <- length(gene)
  gene <- unique(gene[!is.na(gene)])
  
  if (length(gene) != original_gene_count) {
    warning(paste("Removed", original_gene_count - length(gene), 
                  "duplicate or missing gene identifiers"))
  }
  
  if (length(gene) == 0) {
    stop("No valid gene identifiers remaining after cleanup")
  }
  
  # Clean GO database
  go.db_clean <- go.db %>%
    dplyr::filter(!is.na(gene), !is.na(go.id), !is.na(go.term)) %>%
    dplyr::distinct()
  
  if (nrow(go.db_clean) == 0) {
    stop("No valid annotations remain in go.db after removing missing values")
  }
  
  # Check overlap between input genes and GO database
  db_genes <- unique(go.db_clean$gene)
  overlapping_genes <- intersect(gene, db_genes)
  
  if (length(overlapping_genes) == 0) {
    stop("No overlap between input genes and genes in GO database")
  }
  
  if (length(overlapping_genes) < length(gene)) {
    warning(paste("Only", length(overlapping_genes), "of", length(gene), 
                  "input genes found in GO database"))
  }
  
  # Prepare TERM2GENE and TERM2NAME mappings
  term2gene <- go.db_clean %>%
    dplyr::select(go.id, gene) %>%
    dplyr::distinct()
  
  term2name <- go.db_clean %>%
    dplyr::select(go.id, go.term) %>%
    dplyr::distinct() %>%
    # Handle cases where one GO ID has multiple descriptions
    dplyr::group_by(go.id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  
  # Check gene set sizes and filter if necessary
  term_sizes <- term2gene %>%
    dplyr::count(go.id, name = "term_size") %>%
    dplyr::filter(term_size >= min.gene.set, term_size <= max.gene.set)
  
  if (nrow(term_sizes) == 0) {
    stop(paste("No GO terms meet the size criteria (", min.gene.set, "-", max.gene.set, " genes)"))
  }
  
  # Filter TERM2GENE to include only terms meeting size criteria
  term2gene_filtered <- term2gene %>%
    dplyr::filter(go.id %in% term_sizes$go.id)
  
  # Perform enrichment analysis
  tryCatch({
    go_enrich <- clusterProfiler::enricher(
      gene = gene,
      TERM2GENE = term2gene_filtered,
      TERM2NAME = term2name,
      qvalueCutoff = p.adjust,
      pAdjustMethod = pAdjustMethod,
      minGSSize = min.gene.set,
      maxGSSize = max.gene.set
    )
  }, error = function(e) {
    stop(paste("Error in clusterProfiler::enricher:", e$message))
  })
  
  # Extract and process results
  if (nrow(go_enrich@result) == 0) {
    warning("No significantly enriched GO terms found with current parameters")
    # Return empty data frame with expected column structure
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
      total.genes = integer(0)
    )
    return(empty_result)
  }
  
  # Process results
  go_results <- go_enrich@result
  
  # Convert GeneRatio from character to numeric
  go_results$gene.count <- sapply(go_results$GeneRatio, function(x) {
    as.numeric(strsplit(as.character(x), "/")[[1]][1])
  })
  
  go_results$total.genes <- sapply(go_results$GeneRatio, function(x) {
    as.numeric(strsplit(as.character(x), "/")[[1]][2])
  })
  
  # Calculate numeric GeneRatio
  go_results$GeneRatio <- go_results$gene.count / go_results$total.genes
  
  # Reorder columns for better readability
  go_results <- go_results %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, 
                  qvalue, geneID, Count, gene.count, total.genes, 
                  dplyr::everything()) %>%
    dplyr::arrange(p.adjust, pvalue)
  
  # Add metadata as attributes
  attr(go_results, "input_genes") <- length(gene)
  attr(go_results, "genes_in_database") <- length(overlapping_genes)
  attr(go_results, "total_go_terms_tested") <- nrow(term_sizes)
  attr(go_results, "adjustment_method") <- pAdjustMethod
  attr(go_results, "significance_threshold") <- p.adjust
  
  # Print summary if interactive
  if (interactive() && nrow(go_results) > 0) {
    cat("GO Enrichment Analysis Summary:\n")
    cat("Input genes:", length(gene), "\n")
    cat("Genes found in database:", length(overlapping_genes), "\n")
    cat("GO terms tested:", nrow(term_sizes), "\n")
    cat("Significantly enriched terms:", nrow(go_results), "\n")
    cat("Adjustment method:", pAdjustMethod, "\n")
    cat("Significance threshold:", p.adjust, "\n\n")
    
    if (nrow(go_results) > 0) {
      cat("Top 5 enriched terms:\n")
      top_terms <- head(go_results, 5)
      for (i in 1:nrow(top_terms)) {
        cat(sprintf("%s: %s (p.adj = %.2e, genes = %d)\n",
                    top_terms$ID[i], 
                    substr(top_terms$Description[i], 1, 50),
                    top_terms$p.adjust[i],
                    top_terms$gene.count[i]))
      }
    }
  }
  
  return(go_results)
}