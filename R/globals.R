#' Global Variable Declarations for bioRtools Package
#'
#' This file declares global variables used throughout the bioRtools package
#' to avoid R CMD check warnings. These variables are primarily used in
#' non-standard evaluation (NSE) contexts with dplyr, ggplot2, and other
#' tidyverse packages.
#'
#' @name .onLoad
NULL

# =============================================================================
# Global Variable Declarations
# =============================================================================

utils::globalVariables(c(
  # -------------------------------------------------------------------------
  # Core Data Manipulation Variables
  # -------------------------------------------------------------------------

  # Standard data frame columns and identifiers
  "sample", "value", "group", "color", "shape", "size", "alpha",
  "feature", "feature_id", "gene", "taxa", "compound", "pathway",
  "id", "name", "label", "category", "type", "class", "status",

  # Row and column identifiers
  "rowname", "colname", "row_id", "col_id", "index",

  # Generic data manipulation variables
  "x", "y", "z", "n", "count", "total", "mean", "median",
  "min", "max", "sd", "var", "se", "ci", "range",

  # Special tidyverse symbols
  ".", ".data", ".env",

  # -------------------------------------------------------------------------
  # Statistical Analysis Variables
  # -------------------------------------------------------------------------

  # Significance and p-values
  "pvalue", "p.value", "padj", "p.adj", "p.adjust", "pval",
  "qvalue", "q.value", "fdr", "bonferroni",
  "statistic", "stat", "test.stat", "z.score", "t.stat",
  "significance", "sig", "significant", "ns",

  # Effect sizes and estimates
  "estimate", "effect", "effect.size", "coefficient", "coef",
  "beta", "slope", "intercept", "r.squared", "R2", "adj.r.squared",
  "correlation", "cor", "rho", "tau", "kendall", "pearson", "spearman",

  # Confidence intervals and errors
  "lower", "upper", "ci.lower", "ci.upper", "conf.low", "conf.high",
  "std.error", "se", "stderr", "lfcSE", "error", "residual",

  # Model components
  "model", "fit", "fitted", "predicted", "observed", "expected",
  "baseline", "treatment", "control", "reference",

  # -------------------------------------------------------------------------
  # Biological Data Variables
  # -------------------------------------------------------------------------

  # Sample and experimental design
  "Species", "species", "condition", "treatment", "timepoint", "time",
  "batch", "plate", "well", "replicate", "biological.replicate",
  "technical.replicate", "subject", "patient", "individual",
  "age", "sex", "gender", "population", "strain", "genotype",

  # Genomics and transcriptomics
  "Gene", "gene", "gene.name", "gene.id", "ensembl.id", "entrez.id",
  "symbol", "description", "biotype", "chromosome", "start", "end",
  "strand", "length", "gc.content",
  "expression", "counts", "reads", "rpkm", "fpkm", "tpm", "cpm",
  "log2FoldChange", "log2FC", "lfc", "fold.change", "fc",
  "baseMean", "base.mean", "AveExpr", "logCPM",
  "regulation", "up", "down", "upregulated", "downregulated",

  # Functional analysis
  "GO", "go.id", "go.term", "go.ontology", "go.category",
  "KEGG", "kegg.id", "kegg.term", "kegg.pathway", "kegg.category",
  "pathway.name", "term", "description", "ontology", "category",
  "enrichment", "enriched", "depleted", "over.represented",
  "gene.ratio", "bg.ratio", "gene.count", "background.count",

  # Proteomics
  "protein", "protein.id", "peptide", "abundance", "intensity",
  "ratio", "modification", "charge", "mass", "mz",

  # Metabolomics
  "Compound", "compound", "metabolite", "formula", "adduct",
  "rt", "retention.time", "mass.error", "isotope",
  "VIP", "Score", "loading", "contribution",

  # Microbiome and ecology
  "OTU", "otu", "ASV", "asv", "Taxa", "taxa", "taxonomy",
  "phylum", "class", "order", "family", "genus", "species.name",
  "Abundance", "abundance", "relative.abundance", "count",
  "diversity", "richness", "evenness", "shannon", "simpson", "chao1",
  "Group", "community", "microbiome", "sample.type",

  # -------------------------------------------------------------------------
  # Dimensionality Reduction and Visualization
  # -------------------------------------------------------------------------

  # Principal Component Analysis (PCA)
  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
  "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10",
  "Dim1", "Dim2", "Dim3", "Dim4", "Dim5",
  "Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5",
  "variance", "variance.percent", "prop.var", "cumulative.var",
  "eigenvalue", "eigen.value", "contrib", "contribution", "cos2",

  # Other dimensionality reduction methods
  "UMAP1", "UMAP2", "umap1", "umap2",
  "tSNE1", "tSNE2", "tsne1", "tsne2", "tsne.1", "tsne.2",
  "MDS1", "MDS2", "mds1", "mds2",
  "component", "factor", "loading", "rotation",

  # Distance and clustering
  "distance", "dist", "similarity", "dissimilarity",
  "cluster", "clust", "k", "silhouette", "within.ss", "between.ss",

  # -------------------------------------------------------------------------
  # Plotting and Visualization Variables
  # -------------------------------------------------------------------------

  # Aesthetic mappings
  "fill", "stroke", "width", "height", "angle", "hjust", "vjust",
  "family", "fontface", "linewidth", "linetype", "pointsize",

  # Plot elements and coordinates
  "xmin", "xmax", "ymin", "ymax", "xend", "yend",
  "text", "annotation", "tooltip", "hover",

  # Faceting and grouping
  "facet", "panel", "wrap", "grid", "scales",

  # Color and theme variables
  "palette", "scale", "guide", "legend", "axis",

  # -------------------------------------------------------------------------
  # Quality Control and Filtering
  # -------------------------------------------------------------------------

  # Outlier detection
  "outlier", "outliers", "is.outlier", "outlier.status",
  "lower.fence", "upper.fence", "iqr", "mad", "zscore", "z.score",

  # Quality metrics
  "quality", "qc", "pass", "fail", "filter", "filtered",
  "threshold", "cutoff", "limit", "boundary",

  # Missing data
  "missing", "na", "complete", "incomplete", "coverage",

  # -------------------------------------------------------------------------
  # Time Series and Longitudinal Analysis
  # -------------------------------------------------------------------------

  "time", "timepoint", "visit", "day", "week", "month", "year",
  "baseline", "followup", "change", "delta", "slope", "trend",

  # -------------------------------------------------------------------------
  # Specialized Analysis Variables
  # -------------------------------------------------------------------------

  # Differential abundance/expression
  "DA", "DE", "DAM", "DEG", "differential", "significant.features",
  "log2fc.threshold", "padj.threshold", "lfc.threshold",

  # Enrichment analysis
  "enrichment.score", "nes", "fdr.qval", "nominal.pvalue",
  "leading.edge", "core.enrichment",

  # Network analysis
  "node", "edge", "degree", "betweenness", "closeness", "modularity",
  "hub", "cluster.id", "community",

  # Survival analysis
  "time.to.event", "event", "censored", "survival", "hazard.ratio",
  "risk.score", "risk.group",

  # -------------------------------------------------------------------------
  # File and Data Import/Export
  # -------------------------------------------------------------------------

  "file", "path", "filename", "extension", "header", "separator",
  "sheet", "tab", "row.names", "col.names",

  # -------------------------------------------------------------------------
  # Miscellaneous Variables
  # -------------------------------------------------------------------------

  # Logical and status indicators
  "TRUE", "FALSE", "yes", "no", "positive", "negative", "neutral",
  "present", "absent", "detected", "not.detected",

  # Iteration and indexing
  "i", "j", "k", "iter", "iteration", "step", "round",

  # Temporary and intermediate variables
  "temp", "tmp", "result", "output", "summary", "info",
  "meta", "metadata", "data.new", "df.temp",

  # Units and measurements
  "unit", "units", "scale.factor", "normalization.factor",
  "concentration", "molarity", "volume", "weight", "area"
))

# =============================================================================
# Additional Package-Specific Variables
# =============================================================================

# Variables specific to bioRtools functions that may be added dynamically
# These are declared to prevent future R CMD check warnings

utils::globalVariables(c(
  # Function-specific temporary variables
  "group.anova", "value.temp", "data.temp", "result.temp",
  "term2gene", "term2name", "gene2term", "name2term",

  # Analysis result columns that may be generated
  "rank", "ranking", "score", "weight", "importance", "priority",
  "magnitude", "direction", "trend.direction", "regulation.direction",

  # Comparison and contrast variables
  "contrast", "comparison", "versus", "vs", "compared.to", "reference.group",
  "fold.change.direction", "change.direction", "effect.direction",

  # Advanced statistical variables
  "leverage", "influence", "cooks.distance", "hat.values",
  "standardized.residuals", "studentized.residuals",

  # Bootstrap and resampling
  "bootstrap", "resample", "permutation", "iteration.number",
  "confidence.interval", "prediction.interval",

  # Cross-validation and model evaluation
  "fold", "cv.fold", "train", "test", "validation",
  "accuracy", "precision", "recall", "f1.score", "auc", "roc"
))

# =============================================================================
# Documentation Note
# =============================================================================

#' @note
#' This file serves several important purposes:
#'
#' 1. **R CMD check compliance**: Prevents "no visible binding for global variable"
#'    warnings during package checking
#'
#' 2. **Code clarity**: Documents all variables used in non-standard evaluation
#'    contexts throughout the package
#'
#' 3. **Maintenance**: Provides a centralized location for managing global
#'    variable declarations
#'
#' 4. **Extensibility**: Makes it easy to add new variables as the package grows
#'
#' Variables are organized by functional categories to improve readability and
#' maintenance. When adding new functions that use NSE, please add any new
#' global variables to the appropriate section above.
