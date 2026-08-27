#' @keywords internal
"_PACKAGE"

# NSE column names / globals used by dplyr/ggplot2 expressions - declared so
# R CMD check does not flag them as "no visible binding for global variable".
utils::globalVariables(c(
  "..count..",   "..scaled..",   ".pair",   "abs_log2fc",   "abs_pos",   "anova_signif",
  "beta_mntd",   "beta_nti",   "bg",   "BgRatio",   "bin",   "BioRep",
  "chr_color",   "chr_factor",   "ci_lower",   "ci_upper",   "col_var",   "Correlation",
  "Cq",   "cum_pos_end",   "cum_pos_start",   "Description",   "dist_bin_kb",   "dist_kb",
  "domain",   "Eff",   "enrichment_score",   "eq.label",   "gene_id",   "gene_length",
  "gene_mean",   "gene_ratio",   "gene_sd",   "geneID",   "group_anova",   "group_biorep",
  "group.x",   "group.y",   "groups",   "half_decay_kb",   "ID",   "indicator",
  "isTip",   "max_coord",   "max_pos",   "mean_abundance",   "mean_connectivity",   "mean_r2",
  "mean_r2_bin",   "mean_r2_fitted",   "mean_sample",   "mean_val",   "mean_value",   "mean.k.",
  "metric",   "min_coord",   "module",   "Module",   "module_color",   "module_number",
  "N",   "n_obs",   "n_samples",   "n_snp",   "neg_log10_p",   "no visible binding for global variable",
  "overall_pvalue",   "overall_pvalue_formatted",   "p_adj",   "p_adjust",   "p_value",   "path_length",
  "path_string",   "path_to_root",   "pathway_size",   "PC",   "pd",   "peak_count",
  "percentage",   "percentage_label",   "phenotype",   "plot_label",   "pos",   "pos_cum",
  "Power",   "Pr(>F)",   "process",   "proportion",   "r",   "r_squared",
  "r2",   "r2_fit",   "rd",   "RDA1_scaled",   "RDA2_scaled",   "recommendation",
  "rect_id",   "relative_abundance",   "relative_eig",   "row_order",   "row_sum",   "row_var",
  "rr.label",   "run",   "sample_id",   "scale_free_r2",   "sd_abundance",   "SFT.R.sq",
  "stats",   "sum_value",   "tax",   "tax_group",   "taxon",   "TechRep",
  "term_size",   "trait",   "tukey_signif",   "V3",   "V4",   "V5",
  "V7",   "val",   "value_var",   "variable",   "variance_percent",   "vip",
  "x_end",   "x_start",   "x_var",   "y_end",   "y_max",   "y_min",
  "y_num",   "y_start",   "y_var"
))





## usethis namespace: start
#' @importFrom agricolae duncan.test
#' @importFrom ape pcoa
#' @importFrom ape read.tree
#' @importFrom broom tidy
#' @importFrom clusterProfiler enricher
#' @importFrom data.table fread
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 results
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr count
#' @importFrom dplyr inner_join
#' @importFrom dplyr lag
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr slice_max
#' @importFrom factoextra get_eigenvalue
#' @importFrom FactoMineR PCA
#' @importFrom forcats as_factor
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 arrow
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 position_fill
#' @importFrom ggplot2 rel
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 unit
#' @importFrom ggprism guide_prism_offset_minor
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggsci scale_color_d3
#' @importFrom ggtext element_markdown
#' @importFrom ggthemes theme_foundation
#' @importFrom ggtree geom_tiplab
#' @importFrom ggtree ggtree
#' @importFrom janitor clean_names
#' @importFrom lefser lefser
#' @importFrom magrittr set_names
#' @importFrom mixOmics splsda
#' @importFrom multcomp cld
#' @importFrom multcomp glht
#' @importFrom multcomp mcp
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom rlang sym
#' @importFrom ropls opls
#' @importFrom rstatix wilcox_test
#' @importFrom scales percent
#' @importFrom scales pretty_breaks
#' @importFrom stats anova
#' @importFrom stats aov
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats IQR
#' @importFrom stats pf
#' @importFrom stats quantile
#' @importFrom stats reorder
#' @importFrom stringr str_replace
#' @importFrom stringr str_split
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom vegan adonis
#' @importFrom vegan decostand
#' @importFrom vegan envfit
#' @importFrom vegan rda
#' @importFrom vegan scores
#' @importFrom vegan vegdist
#' @importFrom WGCNA corAndPvalue
## usethis namespace: end
NULL
