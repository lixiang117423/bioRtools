utils::globalVariables(c(
  # pca_in_one function
  "%>%",
  ".",
  "variance.percent",
  "pc",
  "plot.pca",
  "sym",
  "return.list",
  "point.data",
  # cor_and_plot function
  "result.cor",
  "df.cor",
  "df.pvalue",
  "from",
  "to",
  "pvalue",
  "cor",
  # lm_and_plot function
  "df.lm",
  "plot.lm",
  "result.lm",
  "anova",
  "model",
  "R2",
  # plot_theme
  "base_size",
  "base_family",
  "mytheme",
  # reorder2heatmap
  "my_data_peak_values",
  "number_of_peaks",
  "my_data_peak_values_reordered",
  "my_data_reordered",
  "peaked_at",
  "reorder",
  "order_rows",
  "n",
  "data",

  # call_DEGs_DESeq2
  "log2FoldChange",
  "padj",
  "degs",

  # enrich_GO
  "go.id",
  "go.term",
  "go.rich",
  "df.term2gene",
  "df.term2name",
  # enrich_KEGGG
  "kegg.id",
  "kegg.term",

  # PCoA
  "pcoa.res",
  "pcoa.weight",
  "df.pcoa.point",
  "pcoa.weight.plot",
  "df.x.label",
  "df.y.label",
  "p",

  # top_10
  "df.all",
  "df.top9",
  "df.all.new",
  "df.stat.new",
  "signif",
  "OTU",
  "Relative_eig",
  "group",
  "group2",
  "label",
  "position_fill",
  "relative.eig",
  "value",

  # RDA
  "data.new",
  "rda.res",
  "df.chem.point",
  "df.sample.point",
  "RDA1",
  "RDA2",
  "envfit.res",
  "chem",
  "PCoA",

  # Anova
  "data.new",
  "fit",
  "res",
  "anova.res",
  "anova.pvalue",
  "value",
  "group",
  "group.anova",
  "id",

  # SPLSDA
  "splsda.result",
  "splsda.result.point",
  "splsda.result.prop_expl_var",
  "comp",
  "prop",

  # OPLSDA
  "meta",
  "VIP",
  "var",

  # GWAS
  "df.chr.posi",
  "df.gwas.new",
  "chr",
  "posi",
  "max_posi",
  "posi_add",
  "posi_cum"
))
