# Migrated from tests/test_synteny.R (bare script).
# plot_synteny() with the package example data. df.synteny.link carries
# curve-style links (2 rows per group), so link_type = "curve" is used.

df.synteny.gene <- bioRtools::df.synteny.gene
df.synteny.link <- bioRtools::df.synteny.link

test_that("basic usage returns the synteny plot and its summary", {
  result <- plot_synteny(
    gene_data = df.synteny.gene,
    syntenic_data = df.synteny.link,
    link_type = "curve"
  )

  expect_true(ggplot2::is_ggplot(result$plot_synteny))
  expect_s3_class(result$data_summary, "data.frame")
  expect_true(all(c("species", "n_genes", "total_links") %in%
                    names(result$data_summary)))
  expect_equal(sum(result$data_summary$n_genes), nrow(df.synteny.gene))
  expect_equal(unique(result$data_summary$total_links),
               length(unique(df.synteny.link$group)))
})

test_that("custom species labels replace the y-axis labels", {
  result <- plot_synteny(
    gene_data = df.synteny.gene,
    syntenic_data = df.synteny.link,
    link_type = "curve",
    species_labels = c("Species A", "Species B", "Species C", "Species D")
  )

  y_labels <- result$plot_synteny$scales$get_scales("y")$labels
  expect_equal(as.character(y_labels),
               c("Species A", "Species B", "Species C", "Species D"))
})

test_that("custom fill colors reach the gene fill scale", {
  result <- plot_synteny(
    gene_data = df.synteny.gene,
    syntenic_data = df.synteny.link,
    link_type = "curve",
    species_labels = c("A", "B", "C", "D"),
    fill_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
  )

  palette <- result$plot_synteny$scales$get_scales("fill")$palette(4)
  expect_equal(palette, c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))
})

test_that("title and subtitle are attached to the plot", {
  result <- plot_synteny(
    gene_data = df.synteny.gene,
    syntenic_data = df.synteny.link,
    link_type = "curve",
    species_labels = c("A", "B", "C", "D"),
    title = "Gene Synteny Visualization",
    subtitle = "Example from bioRtools package"
  )

  expect_identical(result$plot_synteny$labels$title,
                   "Gene Synteny Visualization")
  expect_identical(result$plot_synteny$labels$subtitle,
                   "Example from bioRtools package")
})

test_that("show_gene_labels = FALSE drops the gene-name text layer", {
  result <- plot_synteny(
    gene_data = df.synteny.gene,
    syntenic_data = df.synteny.link,
    link_type = "curve",
    species_labels = c("A", "B", "C", "D"),
    show_gene_labels = FALSE
  )

  layer_classes <- vapply(result$plot_synteny$layers,
                          function(l) class(l$geom)[1], character(1))
  expect_false("GeomText" %in% layer_classes)
})

test_that("link_alpha controls the curve link transparency", {
  result <- plot_synteny(
    gene_data = df.synteny.gene,
    syntenic_data = df.synteny.link,
    link_type = "curve",
    species_labels = c("A", "B", "C", "D"),
    link_alpha = 0.8
  )

  curve_layer <- result$plot_synteny$layers[[2]]
  expect_equal(curve_layer$aes_params$alpha, 0.8)
})
