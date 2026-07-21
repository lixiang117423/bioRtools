#' Test data.
#'
#' @format ## `test_data`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{Sepal.Length}{length of sepal}
#'   \item{Sepal.Width}{width of sepal}
#'   \item{Petal.Length}{length of petal}
#'   \item{Petal.Width}{width of petal}
#'   \item{Species}{species group}
#' }
#' @source <https://archive.ics.uci.edu/dataset/53/iris>
"test_data"

#' Test data.
#'
#' @format ## `df.reorder2heatmap`
#' A data frame with 4,270 rows and 3 columns:
#' \describe{
#'   \item{sample}{sample code}
#'   \item{meta}{metabolite name}
#'   \item{value}{metabolite abundance}
#' }
#' @source Xiang Li
"df.reorder2heatmap"

#' Test data.
#'
#' @format ## `df.rnaseq.gene`
#' A data frame with 57,359 rows and 6 columns:
#' \describe{
#'   \item{SRR12580270}{sample name}
#'   \item{SRR12580269}{sample name}
#'   \item{SRR12580268}{sample name}
#'   \item{SRR12580258}{sample name}
#'   \item{SRR12580257}{sample name}
#'   \item{SRR12580256}{sample name}
#' }
#' @source Xiang Li
"df.rnaseq.gene"

#' Test data.
#'
#' @format ## `df.rnaseq.sample`
#' A data frame with 6 rows and 1 columns:
#' \describe{
#'   \item{group}{sample group}
#' }
#' @source Xiang Li
"df.rnaseq.sample"

#' Test data.
#'
#' @format ## `df.rnaseq.go`
#' A data frame with 6 rows and 1 columns:
#' \describe{
#'   \item{gene}{gene name}
#'   \item{go_id}{GO term id}
#'   \item{go_term}{GO term name}
#'   \item{go_ontology}{GO term ontology}
#' }
#' @source Xiang Li
"df.rnaseq.go"

#' Test data.
#'
#' @format ## `df.rnaseq.degs`
#' A data frame with 6 rows and 1 columns:
#' \describe{
#'   \item{gene}{gene name}
#' }
#' @source Xiang Li
"df.rnaseq.degs"

#' Test data.
#'
#' @format ## `df.rnaseq.kegg`
#' A data frame with 6 rows and 1 columns:
#' \describe{
#'   \item{gene}{gene name}
#'   \item{kegg_id}{KEGG term id}
#'   \item{kegg_term}{KEGG term name}
#' }
#' @source Xiang Li
"df.rnaseq.kegg"

#' Test data.
#'
#' @format ## `df.rnaseq.plot.volcano`
#' A data frame with 57,359 rows and 8 columns:
#' \describe{
#'   \item{gene}{gene name}
#'   \item{baseMea}{baseMea}
#'   \item{log2FoldChange}{log2FoldChange}
#'   \item{lfcSE}{lfcSE}
#'   \item{stat}{stat}
#'   \item{pvalue}{pvalue}
#'   \item{padj}{padj}
#'   \item{group}{group}
#' }
#' @source Xiang Li
"df.rnaseq.plot.volcano"

#' Test data.
#'
#' @format ## `df.pcoa.otu`
#' A data frame with 9 rows and 11,267 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.pcoa.otu"

#' Test data.
#'
#' @format ## `df.pcoa.sample`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#'   \item{sample}{sample name}
#'   \item{group}{sample group}
#' }
#' @source Xiang Li
"df.pcoa.sample"

#' Test data.
#'
#' @format ## `df.top10.otu`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.top10.otu"

#' Test data.
#'
#' @format ## `df.top10.sample`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#'   \item{sample}{sample name}
#'   \item{group}{sample group}
#' }
#' @source Xiang Li
"df.top10.sample"

#' Test data.
#'
#' @format ## `df.top10.class`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.top10.class"

#' Test data.
#'
#' @format ## `df.call_dams_lefse.otu`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.call_dams_lefse.otu"

#' Test data.
#'
#' @format ## `df.call_dams_lefse.sample`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.call_dams_lefse.sample"

#' Test data.
#'
#' @format ## `df.rda.chem`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.rda.chem"

#' Test data.
#'
#' @format ## `df.rda.otu`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.rda.otu"

#' Test data.
#'
#' @format ## `df.permanova.sample`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.permanova.sample"

#' Test data.
#'
#' @format ## `df.permanova.otu`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.permanova.otu"

#' Test data.
#'
#' @format ## `df.splsda.meta`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.splsda.meta"

#' Test data.
#'
#' @format ## `df.splsda.sample`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.splsda.sample"

#' Test data for synteny plot.
#'
#' @name df.synteny.gene
#' @title Gene Synteny Example Data
#' @usage df.synteny.gene
#' @format A data frame with 15 rows and 7 columns containing gene information
#'   for synteny visualization:
#' \describe{
#'   \item{start}{Gene start position (numeric)}
#'   \item{end}{Gene end position (numeric)}
#'   \item{y}{Y-axis position for arranging genes on different tracks (numeric)}
#'   \item{gene_name}{Gene identifier (character)}
#'   \item{species}{Species name (character)}
#'   \item{direction}{Gene strand direction, "+" or "-" (character)}
#'   \item{gene_type}{Gene category for color coding (character)}
#' }
#' @source Simulated example data inspired by comparative genomics studies
#' @description{
#' Test data for synteny plot containing gene information for 4 species
#' with 15 genes total.
#' }
#' @examples
#' library(bioRtools)
#' data(df.synteny.gene)
#' head(df.synteny.gene)
"df.synteny.gene"

#' Test data for synteny plot links (curve style).
#'
#' @name df.synteny.link
#' @title Synteny Link Example Data (Curve Style)
#' @usage df.synteny.link
#' @format A data frame with 14 rows and 3 columns containing synteny relationship
#'   information for connecting homologous genes with simple curves:
#' \describe{
#'   \item{x}{X-coordinate for the link connection point (numeric)}
#'   \item{y}{Y-coordinate for the link connection point (numeric)}
#'   \item{group}{Grouping variable that identifies which genes should
#'     be connected across species (numeric/character)}
#' }
#' @details
#' Each group has exactly 2 rows representing the start and end points of a curve.
#' Use this data with link_type = "curve" in plot_synteny().
#' @source Simulated example data demonstrating synteny relationships
#' @description{
#' Test data for synteny plot containing 7 synteny links connecting
#' homologous genes across 4 species. Each link has 2 points (start and end).
#' }
#' @examples
#' data(df.synteny.link)
#' head(df.synteny.link)
"df.synteny.link"

#' Test data for synteny plot links (ribbon style).
#'
#' @name df.synteny.link.ribbon
#' @title Synteny Link Example Data (Ribbon Style)
#' @usage df.synteny.link.ribbon
#' @format A data frame with 28 rows and 3 columns containing synteny relationship
#'   information for connecting homologous genes with filled ribbons:
#' \describe{
#'   \item{x}{X-coordinate for the link region corner (numeric)}
#'   \item{y}{Y-coordinate for the link region corner (numeric)}
#'   \item{group}{Grouping variable that identifies which genes should
#'     be connected across species (numeric/character)}
#' }
#' @details
#' Each group has exactly 4 rows representing the corners of a ribbon region
#' connecting two gene intervals. Use this data with link_type = "ribbon" in plot_synteny().
#' @source Simulated example data demonstrating synteny relationships
#' @description{
#' Test data for synteny plot containing 7 synteny links connecting
#' homologous genes across 4 species. Each link has 4 points defining
#' the corners of a ribbon-style connection.
#' }
#' @examples
#' data(df.synteny.link.ribbon)
#' head(df.synteny.link.ribbon)
"df.synteny.link.ribbon"

#' Human karyotype (test data for plot_ideogram)
#'
#' Human chromosomes with centromere positions, from the RIdeogram package.
#'
#' @format ## `df.ideo.karyotype`
#' A data frame with 24 rows and 5 columns:
#' \describe{
#'   \item{Chr}{chromosome identifier}
#'   \item{Start}{chromosome start (bp)}
#'   \item{End}{chromosome length / end (bp)}
#'   \item{CE_start}{centromere start (bp)}
#'   \item{CE_end}{centromere end (bp)}
#' }
#' @source RIdeogram (Hao et al. 2020, PeerJ Comput. Sci. 6:e251), Artistic-2.0.
"df.ideo.karyotype"

#' Gene density per 1 Mb window (test data for plot_ideogram overlaid heatmap)
#'
#' @format ## `df.ideo.gene_density`
#' A data frame with 3102 rows and 4 columns:
#' \describe{
#'   \item{Chr}{chromosome identifier}
#'   \item{Start}{window start (bp)}
#'   \item{End}{window end (bp)}
#'   \item{Value}{gene count in the window}
#' }
#' @source RIdeogram (Hao et al. 2020), Artistic-2.0.
"df.ideo.gene_density"

#' LTR density per 1 Mb window (test data for plot_ideogram heatmap label track)
#'
#' @format ## `df.ideo.ltr_density`
#' A data frame with 3102 rows and 4 columns:
#' \describe{
#'   \item{Chr}{chromosome identifier}
#'   \item{Start}{window start (bp)}
#'   \item{End}{window end (bp)}
#'   \item{Value}{LTR count in the window}
#' }
#' @source RIdeogram (Hao et al. 2020), Artistic-2.0.
"df.ideo.ltr_density"

#' Random RNA gene markers (test data for plot_ideogram marker track)
#'
#' @format ## `df.ideo.rna_marker`
#' A data frame with 500 rows and 6 columns:
#' \describe{
#'   \item{Type}{marker category (legend key)}
#'   \item{Shape}{one of \code{triangle}, \code{box}, \code{circle}}
#'   \item{Chr}{chromosome identifier}
#'   \item{Start}{feature start (bp)}
#'   \item{End}{feature end (bp)}
#'   \item{color}{fill colour, hex without \code{#}}
#' }
#' @source RIdeogram (\code{Random_RNAs_500}; Hao et al. 2020), Artistic-2.0.
"df.ideo.rna_marker"

#' Dual-species karyotype (test data for plot_ideogram_synteny)
#'
#' Chromosomes of \emph{Vitis vinifera} (Grape) and \emph{Populus} compared in
#' a synteny idiogram.
#'
#' @format ## `df.ideo.synteny_karyotype`
#' A data frame with 38 rows and 7 columns:
#' \describe{
#'   \item{Chr}{chromosome identifier}
#'   \item{Start}{start (bp)}
#'   \item{End}{chromosome length (bp)}
#'   \item{fill}{chromosome body colour, hex without \code{#}}
#'   \item{species}{species name}
#'   \item{size}{species-label font size}
#'   \item{color}{species-label colour, hex without \code{#}}
#' }
#' @source RIdeogram (\code{karyotype_dual_comparison}; Hao et al. 2020), Artistic-2.0.
"df.ideo.synteny_karyotype"

#' Dual-species synteny links (test data for plot_ideogram_synteny)
#'
#' @format ## `df.ideo.synteny`
#' A data frame with 2483 rows and 7 columns:
#' \describe{
#'   \item{Species_1}{1-based chromosome index in the first species}
#'   \item{Start_1}{syntenic block start on species 1 (bp)}
#'   \item{End_1}{syntenic block end on species 1 (bp)}
#'   \item{Species_2}{1-based chromosome index in the second species}
#'   \item{Start_2}{syntenic block start on species 2 (bp)}
#'   \item{End_2}{syntenic block end on species 2 (bp)}
#'   \item{fill}{ribbon colour, hex without \code{#}}
#' }
#' @source RIdeogram (\code{synteny_dual_comparison}; Hao et al. 2020), Artistic-2.0.
"df.ideo.synteny"
