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
#'   \item{go.id}{GO term id}
#'   \item{go.term}{GP term name}
#'   \item{go.ontology}{GP term  ontology}
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
#'   \item{kegg.id}{KEGG term id}
#'   \item{kegg.term}{KEGG term name}
#' }
#' @source Xiang Li
"df.rnaseq.kegg"

#' Test data.
#'
#' @format ## `df.rnaseq.plot_volcano`
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
"df.rnaseq.plot_volcano"

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
#' @format ## `df.top10.calss`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.top10.class"

#' Test data.
#'
#' @format ## `df.call_DAMs_LEfSe.otu`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.call_DAMs_LEfSe.otu"

#' Test data.
#'
#' @format ## `df.call_DAMs_LEfSe.sample`
#' A data frame with 9 rows and 2 columns:
#' \describe{
#' }
#' @source Xiang Li
"df.call_DAMs_LEfSe.sample"

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
