---
title: "bioRtools: reproducible analysis-to-figure workflows for omics data"
tags:
  - R
  - bioinformatics
  - microbiome
  - metabolomics
  - multivariate analysis
  - reproducible research
  - data visualization
authors:
  - name: Xiang LI
    affiliation: "1"
    corresponding: true
affiliations:
  - name: Yazhouwan National Laboratory, No. 8 Huanjin Road, Yazhou District, Sanya, Hainan 572024, China
    index: 1
date: 30 October 2026
bibliography: paper.bib
---

# Summary

Modern biological studies often measure thousands of genes, metabolites, or
microbial taxa in each sample. Turning those measurements into defensible
results requires several linked steps: fitting a statistical model, checking
the fitted object, extracting comparable tables, and producing figures that can
be revised for publication. The individual statistical methods are well
implemented in mature R packages, but their inputs, object structures, and
output conventions differ. Analysts therefore spend substantial effort moving
results between packages, matching sample metadata, renaming columns, and
redrawing similar figures.

`bioRtools` [@biortools2026] is an R package [@RCore2026] that provides a shared
analysis-to-figure layer for these workflows. It calls established engines such
as DESeq2 [@love2014deseq2], edgeR [@robinson2010edger], vegan
[@oksanen2025vegan], mixOmics [@rohart2017mixomics], and ropls
[@thevenot2015ropls] rather than reimplementing their statistical methods. Each
core analysis returns the native model together with named result tables,
editable `ggplot2` figures [@wickham2016ggplot2], and provenance describing the
input dimensions, engine versions, function call, and seed. The package covers
differential expression and abundance analysis, multivariate ordination and
classification, diversity analysis, enrichment, and common genomics figures.
It is designed for scripted, local, and reproducible analysis rather than as a
replacement for the underlying statistical engines.

# Statement of need

The target users are R-using researchers, students, and bioinformatics support
staff who repeatedly analyze plant, microbial, or other omics datasets. In
practice, a single project may combine DESeq2 for count data, vegan for
community ecology, mixOmics or ropls for multivariate metabolomics, and
`ggplot2` for final figures. Each package uses a different object model and
different conventions for sample identifiers, feature identifiers, missing
values, and result extraction. A seemingly simple error, such as metadata
sorted in a different order from the response matrix, can silently change a
test result.

Existing packages solve important parts of this problem, but they leave a gap
between a fitted model and a reusable analysis product. `bioRtools` addresses
that gap by enforcing a common result contract, checking sample alignment
before model fitting, preserving the native statistical object, and returning
the same table and plot vocabulary across methods. The contract is intended to
make downstream code predictable without hiding the statistical engine. A
thin syntax wrapper would not be sufficient: the research value lies in the
combination of identity checks, stable table semantics, provenance, and
composable figures across methods.

# State of the field

Several research-software projects make biological analysis more accessible.
BLIT integrates command-line bioinformatics tools into R workflows
[@ding2026blit]. cOmicsArt provides a customizable environment for omics
analysis and reporting [@seep2025comicsart]. TBtools offers an interactive
toolkit for large biological datasets [@chen2020tbtools]. Hiplot and ImageGP
provide web-based services for publication-ready biomedical and scientific
visualization [@li2022hiplot; @chen2022imagegp].

These tools differ from `bioRtools` in architecture and intended integration.
Web services and graphical toolkits are convenient for interactive use but are
not primarily R libraries for scripted model fitting. BLIT connects command-line
programs to R but does not define a shared result contract across statistical
models. Contributing a small feature to one of these projects would not provide
the cross-engine behavior needed here. `bioRtools` was therefore built as an
R-native layer that delegates estimation to established engines, keeps their
objects available for inspection, and standardizes the surrounding analysis
and figure interface. This complements rather than competes with the underlying
methods and visualization platforms.

# Software design

The central design decision is to separate three outputs that are often mixed
inside ad hoc scripts. A `bioRtools_fit` object contains `model`, the untouched
object returned by the statistical engine; `tables`, a named list of
analysis-ready data frames; `plots`, a named list of editable `ggplot2`
objects; and `provenance`, which records the analysis type, engine and package
versions, R version, input dimensions, originating call, and optional seed.
Historical list elements remain available for backward compatibility, while
`print()`, `summary()`, and `as_tibble()` provide predictable inspection and
extraction. The package is released under the MIT license.

Sample identity is treated as part of the statistical model, not as a
formatting detail. The multivariate workflows resolve explicit sample-ID
columns or non-default row names, require complete and unique identifiers, and
reorder metadata to the response before fitting. Missing, extra, or duplicated
identifiers produce an error instead of a silent positional match. This
alignment layer is shared by principal component analysis, principal coordinate
analysis, redundancy analysis, PERMANOVA, sparse PLS-DA, and OPLS-DA.

The package also separates core dependencies from optional engines. Core
workflows use a comparatively small set of imports; packages needed for
specific methods or visualizations are checked at call time and listed as
suggested dependencies. The trade-off is that some functions require an
additional installation, but users are not forced to install every supported
omics engine. The public vignettes exercise microbiome and metabolomics
workflows with `mia::GlobalPatterns` [@caporaso2011globalpatterns] and
`mixOmics::nutrimouse` [@martin2007nutrimouse], respectively. Tests compare
outputs with the native engines and cover metadata shuffling, duplicate and
missing identifiers, adjusted-P preservation, and end-to-end public-data
workflows.

# Research impact statement

`bioRtools` has been publicly developed on GitHub since 28 July 2025. At the
time of submission, the repository contains more than 300 commits distributed
over 15 months, more than 15 version tags, a maintained changelog, automated tests,
continuous-integration checks, reference documentation, and public vignettes.
This history shows iterative development rather than a one-time release.

The software is also used in the author's active research program. A local
audit found 944 R, R Markdown, and Quarto files across 22 research projects
calling exported `bioRtools` functions. These projects include soybean pathogen
pangenome and effector analyses, soybean GWAS and NLR work, peanut BSA and
RNA-seq analyses, rice metagenome analysis, and assembly-benchmarking studies.
That use has directly motivated the sample-alignment checks, preservation of
unavailable adjusted P values, optional-dependency guards, and reproducible
result contract described above. External peer-reviewed citations have not yet
accrued, so this submission relies on documented developer use and reproducible
public-data examples rather than claiming a mature external user community.

# AI usage disclosure

Generative AI assistance was used during software development, documentation,
and manuscript preparation. Git commit trailers record Anthropic Claude
assistance on 394 commits between 9 February and 29 August 2026, principally
for code generation and refactoring, test scaffolding, documentation, and
debugging; the exact Claude model version was not recorded. OpenAI Codex, based
on the GPT-5 family, was used in sessions on 28-30 October 2026 to audit the
package, reduce hard dependencies, repair tests and package APIs, develop the
shared result contract, and draft and edit this paper.

The human author reviewed and validated all AI-assisted outputs, ran the test
and package-check suites, and retained responsibility for problem framing,
statistical scope, API and architecture decisions, and final acceptance. All
changes were version-controlled and tested before release. No AI system is
listed as an author.

# Acknowledgements

The author thanks the R and Bioconductor communities and the maintainers of the
statistical engines integrated by `bioRtools`. No specific funding was received
for this work. The author declares no competing interests.

# References
