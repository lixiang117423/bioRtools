#' Analyze Admixture Results with a Phylogenetic Tree
#'
#' @description
#' Process Admixture .Q files and create population structure plots with
#' integrated phylogenetic tree visualization using ggplot2 and ggtree. This
#' function creates a rotated layout with samples on the Y-axis and K values
#' arranged horizontally, and optionally combines the plot with a
#' phylogenetic tree.
#'
#' @param admixture_path Character string specifying the directory containing
#'   Admixture .Q files and the .fam file.
#' @param tree_file Optional path to a Newick tree file, or a \code{phylo}
#'   object (e.g. from \code{ape::read.tree()} or \code{\link{read_data}}).
#'   If provided, the tree is drawn and combined with the admixture plot.
#'   If \code{NULL}, only the admixture plot is returned. Default is \code{NULL}.
#' @param sample_group Optional data frame with columns \code{label} (sample
#'   names), \code{group}, and optionally \code{color}, or a file path to a
#'   tab-delimited file with the same columns. When provided together with a
#'   tree, tree tips are colored by group. Default is \code{NULL}.
#' @param k_range Integer vector of K values to include (e.g. \code{2:5} or
#'   \code{c(2, 3, 5)}). If \code{NULL}, all available K values are used.
#'   Default is \code{NULL}.
#' @param fam_id_col Column of the .fam file holding the sample ID
#'   (default: 2, the PLINK IID column). Use 1 to take the FID instead.
#' @param color_scale Color scale for the admixture bars. One of
#'   \code{"default"} (optimized base palette), \code{"aaas"}, \code{"npg"},
#'   \code{"lancet"}, \code{"jco"}, \code{"ucscgb"}, \code{"uchicago"},
#'   \code{"simpsons"}, \code{"rickandmorty"}, or \code{"manual"}.
#'   Default is \code{"aaas"}.
#' @param tree_color_scale Color scale for tree tip groups. One of
#'   \code{"default"}, \code{"d3"}, \code{"aaas"}, \code{"npg"},
#'   \code{"lancet"}, \code{"jco"}, \code{"ucscgb"}, \code{"uchicago"},
#'   \code{"simpsons"}, or \code{"rickandmorty"}. Default is \code{"d3"}.
#' @param manual_colors Character vector of colors for the admixture scale;
#'   used only when \code{color_scale = "manual"}. Default is \code{NULL}.
#' @param tree_width Numeric relative width of the tree panel
#'   (0-1, default: 0.2).
#' @param tree_linewidth Numeric line width for tree branches (default: 0.2).
#' @param tip_point_size Numeric size of tree tip points (default: 0.1).
#' @param show_branch_length Logical; show branch lengths (default: FALSE).
#' @param ladderize Logical; reorder the tree to a ladderized layout,
#'   matching \code{ggtree}'s default (default: TRUE). The admixture bar
#'   order always follows the drawn tree.
#' @param save Optional output file path. When given, the main plot (the
#'   combined plot in tree mode, otherwise the admixture plot) is saved via
#'   \code{\link{ggsave2}} (8 x 6 in, 600 dpi; the directory is created
#'   when missing). Default is \code{NULL} (no file written).
#' @param verbose Logical; print progress messages (default: TRUE).
#'
#' @return An object of class \code{admixture_phylo_results} (a list) with:
#'   \item{combined_plot}{aplot object (NULL without tree)}
#'   \item{admixture_plot}{ggplot2 structure plot}
#'   \item{tree_plot}{ggtree object (NULL without tree)}
#'   \item{samples}{character vector of sample IDs from the .fam file}
#'   \item{k_values}{numeric vector of K values used}
#'   \item{matched_samples}{character vector of samples present in both the
#'   .Q files and the tree}
#'   A \code{print} method shows a summary.
#'
#' @details
#' Two modes exist depending on whether a tree is provided:
#'
#' \strong{Mode 1: without tree (traditional layout)}
#' \itemize{
#'   \item Samples on the X-axis, K values as vertical facets
#'   \item Strip labels on the left
#'   \item Returns the admixture plot only
#' }
#'
#' \strong{Mode 2: with tree (rotated layout)}
#' \itemize{
#'   \item Samples on the Y-axis in tree tip order, K values as horizontal
#'     facets
#'   \item Strip labels on top
#'   \item Optional tip coloring from \code{sample_group}
#'   \item Returns the combined plot, admixture plot, and tree plot
#' }
#'
#' The .fam file may be named arbitrarily but must end in \code{.fam}. By
#' default sample IDs are taken from column 2 (PLINK IID); pass
#' \code{fam_id_col = 1} to use the FID column. Each .Q file must have
#' exactly as many rows as the .fam file. Samples present in only one of the
#' tree and the .Q files are reported with a warning (not silently dropped).
#'
#' @note
#' Sample labels are matched after removing surrounding single quotes, so
#' \code{'S1'} in a tree matches \code{S1} in the .fam file. No files are
#' saved unless \code{save} is given.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # Mode 1: without tree (traditional layout)
#' result <- admixture_phylo_analysis("path/to/admixture/results/")
#' result$admixture_plot
#'
#' # Mode 2: with tree, tree read beforehand and reused
#' tree <- read_data("path/to/phylogenetic_tree.nwk")  # phylo object
#' result <- admixture_phylo_analysis(
#'   "path/to/admixture/results/",
#'   tree_file = tree,
#'   k_range = 2:6
#' )
#' result$combined_plot
#'
#' # Mode 2 with sample grouping
#' sample_group <- data.frame(
#'   label = c("S1", "S2", "S3"),
#'   group = c("A", "A", "B")
#' )
#' result <- admixture_phylo_analysis(
#'   "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   sample_group = sample_group,
#'   color_scale = "npg",
#'   tree_color_scale = "aaas",
#'   tree_width = 0.3,
#'   save = "figures/admixture_tree.png"
#' )
#' }
admixture_phylo_analysis <- function(admixture_path,
                                     tree_file = NULL,
                                     sample_group = NULL,
                                     k_range = NULL,
                                     fam_id_col = 2,
                                     color_scale = c("default", "aaas", "npg", "lancet", "jco", "ucscgb", "uchicago", "simpsons", "rickandmorty", "manual"),
                                     tree_color_scale = c("default", "d3", "aaas", "npg", "lancet", "jco", "ucscgb", "uchicago", "simpsons", "rickandmorty"),
                                     manual_colors = NULL,
                                     tree_width = 0.2,
                                     tree_linewidth = 0.2,
                                     tip_point_size = 0.1,
                                     show_branch_length = FALSE,
                                     ladderize = TRUE,
                                     save = NULL,
                                     verbose = TRUE) {
  # ---- 1. Input validation ----
  if (!dir.exists(admixture_path)) {
    stop("Admixture directory not found: ", admixture_path)
  }
  if (!is.numeric(fam_id_col) || length(fam_id_col) != 1L || is.na(fam_id_col) ||
      fam_id_col < 1 || fam_id_col != floor(fam_id_col)) {
    stop("'fam_id_col' must be a single positive integer")
  }
  if (!is.null(k_range) &&
      (!is.numeric(k_range) || length(k_range) == 0L || any(!is.finite(k_range)))) {
    stop("'k_range' must be NULL or a numeric vector of K values")
  }
  if (!is.numeric(tree_width) || length(tree_width) != 1L || is.na(tree_width) ||
      tree_width <= 0 || tree_width >= 1) {
    stop("'tree_width' must be a single number between 0 and 1")
  }
  if (!is.numeric(tree_linewidth) || length(tree_linewidth) != 1L ||
      is.na(tree_linewidth) || tree_linewidth <= 0) {
    stop("'tree_linewidth' must be a single positive number")
  }
  if (!is.numeric(tip_point_size) || length(tip_point_size) != 1L ||
      is.na(tip_point_size) || tip_point_size <= 0) {
    stop("'tip_point_size' must be a single positive number")
  }
  if (!is.logical(show_branch_length) || length(show_branch_length) != 1L ||
      is.na(show_branch_length)) {
    stop("'show_branch_length' must be TRUE or FALSE")
  }
  if (!is.logical(ladderize) || length(ladderize) != 1L || is.na(ladderize)) {
    stop("'ladderize' must be TRUE or FALSE")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE")
  }
  if (!is.null(save) &&
      (!is.character(save) || length(save) != 1L || !nzchar(save))) {
    stop("'save' must be NULL or a single non-empty file path")
  }
  if (!is.null(manual_colors) && !is.character(manual_colors)) {
    stop("'manual_colors' must be NULL or a character vector")
  }

  color_scale <- match.arg(color_scale)
  tree_color_scale <- match.arg(tree_color_scale)

  if (verbose) {
    message("Starting admixture phylogenetic analysis...")
    message("Admixture path: ", admixture_path)
    message("Tree file: ", if (is.null(tree_file)) "none (admixture plot only)" else "provided")
    message("Color scale: ", color_scale)
  }

  # ---- 2. Read the tree once and extract tip order ----
  phylo_tree <- NULL
  df_order_tree <- NULL

  if (!is.null(tree_file)) {
    if (inherits(tree_file, "phylo")) {
      phylo_tree <- tree_file
    } else if (is.character(tree_file) && length(tree_file) == 1L &&
               file.exists(tree_file)) {
      phylo_tree <- tryCatch(
        ape::read.tree(tree_file),
        error = function(e) stop("Failed to read tree file: ", rlang::cnd_message(e))
      )
    } else {
      stop("'tree_file' must be a phylo object or a path to an existing Newick file")
    }

    # ladderize here exactly like ggtree's default (right = FALSE), so the
    # extracted tip order matches the tree as it will be drawn
    if (ladderize) {
      phylo_tree <- ape::ladderize(phylo_tree, right = FALSE)
    }

    # clean tip labels once, so the order table and the drawn tree agree
    phylo_tree$tip.label <- .clean_tip_labels(phylo_tree$tip.label)
  }

  # ---- 3. Read the .fam file ----
  if (verbose) message("Reading .fam file...")

  fam_files <- list.files(admixture_path, pattern = "\\.fam$", full.names = TRUE)
  if (length(fam_files) == 0L) {
    stop("No .fam file found in: ", admixture_path)
  }
  if (length(fam_files) > 1L) {
    warning("Multiple .fam files found; using the first: ", fam_files[1])
  }

  df_sample <- tryCatch(
    {
      fam_raw <- readr::read_delim(fam_files[1], col_names = FALSE,
                                   show_col_types = FALSE,
                                   col_types = readr::cols(.default = readr::col_character()))
      if (fam_id_col > ncol(fam_raw)) {
        stop("'fam_id_col' (", fam_id_col, ") exceeds the .fam column count (", ncol(fam_raw), ")")
      }
      data.frame(sample = fam_raw[[fam_id_col]], stringsAsFactors = FALSE)
    },
    error = function(e) stop("Failed to read .fam file: ", rlang::cnd_message(e))
  )

  if (verbose) message("Read ", nrow(df_sample), " samples from .fam file")

  # default order (no tree) follows the .fam file
  if (is.null(df_order_tree)) {
    df_order_tree <- data.frame(
      sample = df_sample$sample,
      order = seq_len(nrow(df_sample)),
      stringsAsFactors = FALSE
    )
  }

  # ---- 4. Read .Q files ----
  if (verbose) message("Reading .Q files...")

  q_files <- list.files(admixture_path, pattern = "\\.Q$", full.names = TRUE)
  if (length(q_files) == 0L) {
    stop("No .Q files found in: ", admixture_path)
  }

  df_admixture <- tryCatch(
    purrr::map_dfr(q_files, function(f) {
      k_value <- as.numeric(stringr::str_extract(basename(f), "(?<=\\.)\\d+(?=\\.Q$)"))
      if (is.na(k_value)) {
        stop("Cannot extract K from file name: ", basename(f),
             " (expected a name like 'project.3.Q')")
      }
      q <- readr::read_delim(f, col_names = FALSE, show_col_types = FALSE)
      if (nrow(q) != nrow(df_sample)) {
        stop("Row mismatch in ", basename(f), ": ", nrow(q),
             " rows in .Q vs ", nrow(df_sample), " samples in .fam")
      }
      names(q) <- paste0("Cluster", seq_len(ncol(q)))
      q$sample <- df_sample$sample
      q %>%
        tidyr::pivot_longer(cols = dplyr::starts_with("Cluster"),
                            names_to = "name", values_to = "value") %>%
        dplyr::mutate(k_value = k_value, k = paste0("K = ", k_value))
    }),
    error = function(e) stop("Failed to read .Q files: ", rlang::cnd_message(e))
  )

  available_k <- sort(unique(df_admixture$k_value))
  if (verbose) message("Found K values: ", paste(available_k, collapse = ", "))

  # ---- 5. Filter K values ----
  if (!is.null(k_range)) {
    missing_k <- setdiff(k_range, available_k)
    if (length(missing_k) > 0L) {
      warning("Requested K value(s) not found in .Q files: ",
              paste(missing_k, collapse = ", "))
    }
    df_admixture <- df_admixture %>%
      dplyr::filter(.data$k_value %in% k_range)
    if (nrow(df_admixture) == 0L) {
      stop("No data found for specified K values: ", paste(k_range, collapse = ", "))
    }
  }
  available_k <- sort(unique(df_admixture$k_value))
  df_admixture <- df_admixture %>%
    dplyr::mutate(k = factor(.data$k, levels = paste0("K = ", available_k)))

  # ---- 5b. Build the tree plot first, so the bar order can follow it ----
  p_tree <- NULL
  if (!is.null(phylo_tree)) {
    if (verbose) message("Creating phylogenetic tree plot...")

    # the tree was already ladderized above, so drawing must not reorder again
    p_tree <- ggtree::ggtree(phylo_tree,
      branch.length = if (show_branch_length) NULL else "none",
      ladderize = FALSE,
      linewidth = tree_linewidth)

    # take the tip order from the drawn tree itself, bottom to top
    df_order_tree <- p_tree$data %>%
      dplyr::filter(.data$isTip) %>%
      dplyr::arrange(.data$y) %>%
      dplyr::transmute(sample = .data$label, order = .data$y) %>%
      as.data.frame()

    if (verbose) message("Found ", nrow(df_order_tree), " samples in tree")
  }

  # ---- 6. Check sample matching ----
  matched_samples <- intersect(df_admixture$sample, df_order_tree$sample)

  if (!is.null(phylo_tree)) {
    if (length(matched_samples) == 0L) {
      stop("No samples match between .Q files and phylogenetic tree.\n",
           "  .Q/.fam samples (first 10): ",
           paste(head(unique(df_admixture$sample), 10), collapse = ", "), "\n",
           "  Tree samples (first 10): ",
           paste(head(df_order_tree$sample, 10), collapse = ", "))
    }
    only_q <- setdiff(df_admixture$sample, df_order_tree$sample)
    only_tree <- setdiff(df_order_tree$sample, df_admixture$sample)
    if (length(only_q) > 0L) {
      warning(length(only_q), " sample(s) in .Q but not in the tree; dropped from the plot: ",
              paste(head(only_q, 10), collapse = ", "))
    }
    if (length(only_tree) > 0L) {
      warning(length(only_tree), " sample(s) in the tree but not in .Q; shown as blank rows: ",
              paste(head(only_tree, 10), collapse = ", "))
    }
    if (verbose) message("Successfully matched ", length(matched_samples), " samples with tree")
  } else if (verbose) {
    message("Using ", length(matched_samples), " samples in default order")
  }

  # ---- 7. Base color palette ----
  total_clusters <- length(unique(df_admixture$name))
  if (verbose) message("Total clusters needed: ", total_clusters)

  base_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
    "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
    "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5")

  # ---- 8. Admixture plot ----
  if (verbose) message("Creating admixture plot...")

  if (is.null(phylo_tree)) {
    # Mode 1: traditional layout (samples on X-axis, K facets stacked vertically)
    p_admixture <- df_admixture %>%
      ggplot2::ggplot(ggplot2::aes(.data$sample, .data$value, fill = .data$name)) +
      ggplot2::geom_bar(stat = "identity", width = 1,
                        position = ggplot2::position_fill()) +
      ggplot2::scale_x_discrete(limits = df_order_tree$sample) +
      ggplot2::facet_wrap(~ k, ncol = 1, strip.position = "left") +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.spacing.y = ggplot2::unit(0, "lines")
      )
  } else {
    # Mode 2: rotated layout (samples on Y-axis in tree order, K facets horizontal)
    p_admixture <- df_admixture %>%
      ggplot2::ggplot(ggplot2::aes(.data$value, .data$sample, fill = .data$name)) +
      ggplot2::geom_bar(stat = "identity", width = 1,
                        position = ggplot2::position_fill()) +
      ggplot2::scale_y_discrete(limits = df_order_tree$sample) +
      ggplot2::facet_wrap(~ k, ncol = length(available_k), strip.position = "top") +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_text(angle = 90),
        strip.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.spacing.x = ggplot2::unit(0, "lines")
      )
  }

  # ---- 9. Apply color scale ----
  if (color_scale == "default") {
    p_admixture <- p_admixture +
      ggplot2::scale_fill_manual(values = .extend_colors(total_clusters, base_colors))
  } else if (color_scale == "manual") {
    if (is.null(manual_colors)) {
      warning("manual_colors not provided, using default colors")
      p_admixture <- p_admixture +
        ggplot2::scale_fill_manual(values = .extend_colors(total_clusters, base_colors))
    } else if (length(manual_colors) < total_clusters) {
      warning("Not enough manual colors provided (", length(manual_colors),
        ") for clusters needed (", total_clusters, "). Extending with default colors.")
      p_admixture <- p_admixture +
        ggplot2::scale_fill_manual(
          values = c(manual_colors,
                     .extend_colors(total_clusters - length(manual_colors), base_colors)))
    } else {
      p_admixture <- p_admixture +
        ggplot2::scale_fill_manual(values = manual_colors[seq_len(total_clusters)])
    }
  } else {
    p_admixture <- switch(color_scale,
      "aaas" = p_admixture + ggsci::scale_fill_aaas(),
      "npg" = p_admixture + ggsci::scale_fill_npg(),
      "lancet" = p_admixture + ggsci::scale_fill_lancet(),
      "jco" = p_admixture + ggsci::scale_fill_jco(),
      "ucscgb" = p_admixture + ggsci::scale_fill_ucscgb(),
      "uchicago" = p_admixture + ggsci::scale_fill_uchicago(),
      "simpsons" = p_admixture + ggsci::scale_fill_simpsons(),
      "rickandmorty" = p_admixture + ggsci::scale_fill_rickandmorty()
    )
  }

  # ---- 10. Tree plot + combination ----
  p_tree <- NULL
  combined_plot <- NULL

  if (!is.null(phylo_tree)) {
    sample_group_data <- NULL
    if (!is.null(sample_group)) {
      if (is.character(sample_group) && length(sample_group) == 1L) {
        if (verbose) message("Reading sample group information from file: ", sample_group)
        sample_group_data <- tryCatch(
          readr::read_delim(sample_group, delim = "\t", show_col_types = FALSE),
          error = function(e) {
            warning("Failed to read sample group file: ", rlang::cnd_message(e))
            NULL
          })
      } else if (is.data.frame(sample_group)) {
        sample_group_data <- sample_group
      } else {
        warning("sample_group must be a file path or data frame; ignoring it")
      }
    }

    p_tree <- ggtree::ggtree(phylo_tree,
      branch.length = if (show_branch_length) NULL else "none",
      ladderize = FALSE,
      linewidth = tree_linewidth)

    if (!is.null(sample_group_data)) {
      required_cols <- c("label", "group")
      if (!all(required_cols %in% colnames(sample_group_data))) {
        warning("sample_group must contain 'label' and 'group' columns. Skipping group coloring.")
      } else {
        if (!"color" %in% colnames(sample_group_data)) {
          sample_group_data$color <- sample_group_data$group
        }
        p_tree <- ggtree::`%<+%`(p_tree, sample_group_data) +
          ggtree::geom_tippoint(ggplot2::aes(color = .data$color),
                                size = tip_point_size)
        p_tree <- switch(tree_color_scale,
          "default" = p_tree + ggplot2::scale_color_discrete(),
          "d3" = p_tree + ggsci::scale_color_d3(),
          "aaas" = p_tree + ggsci::scale_color_aaas(),
          "npg" = p_tree + ggsci::scale_color_npg(),
          "lancet" = p_tree + ggsci::scale_color_lancet(),
          "jco" = p_tree + ggsci::scale_color_jco(),
          "ucscgb" = p_tree + ggsci::scale_color_ucscgb(),
          "uchicago" = p_tree + ggsci::scale_color_uchicago(),
          "simpsons" = p_tree + ggsci::scale_color_simpsons(),
          "rickandmorty" = p_tree + ggsci::scale_color_rickandmorty()
        )
      }
    }

    combined_plot <- p_admixture %>%
      aplot::insert_left(p_tree, width = tree_width)
  }

  # ---- 11. Optional saving ----
  if (!is.null(save)) {
    main_plot <- if (!is.null(combined_plot)) combined_plot else p_admixture
    out_dir <- dirname(save)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    ggsave2(filename = basename(save), plot = main_plot, path = out_dir)
    if (verbose) message("Saved main plot to: ", save)
  }

  if (verbose) message("Analysis completed! Returning plot objects.")

  structure(
    list(
      combined_plot = combined_plot,
      admixture_plot = p_admixture,
      tree_plot = p_tree,
      samples = df_sample$sample,
      k_values = available_k,
      matched_samples = matched_samples
    ),
    class = "admixture_phylo_results"
  )
}

#' Print ADMIXTURE Phylogenetic Analysis Results
#'
#' @param x An object of class \code{admixture_phylo_results} from
#'   \code{\link{admixture_phylo_analysis}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- admixture_phylo_analysis("path/to/admixture/results/")
#' res   # auto-prints via print.admixture_phylo_results
#' }
print.admixture_phylo_results <- function(x, ...) {
  cat("ADMIXTURE phylogenetic analysis\n")
  cat("  Samples:    ", length(x$samples), "\n")
  cat("  K values:   ", paste(x$k_values, collapse = ", "), "\n")
  if (!is.null(x$tree_plot)) {
    cat("  Tree:       yes (", length(x$matched_samples), " matched samples)\n")
    cat("  Plots:      combined_plot, admixture_plot, tree_plot\n")
  } else {
    cat("  Tree:       no\n")
    cat("  Plots:      admixture_plot\n")
  }
  invisible(x)
}

#' Extend the base palette to n colors
#'
#' Fills the tail of the palette with ggplot2-style hues when n exceeds the
#' base palette length.
#'
#' @param n Number of colors needed.
#' @param base_colors Character vector of base colors.
#' @return A character vector of length \code{n}.
#' @keywords internal
#' @noRd
.extend_colors <- function(n, base_colors) {
  if (n <= length(base_colors)) {
    return(base_colors[seq_len(n)])
  }
  c(base_colors, scales::hue_pal()(n - length(base_colors)))
}

#' Remove surrounding quotes from tree labels
#'
#' @param x Character vector of tip labels.
#' @return Labels with single quotes removed.
#' @keywords internal
#' @noRd
.clean_tip_labels <- function(x) {
  stringr::str_remove_all(x, "\\'")
}