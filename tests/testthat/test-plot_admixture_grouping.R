make_group_tree <- function(extra_tip = NULL) {
  newick <- if (is.null(extra_tip)) {
    "((A:1,B:1):1,(C:1,D:1):1,E:2);"
  } else {
    paste0("((A:1,B:1):1,(C:1,D:1):1,E:2,", extra_tip, ":2);")
  }
  ape::read.tree(text = newick)
}

make_group_df <- function(extra_sample = NULL) {
  df <- data.frame(
    FID = c("A", "B", "C", "D", "E"),
    IID = c("A", "B", "C", "D", "E"),
    K8 = c(1, 2, 1, 2, 3),
    stringsAsFactors = FALSE
  )
  if (!is.null(extra_sample)) {
    df <- rbind(df, data.frame(FID = extra_sample, IID = extra_sample,
                               K8 = 4, stringsAsFactors = FALSE))
  }
  df
}

tree_order_rev <- function(tree) {
  rev(ggtree::get_taxa_name(ggtree::ggtree(tree)))
}

test_that("x-axis sample order follows the tree tip order (reversed)", {
  tr <- make_group_tree()
  p <- plot_admixture_grouping(tree = tr, data = make_group_df(), k = 8)

  expect_s3_class(p, "ggplot")
  x_limits <- p$scales$get_scales("x")$limits
  expect_equal(x_limits, tree_order_rev(tr))
})

test_that("accepts tree and data file paths", {
  dir <- tempfile("grp")
  dir.create(dir)
  writeLines("((A:1,B:1):1,(C:1,D:1):1,E:2);",
             file.path(dir, "tree.nwk"))
  utils::write.csv(make_group_df(), file.path(dir, "grouping.csv"),
                   row.names = FALSE)

  p <- plot_admixture_grouping(tree = file.path(dir, "tree.nwk"),
                               data = file.path(dir, "grouping.csv"), k = 8)
  expect_s3_class(p, "ggplot")
  expect_equal(p$scales$get_scales("x")$limits,
               tree_order_rev(ape::read.tree(file.path(dir, "tree.nwk"))))
})

test_that("character k is used as the column name directly", {
  tr <- make_group_tree()
  df <- make_group_df()

  p_num <- plot_admixture_grouping(tr, df, k = 8)
  p_char <- plot_admixture_grouping(tr, df, k = "K8")
  expect_equal(p_char$scales$get_scales("x")$limits,
               p_num$scales$get_scales("x")$limits)
  expect_equal(p_char$data$cluster, p_num$data$cluster)
})

test_that("numeric k past 9 selects K10 with numerically sorted levels", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$K10 <- c(1, 10, 2, 9, 8)

  p <- plot_admixture_grouping(tr, df, k = 10)
  expect_equal(levels(p$data$cluster), c("1", "2", "8", "9", "10"))
})

test_that("id_col selects the sample ID column", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$IID <- c("A", "B", "C", "D", "E")  # unchanged; FID is the default

  p <- plot_admixture_grouping(tr, df, k = 8, id_col = "IID")
  expect_equal(p$scales$get_scales("x")$limits, tree_order_rev(tr))
})

test_that("unmatched samples raise warnings, not silence", {
  tr <- make_group_tree(extra_tip = "X5")
  df <- make_group_df(extra_sample = "F")

  expect_warning(
    p <- plot_admixture_grouping(tr, df, k = 8),
    "in data but not in the tree"
  )
  expect_warning(
    p <- plot_admixture_grouping(tr, df, k = 8),
    "in the tree but not in data"
  )
  x_limits <- p$scales$get_scales("x")$limits
  expect_true("X5" %in% x_limits)  # tree-only sample keeps an empty slot
  expect_false("F" %in% x_limits)  # data-only sample is dropped
})

test_that("NA cluster values are dropped with a warning", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$K8[1] <- NA

  expect_warning(p <- plot_admixture_grouping(tr, df, k = 8),
                 "no cluster value")
  # A stays as an empty slot on the axis
  expect_true("A" %in% p$scales$get_scales("x")$limits)
})

test_that("invalid inputs raise errors", {
  tr <- make_group_tree()
  df <- make_group_df()

  expect_error(plot_admixture_grouping(tree = 42, data = df, k = 8),
               "'tree' must be a phylo object or a path")
  expect_error(plot_admixture_grouping(tree = tr, data = 42, k = 8),
               "'data' must be a data frame or a path")
  expect_error(plot_admixture_grouping(tr, df, k = 99),
               "Cluster column 'K99_qgroup'/'K99' not found")
  expect_error(plot_admixture_grouping(tr, df, k = "K9"),
               "Cluster column 'K9' not found")
  expect_error(plot_admixture_grouping(tr, df, id_col = "nope"),
               "'id_col' \\(nope\\) not found")
  expect_error(plot_admixture_grouping(tr, df, k = 1.5),
               "'k' must be a single positive integer")
})
test_that("numeric k selects the K<k>_qgroup column", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$K8_qgroup <- c(2, 1, 3, 2, 1)

  p <- plot_admixture_grouping(tr, df, k = 8)
  expect_equal(levels(p$data$cluster), c("1", "2", "3"))
  expect_equal(as.character(p$data$cluster), c("2", "1", "3", "2", "1"))
})

test_that("K<k>_qgroup takes precedence over legacy K<k>", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$K8_qgroup <- c(2, 2, 2, 2, 2)  # both K8 and K8_qgroup present

  p <- plot_admixture_grouping(tr, df, k = 8)
  expect_equal(as.character(p$data$cluster), rep("2", 5))
})

test_that("numeric k past 9 selects K10_qgroup with numerically sorted levels", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$K10_qgroup <- c(1, 10, 2, 9, 8)

  p <- plot_admixture_grouping(tr, df, k = 10)
  expect_equal(levels(p$data$cluster), c("1", "2", "8", "9", "10"))
})

test_that("character k accepts the qgroup column name directly", {
  tr <- make_group_tree()
  df <- make_group_df()
  df$K8_qgroup <- c(2, 1, 3, 2, 1)

  p <- plot_admixture_grouping(tr, df, k = "K8_qgroup")
  expect_equal(as.character(p$data$cluster), c("2", "1", "3", "2", "1"))

  # missing qgroup column is reported among the available K columns
  df2 <- df[, c("FID", "IID", "K8_qgroup")]
  expect_error(plot_admixture_grouping(tr, df2, k = 9),
               "Available K columns: K8_qgroup")
})

test_that("no matching samples raises an error", {
  tr <- make_group_tree(extra_tip = "X5")  # tips A-E + X5
  df <- make_group_df()
  df$FID <- c("P", "Q", "R", "S", "T")     # none shared with the tree

  expect_error(plot_admixture_grouping(tr, df, k = 8), "No samples match")
})
test_that("legend is hidden by default and opt-in via show_legend", {
  tr <- make_group_tree()
  df <- make_group_df()

  legend_drawn <- function(p) {
    g <- ggplot2::ggplotGrob(p)
    idx <- which(grepl("^guide-box", g$layout$name))
    any(vapply(idx, function(i) !inherits(g$grobs[[i]], "zeroGrob"),
               logical(1)))
  }

  p_default <- plot_admixture_grouping(tr, df, k = 8)
  expect_false(legend_drawn(p_default))

  p_legend <- plot_admixture_grouping(tr, df, k = 8, show_legend = TRUE)
  expect_true(legend_drawn(p_legend))

  expect_error(plot_admixture_grouping(tr, df, k = 8, show_legend = "yes"),
               "'show_legend' must be a single logical value")
})
