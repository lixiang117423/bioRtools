#' Microbiome Network Analysis via SpiecEasi
#'
#' Constructs microbial association networks for each group using SpiecEasi,
#' computes network properties (centrality, hubs, clusters), and compares
#' networks across groups.
#'
#' @param data Count matrix or data frame. Rows = ASV/OTU, cols = samples.
#'   First column may be feature IDs (auto-detected if non-numeric).
#' @param sample Data frame of sample metadata. Must contain a column matching
#'   sample names in \code{data} and a \code{groupCol} column for grouping.
#' @param groupCol Column name in \code{sample} for grouping (default: "group").
#' @param taxonomy Optional data frame of feature taxonomy. Row names or first
#'   column should match feature IDs in \code{data}.
#' @param method Network construction method: "mb" (default), "glasso", or "cor".
#'   "mb" and "glasso" use sparse inverse covariance via huge+pulsar (memory-intensive).
#'   "cor" uses correlation with threshold filtering (fast, low memory).
#' @param cor_method Character string, correlation method for \code{method = "cor"}:
#'   "spearman" (default), "pearson", or "kendall".
#' @param minSamples Minimum sample prevalence to keep a feature. If < 1,
#'   treated as proportion; if >= 1, treated as count (default: 0.1 = 10\%).
#' @param minReads Minimum total reads to keep a feature (default: 10).
#' @param hubQuant Quantile threshold for hub detection (default: 0.95).
#' @param clustMethod igraph clustering method (default: "cluster_fast_greedy").
#' @param pulsar.params List passed to \code{SpiecEasi::spiec.easi} pulsar
#'   (default: \code{list(rep.num = 20)}).
#' @param verbose Print progress messages (default: TRUE).
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{networks}}{Named list of igraph objects, one per group.}
#'     \item{\code{adjaMats}}{Named list of adjacency matrices.}
#'     \item{\code{nodeProps}}{Data frame of node properties for all groups,
#'       with columns: feature, group, degree, betweenness, closeness,
#'       eigenvector, is_hub, cluster, and taxonomy columns if provided.}
#'     \item{\code{globalStats}}{Data frame with one row per group: n_nodes,
#'       n_edges, density, avg_path_length, clustering_coef, modularity,
#'       pos_edge_pct, n_hubs.}
#'     \item{\code{compare}}{Data frame of pairwise comparisons (when > 1
#'       group): Jaccard index of hub nodes and differences in global metrics.}
#'     \item{\code{spiecResults}}{Named list of raw spiec.easi output objects.}
#'   }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#' @seealso \code{\link[SpiecEasi]{spiec.easi}}
#'
#' @examples
#' \dontrun{
#' res <- microbiome_net(df.asv, df.sample, groupCol = "group", taxonomy = df.tax)
#' res$globalStats
#' head(res$nodeProps)
#' }
microbiome_net <- function(data, sample, groupCol = "group", taxonomy = NULL,
                           method = "mb", cor_method = "spearman",
                           minSamples = 0.1, minReads = 10,
                           cor.threshold = 0.5, cor.pvalue = 0.05,
                           hubQuant = 0.95, clustMethod = "cluster_fast_greedy",
                           pulsar.params = list(rep.num = 20),
                           verbose = TRUE) {

  # --- Prepare data ---------------------------------------------------------
  data <- as.data.frame(data)
  sample <- as.data.frame(sample)

  # Auto-detect feature ID column (first non-numeric column)
  if (!is.numeric(data[[1]])) {
    feat_ids <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    feat_ids <- rownames(data)
    if (is.null(feat_ids)) feat_ids <- paste0("ASV_", seq_len(nrow(data)))
  }
  data_mat <- as.matrix(data)
  rownames(data_mat) <- feat_ids

  # Auto-detect sample ID column in sample
  rn <- rownames(sample)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(sample))))) {
    for (col in names(sample)) {
      if (is.character(sample[[col]]) && all(colnames(data_mat) %in% sample[[col]])) {
        rn <- sample[[col]]
        rownames(sample) <- rn
        break
      }
    }
  }

  if (!groupCol %in% names(sample)) {
    stop(sprintf("Column '%s' not found in sample data", groupCol))
  }

  if (!cor_method %in% c("pearson", "spearman", "kendall")) {
    stop("'cor_method' must be 'pearson', 'spearman', or 'kendall'")
  }

  # Align samples
  common <- intersect(colnames(data_mat), rownames(sample))
  if (length(common) == 0) stop("No matching samples between data and sample")
  data_mat <- data_mat[, common, drop = FALSE]
  sample <- sample[common, , drop = FALSE]

  # Filter low-frequency features
  n_samp <- ncol(data_mat)
  if (minSamples < 1) {
    minSampCount <- ceiling(minSamples * n_samp)
  } else {
    minSampCount <- minSamples
  }
  prev <- rowSums(data_mat > 0)
  total <- rowSums(data_mat)
  keep <- prev >= minSampCount & total >= minReads
  data_mat <- data_mat[keep, , drop = FALSE]
  feat_ids <- rownames(data_mat)

  if (verbose) {
    message(sprintf("Features retained: %d / %d (prevalence >= %d, total reads >= %d)",
                    nrow(data_mat), length(keep), minSampCount, minReads))
  }

  if (is.null(groupCol)) {
    sample$.all <- "All"
    groupCol <- ".all"
  }
  groups <- unique(as.character(sample[[groupCol]]))

  # Prepare taxonomy
  if (!is.null(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy)
    if (!is.numeric(taxonomy[[1]])) {
      tax_ids <- taxonomy[[1]]
      taxonomy <- taxonomy[, -1, drop = FALSE]
    } else {
      tax_ids <- rownames(taxonomy)
    }
    rownames(taxonomy) <- tax_ids
  }

  # --- Build networks per group ---------------------------------------------
  networks <- list()
  adjaMats <- list()
  spiecResults <- list()
  nodeProps_list <- list()
  globalStats_list <- list()

  for (grp in groups) {
    if (verbose) message(sprintf("\n--- Group: %s ---", grp))

    idx <- sample[[groupCol]] == grp
    d_sub <- data_mat[, idx, drop = FALSE]

    # Remove features with zero variance in this group
    var_check <- apply(d_sub, 1, var)
    d_sub <- d_sub[var_check > 0, , drop = FALSE]

    if (verbose) message(sprintf("  Samples: %d, Features: %d", ncol(d_sub), nrow(d_sub)))

    # Warn if features >> samples (memory risk in huge/pulsar)
    if (nrow(d_sub) > ncol(d_sub) * 100) {
      warning(sprintf(
        "Group '%s': %d features / %d samples (ratio %.0f:1). Sparse estimation may be unstable or crash.\n  Consider increasing minSamples/minReads to reduce features to < %d.",
        grp, nrow(d_sub), ncol(d_sub), round(nrow(d_sub) / ncol(d_sub)),
        ncol(d_sub) * 50
      ))
    }

    # --- Build adjacency matrix ---
    if (method == "cor") {
      # Correlation + threshold (fast, low memory)
      if (verbose) message("  Method: ", cor_method, " correlation")

      # samples in rows, features in columns for corAndPvalue
      cp_result <- WGCNA::corAndPvalue(t(d_sub), method = cor_method)
      cor_mat <- cp_result$cor
      p_mat <- cp_result$p

      # Adjust p-values (upper triangle only, then mirror)
      n_feat <- nrow(cor_mat)
      p_upper <- p_mat[upper.tri(p_mat)]
      p_adj_upper <- stats::p.adjust(p_upper, method = "BH")

      p_adj <- matrix(1, n_feat, n_feat)
      p_adj[upper.tri(p_adj)] <- p_adj_upper
      p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]

      # Threshold: |r| >= threshold AND padj < pvalue
      adja <- cor_mat
      adja[abs(cor_mat) < cor.threshold | p_adj >= cor.pvalue] <- 0
      diag(adja) <- 0
      colnames(adja) <- rownames(adja) <- rownames(d_sub)

      # Store dummy objects for consistent output
      est <- NULL
      pulsar_out <- NULL
      opt_idx <- NA
    } else {
      # CLR transformation (samples x features)
      X <- t(d_sub) + 1  # pseudocount
      X_clr <- scale(log(X), center = TRUE, scale = FALSE)
      X_clr <- X_clr / sqrt(apply(X_clr^2, 1, sum))

      # Sparse inverse covariance via huge + pulsar
      if (verbose) message(sprintf("  Method: %s (huge + pulsar)", method))
      est <- tryCatch({
        huge::huge(X_clr, method = method, verbose = FALSE, cov.output = TRUE)
      }, error = function(e) {
        stop(sprintf("huge estimation failed for group '%s': %s", grp, e$message))
      })

      # Model selection via pulsar StARS
      pp <- pulsar.params
      if (is.null(pp$rep.num)) pp$rep.num <- 20
      if (is.null(pp$thresh)) pp$thresh <- 0.05
      pp$criterion <- "stars"

      pulsar_out <- tryCatch({
        pulsar::pulsar(X_clr,
          fun = function(data, ...) huge::huge(data, method = method, cov.output = TRUE, ...),
          fargs = list(lambda = est$lambda),
          criterion = pp$criterion,
          thresh = pp$thresh,
          rep.num = pp$rep.num)
      }, error = function(e) {
        stop(sprintf("Pulsar model selection failed for group '%s': %s", grp, e$message))
      })

      opt_idx <- pulsar::opt.index(pulsar_out, "stars")
      if (is.null(opt_idx) || is.na(opt_idx)) opt_idx <- length(est$lambda)

      # Extract adjacency matrix at optimal lambda
      if (method == "glasso") {
        icov <- est$icov[[opt_idx]]
        secor <- stats::cov2cor(solve(icov))
        refit <- est$path[[opt_idx]]
        adja <- as.matrix(secor * refit)
      } else {
        beta <- as.matrix(est$beta[[opt_idx]])
        adja <- (abs(beta) + t(abs(beta))) / 2
        sign_mat <- sign(beta + t(beta))
        adja <- adja * sign_mat
      }
      colnames(adja) <- rownames(adja) <- rownames(d_sub)
    }

    # Convert to igraph (absolute weights for structure, sign stored as edge attr)
    abs_adja <- abs(adja)
    g <- igraph::graph_from_adjacency_matrix(abs_adja, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
    # Assign edge signs from original adjacency
    edges <- igraph::as_edgelist(g, names = FALSE)
    edge_signs <- sapply(seq_len(nrow(edges)), function(i) sign(adja[edges[i, 1], edges[i, 2]]))
    igraph::E(g)$sign <- edge_signs

    # --- Analyze network properties -----------------------------------------
    # Dissimilarity for path-based measures
    diss <- 1 / abs_adja
    diag(diss) <- 0
    diss[diss == Inf] <- 0
    g_diss <- igraph::graph_from_adjacency_matrix(diss, mode = "undirected",
                                                   weighted = TRUE, diag = FALSE)

    # Centrality
    deg <- igraph::degree(g, normalized = TRUE)
    betw <- igraph::betweenness(g_diss, normalized = TRUE)

    # For closeness/eigenvector, use largest connected component if disconnected
    comps <- igraph::components(g)
    if (comps$no > 1) {
      lcc_idx <- which.max(comps$csize)
      lcc_nodes <- names(which(comps$membership == lcc_idx))
      g_lcc <- igraph::induced_subgraph(g, lcc_nodes)
      close_lcc <- igraph::closeness(g_lcc, normalized = TRUE)
      eigen_lcc <- igraph::eigen_centrality(g_lcc)$vector
      close_all <- rep(NA_real_, igraph::vcount(g))
      eigen_all <- rep(NA_real_, igraph::vcount(g))
      names(close_all) <- names(close_lcc) # placeholder
      names(close_all) <- igraph::V(g)$name
      names(eigen_all) <- igraph::V(g)$name
      close_all[lcc_nodes] <- close_lcc[lcc_nodes]
      eigen_all[lcc_nodes] <- eigen_lcc[lcc_nodes]
      closeness_vec <- close_all
      eigen_vec <- eigen_all
    } else {
      closeness_vec <- igraph::closeness(g, normalized = TRUE)
      eigen_vec <- igraph::eigen_centrality(g)$vector
    }

    # Hub detection
    eigen_thresh <- quantile(eigen_vec, probs = hubQuant, na.rm = TRUE)
    is_hub <- !is.na(eigen_vec) & eigen_vec >= eigen_thresh

    # Clustering
    clust_fn <- get(clustMethod, envir = asNamespace("igraph"))
    clust_res <- tryCatch(clust_fn(g), error = function(e) NULL)
    membership <- if (!is.null(clust_res)) igraph::membership(clust_res) else rep(NA_integer_, igraph::vcount(g))

    # Node properties
    node_df <- data.frame(
      feature = igraph::V(g)$name,
      group = grp,
      degree = deg,
      betweenness = betw,
      closeness = closeness_vec,
      eigenvector = eigen_vec,
      is_hub = is_hub,
      cluster = as.integer(membership),
      stringsAsFactors = FALSE
    )

    # Add taxonomy
    if (!is.null(taxonomy)) {
      tax_match <- taxonomy[node_df$feature, , drop = FALSE]
      node_df <- cbind(node_df, tax_match)
    }

    # Global stats
    n_edges <- igraph::ecount(g)
    n_nodes <- igraph::vcount(g)
    pos_pct <- if (n_edges > 0) sum(igraph::E(g)$sign > 0) / n_edges * 100 else NA_real_
    mod <- if (!is.null(clust_res)) igraph::modularity(clust_res) else NA_real_

    global_df <- data.frame(
      group = grp,
      n_nodes = n_nodes,
      n_edges = n_edges,
      density = igraph::edge_density(g),
      avg_path_length = if (n_nodes > 1) igraph::mean_distance(g) else NA_real_,
      clustering_coef = igraph::transitivity(g, type = "global"),
      modularity = mod,
      pos_edge_pct = round(pos_pct, 1),
      n_hubs = sum(is_hub),
      stringsAsFactors = FALSE
    )

    networks[[grp]] <- g
    adjaMats[[grp]] <- adja
    spiecResults[[grp]] <- list(est = est, pulsar = pulsar_out, opt_idx = opt_idx)
    nodeProps_list[[grp]] <- node_df
    globalStats_list[[grp]] <- global_df
  }

  if (length(nodeProps_list) == 0) {
    warning("No groups were successfully analyzed. Returning NULL.")
    return(NULL)
  }

  nodeProps <- do.call(rbind, nodeProps_list)
  rownames(nodeProps) <- NULL
  globalStats <- do.call(rbind, globalStats_list)
  rownames(globalStats) <- NULL

  # --- Pairwise comparison --------------------------------------------------
  compare <- NULL
  analyzed_groups <- names(nodeProps_list)
  if (length(analyzed_groups) > 1) {
    pairs <- utils::combn(analyzed_groups, 2, simplify = FALSE)
    compare_rows <- lapply(pairs, function(pair) {
      g1 <- pair[1]; g2 <- pair[2]
      hubs1 <- nodeProps_list[[g1]]$feature[nodeProps_list[[g1]]$is_hub]
      hubs2 <- nodeProps_list[[g2]]$feature[nodeProps_list[[g2]]$is_hub]

      # Jaccard index of hubs
      jacc_hub <- length(intersect(hubs1, hubs2)) / length(union(hubs1, hubs2))
      if (length(union(hubs1, hubs2)) == 0) jacc_hub <- NA_real_

      # Jaccard index of all nodes
      nodes1 <- nodeProps_list[[g1]]$feature
      nodes2 <- nodeProps_list[[g2]]$feature
      jacc_all <- length(intersect(nodes1, nodes2)) / length(union(nodes1, nodes2))

      s1 <- globalStats_list[[g1]]
      s2 <- globalStats_list[[g2]]

      data.frame(
        comparison = paste(g1, "vs", g2),
        jaccard_hubs = round(jacc_hub, 3),
        jaccard_nodes = round(jacc_all, 3),
        diff_n_edges = s1$n_edges - s2$n_edges,
        diff_density = round(s1$density - s2$density, 4),
        diff_modularity = round(s1$modularity - s2$modularity, 4),
        diff_avg_path = round(s1$avg_path_length - s2$avg_path_length, 4),
        stringsAsFactors = FALSE
      )
    })
    compare <- do.call(rbind, compare_rows)
    rownames(compare) <- NULL
  }

  if (verbose) {
    message("\n=== Microbiome Network Analysis Complete ===")
    message(sprintf("Groups analyzed: %d (%s)", length(groups), paste(groups, collapse = ", ")))
    message("\nGlobal statistics:")
    print(globalStats[, c("group", "n_nodes", "n_edges", "density", "modularity", "n_hubs")])
    if (!is.null(compare)) {
      message("\nPairwise comparison:")
      print(compare[, c("comparison", "jaccard_hubs", "diff_n_edges", "diff_density")])
    }
  }

  list(
    networks = networks,
    adjaMats = adjaMats,
    nodeProps = nodeProps,
    globalStats = globalStats,
    compare = compare,
    spiecResults = spiecResults
  )
}
