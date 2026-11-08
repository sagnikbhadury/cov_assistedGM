# Network drawing: Sagnik Bhadury's Network_VisualizationSB() (circular ggraph layout, edge width
# by |weight|, edge colour by sign, node size by connectivity), brought into the package from
# "Network Visualization.R". The data wrangling that used reshape2/dplyr is written in base R
# with the same semantics; the ggraph/ggplot2 drawing code is unchanged apart from the font
# and title arguments. pxg_network_matrix() prepares PxG output for it.

utils::globalVariables(c("Correlation", "node_size", "name", "x", "y"))

#' Circular network plot with edge sign, edge weight and node connectivity
#'
#' Draws a weighted undirected network from a symmetric matrix: every pair whose absolute
#' entry is at least \code{thres} becomes an edge drawn as an arc whose width is the absolute
#' value and whose colour is the sign; node size is the sum of absolute edge weights at the
#' node ("connectivity degree"); nodes without an edge are not drawn. This is the function
#' used for the network figures of this package (and of the author's proteomics work); the
#' same code can display partial correlations, posterior edge probabilities or any other
#' symmetric weight matrix.
#'
#' The matrix is plotted as given. An inverse covariance matrix would show reversed signs, so
#' convert it to partial correlations first (\code{\link{pxg_network_matrix}} does this for
#' PxG output). The diagonal is ignored.
#'
#' @param graph_mtx Symmetric numeric matrix with dimnames; the entries are the edge weights.
#' @param dinverse Unused; kept for compatibility with the original signature.
#' @param thres Edges with \code{abs(weight) < thres} are not drawn.
#' @param regionname Subtitle text (the original use was a slide or patient identifier).
#' @param title Plot title.
#' @param font Font family for all text; the original used "Calibri", which is not available
#'   on every system, so the default is \code{"sans"}.
#' @param label_size Size of the node labels.
#' @param seed Seed set at the start (the original sets 1234; the circular layout itself is
#'   deterministic).
#' @return A list with \code{graph} (igraph object with edge attribute \code{Correlation} and
#'   vertex attribute \code{node_size}), \code{node_range}, \code{edge_abs_range},
#'   \code{node_size}, \code{cor_phi} (all pairs with their weight), \code{cor_list} (the pairs
#'   drawn) and \code{myplot} (the ggplot object).
#' @examples
#' res <- readRDS(system.file("extdata", "brca_pxg_results.rds", package = "covassistedGM"))
#' g1 <- res$cluster_graphs[["Cluster 1"]]
#' M <- pxg_network_matrix(g1$prob, g1$adj)
#' out <- Network_VisualizationSB(M, NULL, thres = 0.5, regionname = "Cluster 1",
#'                                title = "PxG breast cancer network")
#' out$myplot
#' @export
Network_VisualizationSB <- function(graph_mtx, dinverse = NULL, thres, regionname,
                                    title = " Slide Specific Networks", font = "sans",
                                    label_size = 2.6, seed = 1234) {
  set.seed(seed)
  stopifnot(is.matrix(graph_mtx), nrow(graph_mtx) == ncol(graph_mtx), !is.null(colnames(graph_mtx)))
  graph_mtx <- as.matrix(graph_mtx)
  diag(graph_mtx) <- 0
  feature_f <- sort(colnames(graph_mtx))

  # all unordered pairs, one row each (the original melts the matrix and keeps, per pair, the
  # entry with the larger absolute value; in a symmetric matrix the two entries are equal)
  idx <- which(upper.tri(graph_mtx), arr.ind = TRUE)
  f1 <- rownames(graph_mtx)[idx[, 1]]; f2 <- colnames(graph_mtx)[idx[, 2]]
  v12 <- graph_mtx[idx]; v21 <- graph_mtx[idx[, c(2, 1), drop = FALSE]]
  value <- ifelse(abs(v12) >= abs(v21), v12, v21)
  # order the two names as the original does (feature1 the earlier level, feature2 the later)
  lev <- match(f1, feature_f) < match(f2, feature_f)
  cor_phi <- data.frame(Correlation = pmin(pmax(value, -1), 1),
                        feature1 = ifelse(lev, f1, f2), feature2 = ifelse(lev, f2, f1),
                        stringsAsFactors = FALSE)
  cor_phi$feature1 <- factor(cor_phi$feature1, levels = feature_f)
  cor_phi$feature2 <- factor(cor_phi$feature2, levels = feature_f)

  cor_phi_list <- cor_phi[abs(cor_phi$Correlation) >= thres, , drop = FALSE]
  if (nrow(cor_phi_list) == 0) stop("no edge reaches the threshold ", thres)

  # node connectivity: sum of |weight| over the edges at each node; isolated nodes are dropped
  feats <- unique(c(as.character(cor_phi_list$feature1), as.character(cor_phi_list$feature2)))
  feats <- feats[order(match(feats, feature_f))]
  s1 <- tapply(abs(cor_phi_list$Correlation), as.character(cor_phi_list$feature1), sum)
  s2 <- tapply(abs(cor_phi_list$Correlation), as.character(cor_phi_list$feature2), sum)
  node_size <- ifelse(is.na(s1[feats]), 0, s1[feats]) + ifelse(is.na(s2[feats]), 0, s2[feats])
  sum_all <- data.frame(feature = feats, node_size = as.numeric(node_size), id = seq_along(feats),
                        stringsAsFactors = FALSE)

  edges <- data.frame(from = as.character(cor_phi_list$feature1), to = as.character(cor_phi_list$feature2),
                      Correlation = cor_phi_list$Correlation, stringsAsFactors = FALSE)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = sum_all)
  node_range <- c(min(sum_all$node_size), max(sum_all$node_size))
  edge_abs_range <- c(min(abs(edges$Correlation)), max(abs(edges$Correlation)))

  p3 <- ggraph::ggraph(g, layout = "circle") +
    ggraph::geom_edge_arc(ggplot2::aes(color = factor(sign(Correlation), levels = c(1, -1)),
                                       width = abs(Correlation)), strength = 0.2, alpha = 0.75) +
    ggraph::scale_edge_width(name = "Magnitude", breaks = c(0.2, 0.4, 0.6, 0.8), range = 1 * edge_abs_range) +
    ggraph::scale_edge_color_discrete(name = "Sign", labels = c("Positive", "Negative"), drop = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(size = node_size), color = "blue") +
    ggplot2::scale_size(name = "Connectivity Degree", range = 2 * node_range) +
    ggplot2::theme_void() + ggplot2::coord_fixed() +
    ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.15, y = y * 1.15,
                                        angle = ifelse(atan(-(x / y)) * (180 / pi) < 0,
                                                       90 + atan(-(x / y)) * (180 / pi), 270 + atan(-x / y) * (180 / pi)),
                                        hjust = ifelse(x > 0, 0, 1)),
                           size = label_size, fontface = "bold") +
    ggplot2::expand_limits(x = c(-2, 2), y = c(-2, 2)) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1)) +
    ggplot2::theme(text = ggplot2::element_text(family = font, face = "bold"),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   plot.caption = ggplot2::element_text(hjust = 0.5),
                   legend.title = ggplot2::element_text(size = 6, face = "bold"),
                   legend.text = ggplot2::element_text(size = 5), legend.position = "bottom",
                   legend.box = "vertical", legend.margin = ggplot2::margin(0, 0, -0.5, -0.5),
                   legend.spacing.x = ggplot2::unit(0, "mm"), legend.spacing.y = ggplot2::unit(0, "mm")) +
    ggplot2::labs(title = title, subtitle = paste0(" ", regionname, " "))

  list(graph = g, node_range = node_range, edge_abs_range = edge_abs_range, node_size = sum_all$node_size,
       cor_phi = cor_phi, cor_list = cor_phi_list, myplot = p3)
}

#' Weight matrix for drawing a PxG cluster network
#'
#' Prepares the input of \code{\link{Network_VisualizationSB}} from PxG output.
#' \itemize{
#'   \item For the pseudo-likelihood models, pass the symmetrised posterior edge probabilities
#'     \code{prob} and the adjacency \code{adj} from \code{\link{cluster_graphs}} (or
#'     \code{\link{cluster_edge_probabilities}}); the weight of a selected edge is its
#'     posterior probability (positive, so all edges take the "Positive" colour) and
#'     non-selected pairs are set to 0.
#'   \item For the G-Wishart models, pass the posterior mean precision matrix of a sample or
#'     cluster (\code{Omega}, a slice of \code{fit$Osum}) and optionally its edge probability
#'     matrix (a slice of \code{fit$Gsum}); the weight is the signed partial correlation
#'     \eqn{-\Omega_{ij}/\sqrt{\Omega_{ii}\Omega_{jj}}}, masked by \code{Gsum >= cutoff} when
#'     given.
#' }
#'
#' @param prob Symmetric matrix of posterior edge probabilities (pseudo-likelihood route).
#' @param adj Optional 0/1 adjacency; default \code{prob >= cutoff}.
#' @param Omega Posterior mean precision matrix (G-Wishart route); overrides \code{prob}.
#' @param Gsum Optional edge probability matrix matching \code{Omega}.
#' @param cutoff Threshold applied to \code{prob} or \code{Gsum} when no \code{adj} is given.
#' @param node_names Optional dimnames to attach.
#' @return A symmetric weight matrix with unit diagonal and dimnames, ready for
#'   \code{Network_VisualizationSB(M, NULL, thres, ...)} with \code{thres} a small positive
#'   number (every non-zero weight drawn) or \code{cutoff} (the probability route).
#' @export
pxg_network_matrix <- function(prob = NULL, adj = NULL, Omega = NULL, Gsum = NULL, cutoff = 0.5,
                               node_names = NULL) {
  if (!is.null(Omega)) {
    d <- 1 / sqrt(diag(Omega))
    M <- -(diag(d) %*% Omega %*% diag(d))
    if (!is.null(Gsum)) M[Gsum < cutoff] <- 0
    dimnames(M) <- dimnames(Omega)
  } else {
    stopifnot(!is.null(prob))
    if (is.null(adj)) adj <- (prob >= cutoff) * 1
    M <- prob * (adj != 0)
    dimnames(M) <- dimnames(prob)
  }
  if (!is.null(node_names)) dimnames(M) <- list(node_names, node_names)
  if (is.null(colnames(M))) { nm <- paste0("V", seq_len(ncol(M))); dimnames(M) <- list(nm, nm) }
  diag(M) <- 1
  M
}

#' Draw every PxG cluster network with Network_VisualizationSB
#'
#' Convenience wrapper around \code{\link{pxg_network_matrix}} and
#' \code{\link{Network_VisualizationSB}} for the output of \code{\link{cluster_graphs}}.
#'
#' @param graphs The list returned by \code{\link{cluster_graphs}}.
#' @param title Plot title shared by the panels.
#' @param thres Threshold passed on; default the \code{cutoff} used for \code{adj}, so that
#'   exactly the selected edges are drawn.
#' @param ... Further arguments to \code{\link{Network_VisualizationSB}}.
#' @return A named list of ggplot objects, one per cluster.
#' @export
plot_pxg_networks <- function(graphs, title = "PxG cluster network", thres = 0.5, ...) {
  out <- lapply(names(graphs), function(k) {
    g <- graphs[[k]]
    M <- pxg_network_matrix(g$prob, g$adj)
    sub <- sprintf("%s: n = %d, %d edges", k, g$n, sum(g$adj) / 2)
    Network_VisualizationSB(M, NULL, thres = thres, regionname = sub, title = title, ...)$myplot
  })
  names(out) <- names(graphs)
  out
}
