# Helper functions added for the package. Everything in this file is new (not from the
# supplement); each helper reproduces a post-processing step that the authors' scripts perform
# inline, and says so in its documentation.

#' Simulated data of Example 1 (piecewise linear precision matrix, 3 nodes)
#'
#' Generates the data of Section 4.1 of Niu et al. (2024) exactly as \code{Example1.R} does:
#' a scalar covariate \code{x} on three intervals, and a 3-node Gaussian graphical model
#' whose precision matrix is a different linear function of \code{x} on each interval, with a
#' different graph on each interval. \code{Y} and \code{X} are standardised as in the script.
#'
#' @param seed Random seed; the paper uses \code{2020}.
#' @return A list with \code{Y} (300 x 3), \code{X} (300 x 1), the raw covariate \code{x},
#'   and the functions \code{Omega(i, x)} and \code{G0(i, x)} giving the true precision
#'   matrix and graph on interval \code{i}.
#' @export
simulate_example1 <- function(seed = 2020) {
  set.seed(seed)
  x <- rep(0, 300)
  x[1:100]   <- runif(100, -1, -1/3)
  x[101:200] <- runif(100, -1/3, 1/3)
  x[201:300] <- runif(100, 1/3, 1)
  x <- x[order(x)]
  Omega <- function(i, x) {
    Omega <- diag(3) * 1.2
    if (i == 1) { Omega[1,2] <- Omega[2,1] <- 0;              Omega[1,3] <- Omega[3,1] <- -0.75*x+0.25; Omega[2,3] <- Omega[3,2] <-  0.75*x+1.25 }
    if (i == 2) { Omega[1,2] <- Omega[2,1] <-  0.75*x+0.75;   Omega[1,3] <- Omega[3,1] <- 0;            Omega[2,3] <- Omega[3,2] <- -0.75*x+0.75 }
    if (i == 3) { Omega[1,2] <- Omega[2,1] <- -0.75*x+1.25;   Omega[1,3] <- Omega[3,1] <-  0.75*x+0.25; Omega[2,3] <- Omega[3,2] <- 0 }
    Omega
  }
  G0 <- function(i, x) {
    G <- diag(3)
    if (i == 1) { G[1,3] <- G[3,1] <- 1; G[2,3] <- G[3,2] <- 1 }
    if (i == 2) { G[1,2] <- G[2,1] <- 1; G[2,3] <- G[3,2] <- 1 }
    if (i == 3) { G[1,2] <- G[2,1] <- 1; G[1,3] <- G[3,1] <- 1 }
    G
  }
  y <- matrix(0, 300, 3)
  for (i in 1:100)   { temp <- solve(Omega(1, x[i])); y[i,] <- t(t(chol(temp)) %*% rnorm(3, 0, 1)) }
  for (i in 101:200) { temp <- solve(Omega(2, x[i])); y[i,] <- t(t(chol(temp)) %*% rnorm(3, 0, 1)) }
  for (i in 201:300) { temp <- solve(Omega(3, x[i])); y[i,] <- t(t(chol(temp)) %*% rnorm(3, 0, 1)) }
  xm <- matrix(x, 300, 1)
  Y <- apply(y,  2, function(v) (v - mean(v)) / stats::sd(v))
  X <- apply(xm, 2, function(v) (v - mean(v)) / stats::sd(v))
  list(Y = Y, X = X, x = x, Omega = Omega, G0 = G0)
}

#' Simulated data of Example 2 (linear precision matrix, constant chain graph, 5 nodes)
#'
#' Generates the data of Section 4.2 of Niu et al. (2024) exactly as \code{Example2.R} does:
#' \code{x ~ U(-0.8, 0.8)}, a 5-node chain graph 1-2-3-4-5 with all off-diagonal precision
#' entries equal to \code{x} and diagonal 1.4.
#'
#' @param seed Random seed; the paper uses \code{2020}.
#' @return A list with \code{Y} (300 x 5), \code{X} (300 x 1), \code{x}, the function
#'   \code{Omega(x)} and the true adjacency matrix \code{G0}.
#' @export
simulate_example2 <- function(seed = 2020) {
  set.seed(seed)
  x <- runif(300, -0.8, 0.8)
  G0 <- diag(5) - diag(5)
  G0[1,2] <- G0[2,3] <- G0[3,4] <- G0[4,5] <- 1
  G0 <- G0 + t(G0)
  Omega <- function(x) { Om <- x * G0; diag(Om) <- 1.4; Om }
  y <- matrix(0, 300, 5)
  for (i in 1:300) { temp <- solve(Omega(x[i])); y[i,] <- t(t(chol(temp)) %*% rnorm(5, 0, 1)) }
  diag(G0) <- 1
  list(Y = y, X = matrix(x, 300, 1), x = x, Omega = Omega, G0 = G0)
}

#' Simulated data of Example 3 (two clusters, 50 nodes, 10 covariates)
#'
#' Generates the data of Section 4.3 of Niu et al. (2024) exactly as \code{Example3.R} does:
#' 500 samples in two equal clusters whose 10 covariates are \code{N(0, 1)} and \code{N(2, 1)}
#' respectively, and whose 50-node precision matrices are drawn from a G-Wishart distribution
#' over two independent sparse random graphs (edge probability 0.0085).
#'
#' @param seed Random seed; the paper uses \code{2022}.
#' @return A list with \code{Y} (500 x 50), \code{X} (500 x 10), the true cluster labels
#'   \code{z} and the two true adjacency matrices \code{G1}, \code{G2}.
#' @export
simulate_example3 <- function(seed = 2022) {
  set.seed(seed)
  p <- 10; q <- 50; n <- 500
  X <- matrix(0, n, p)
  for (i in 1:(n/2))     X[i,] <- rnorm(p, 0, 1)
  for (i in (n/2+1):n)   X[i,] <- rnorm(p, 2, 1)
  Y <- matrix(0, n, q)
  Omega1 <- matrix(0, q, q); Omega1[upper.tri(Omega1)] <- rbinom(q*(q-1)/2, 1, 0.0085)
  Omega1 <- Omega1 + t(Omega1); diag(Omega1) <- 1; G1 <- Omega1
  Omega2 <- matrix(0, q, q); Omega2[upper.tri(Omega2)] <- rbinom(q*(q-1)/2, 1, 0.0085)
  Omega2 <- Omega2 + t(Omega2); diag(Omega2) <- 1; G2 <- Omega2
  Omega1 <- BDgraph::rgwish(1, Omega1, 3, diag(q))
  Omega2 <- BDgraph::rgwish(1, Omega2, 3, diag(q))
  Sigma1 <- solve(Omega1); Sigma2 <- solve(Omega2)
  for (i in 1:(n/2))   Y[i,] <- t(t(chol(Sigma1)) %*% rnorm(q, 0, 1))
  for (i in (n/2+1):n) Y[i,] <- t(t(chol(Sigma2)) %*% rnorm(q, 0, 1))
  list(Y = Y, X = X, z = rep(1:2, each = n/2), G1 = G1, G2 = G2)
}

#' Symmetrise a neighbourhood-selection matrix
#'
#' \code{\link{Graph_NS}} and the \code{Gsum} output of \code{\link{PxG_PL}} describe a graph
#' by directed neighbourhood regressions, so the \code{(i, j)} and \code{(j, i)} frequencies
#' differ. The authors' scripts combine them with the maximum ("OR") rule and set the
#' diagonal to 1; this function does exactly that.
#'
#' @param G A square numeric matrix of neighbourhood-selection frequencies.
#' @param rule \code{"max"} (the authors' choice), \code{"min"} ("AND" rule) or \code{"mean"}.
#' @param diag_value Value placed on the diagonal (the scripts use 1).
#' @return A symmetric matrix with the same dimnames as \code{G}.
#' @export
symmetrize_graph <- function(G, rule = c("max", "min", "mean"), diag_value = 1) {
  rule <- match.arg(rule)
  stopifnot(is.matrix(G), nrow(G) == ncol(G))
  S <- switch(rule, max = pmax(G, t(G)), min = pmin(G, t(G)), mean = (G + t(G)) / 2)
  diag(S) <- diag_value
  dimnames(S) <- dimnames(G)
  S
}

#' Turn the \code{Gsum} array of \code{PxG_PL} into edge probability matrices
#'
#' \code{\link{PxG_PL}} returns \code{Gsum}, an integer array of dimension
#' \code{(q-1) x q x n} holding, for each sample, the number of kept iterations in which each
#' neighbourhood-regression indicator of that sample's cluster graph was on. This helper
#' divides by \code{nmcmc}, re-inserts the diagonal to obtain \code{q x q} matrices, averages
#' over the samples of each cluster, and symmetrises with \code{\link{symmetrize_graph}}.
#' The authors' scripts do not use \code{Gsum} for their figures (they refit each cluster's
#' graph with \code{\link{Graph_NS}}, see \code{\link{cluster_graphs}}); this function offers
#' the MCMC-averaged alternative.
#'
#' @param fit A list returned by \code{\link{PxG_PL}}.
#' @param nmcmc The number of kept iterations used in the fit.
#' @param labels Cluster labels to average over; default \code{fit$Label}.
#' @param node_names Optional node names for the dimnames.
#' @param rule Symmetrisation rule passed to \code{\link{symmetrize_graph}}.
#' @return A list of symmetric \code{q x q} edge probability matrices, one per cluster, in
#'   increasing order of the label values.
#' @export
cluster_edge_probabilities <- function(fit, nmcmc, labels = fit$Label, node_names = NULL,
                                       rule = "max") {
  Gx <- fit$Gsum
  stopifnot(length(dim(Gx)) == 3, dim(Gx)[1] == dim(Gx)[2] - 1)
  q <- dim(Gx)[2]; n <- dim(Gx)[3]
  stopifnot(length(labels) == n)
  to_square <- function(M) {           # M is (q-1) x q: column i = regression of node i on the others
    S <- matrix(0, q, q)
    for (i in seq_len(q)) S[-i, i] <- M[, i]
    S
  }
  out <- lapply(sort(unique(labels)), function(k) {
    idx <- which(labels == k)
    acc <- matrix(0, q, q)
    for (s in idx) acc <- acc + to_square(Gx[, , s])
    P <- symmetrize_graph(acc / (length(idx) * nmcmc), rule = rule, diag_value = 1)
    if (!is.null(node_names)) dimnames(P) <- list(node_names, node_names)
    P
  })
  names(out) <- paste0("Cluster ", sort(unique(labels)))
  out
}

#' Estimate one graph per cluster with Bayesian neighbourhood selection
#'
#' Reproduces the post-processing step of the authors' \code{Example3.R} and
#' \code{RealDataAnalysis.R}: for each cluster of a PxG fit, run \code{\link{Graph_NS}} on the
#' rows of \code{Y} allocated to that cluster, symmetrise the result with the maximum rule and
#' threshold it at \code{cutoff} to obtain an adjacency matrix. The seed is reset before each
#' cluster, as the scripts do (\code{set.seed(2022)} before every call).
#'
#' @param Y Graphical data, \code{n x q}.
#' @param labels Integer cluster labels of length \code{n}, typically \code{fit$Label}.
#' @param spike,slab,a1,a2,alphaG,burnin,nmcmc Passed to \code{\link{Graph_NS}}. The defaults
#'   are the values used for the breast cancer analysis (paper Figure 12); Example 3 uses
#'   \code{spike = 0.01, alphaG = 0.05, nmcmc = 500}.
#' @param cutoff Posterior probability above which an edge is declared (the paper uses 0.5).
#' @param seed Seed set before each cluster's sampler; \code{NULL} to leave the RNG alone.
#' @param verbose Passed to \code{\link{Graph_NS}}.
#' @return A list with one element per cluster (named \code{"Cluster k"}), each a list with
#'   \code{prob} (symmetrised posterior edge probabilities), \code{adj} (0/1 adjacency with
#'   zero diagonal), \code{n} (cluster size) and \code{raw} (the asymmetric \code{Graph_NS} output).
#' @export
cluster_graphs <- function(Y, labels, spike = 0.1, slab = 5, a1 = 1, a2 = 1, alphaG = 0.5,
                           burnin = 100, nmcmc = 1000, cutoff = 0.5, seed = 2022, verbose = FALSE) {
  stopifnot(nrow(Y) == length(labels))
  ks <- sort(unique(labels))
  out <- lapply(ks, function(k) {
    if (!is.null(seed)) set.seed(seed)
    idx <- which(labels == k)
    raw <- Graph_NS(YData = Y[idx, , drop = FALSE], spike = spike, slab = slab, a1 = a1, a2 = a2,
                    alphaG = alphaG, burnin = burnin, nmcmc = nmcmc, verbose = verbose)
    dimnames(raw) <- list(colnames(Y), colnames(Y))
    prob <- symmetrize_graph(raw, rule = "max", diag_value = 1)
    adj <- (prob >= cutoff) * 1; diag(adj) <- 0
    list(prob = prob, adj = adj, n = length(idx), raw = raw)
  })
  names(out) <- paste0("Cluster ", ks)
  out
}

#' Draw a cluster graph with igraph
#'
#' Plots the adjacency matrix of one cluster graph (as returned by \code{\link{cluster_graphs}})
#' with edge widths proportional to the posterior edge probabilities. A fixed circular layout
#' is used so that the same node sits in the same place in every cluster's plot.
#'
#' @param prob Symmetric matrix of posterior edge probabilities with node names as dimnames.
#' @param adj Optional 0/1 adjacency; default \code{prob >= cutoff}.
#' @param cutoff Threshold used when \code{adj} is missing.
#' @param main Plot title.
#' @param layout An igraph layout matrix; default circular in the column order of \code{prob}.
#' @param vertex_color Fill colour of the nodes.
#' @param edge_color Edge colour.
#' @param vertex_size,label_cex Node size and label scaling passed to \code{igraph::plot.igraph}.
#' @param ... Further arguments to \code{igraph::plot.igraph}.
#' @return The igraph object, invisibly.
#' @export
plot_cluster_network <- function(prob, adj = NULL, cutoff = 0.5, main = NULL, layout = NULL,
                                 vertex_color = "#DDEBF7", edge_color = "#1F4E79",
                                 vertex_size = 26, label_cex = 0.8, ...) {
  if (is.null(adj)) { adj <- (prob >= cutoff) * 1; diag(adj) <- 0 }
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  el <- igraph::as_edgelist(g, names = FALSE)
  w <- if (nrow(el) > 0) prob[cbind(el[, 1], el[, 2])] else numeric(0)
  if (is.null(layout)) layout <- igraph::layout_in_circle(g)
  igraph::plot.igraph(g, layout = layout, main = main,
                      vertex.color = vertex_color, vertex.frame.color = edge_color,
                      vertex.size = vertex_size, vertex.label.cex = label_cex,
                      vertex.label.color = "black", vertex.label.family = "sans",
                      edge.color = edge_color, edge.width = 1 + 4 * (w - cutoff) / (1 - cutoff), ...)
  invisible(g)
}

#' Heat map of posterior edge probabilities
#'
#' Base-graphics version of the authors' tile plots (paper Figures 8, 9 and 12): grey scale
#' from white (0) to black (1), with a dot marking each edge whose probability exceeds
#' \code{cutoff}. Rows are drawn top to bottom in the column order of \code{prob}.
#'
#' @param prob Symmetric matrix of posterior edge probabilities with dimnames.
#' @param cutoff Threshold for the dots; \code{NULL} for no dots.
#' @param main Plot title.
#' @param dot_col Dot colour (the authors use \code{"#EF6F6A"}).
#' @param cex_axis Axis label size.
#' @return \code{NULL}, invisibly.
#' @export
plot_edge_heatmap <- function(prob, cutoff = 0.5, main = NULL, dot_col = "#EF6F6A", cex_axis = 0.7) {
  q <- nrow(prob); nm <- colnames(prob); if (is.null(nm)) nm <- as.character(seq_len(q))
  M <- prob; diag(M) <- 1
  op <- par(mar = c(5, 5, 3, 1)); on.exit(par(op))
  image(x = seq_len(q), y = seq_len(q), z = t(M[q:1, ]), col = gray(seq(1, 0, length.out = 64)),
        zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = main)
  axis(1, at = seq_len(q), labels = nm, las = 2, cex.axis = cex_axis)
  axis(2, at = seq_len(q), labels = rev(nm), las = 1, cex.axis = cex_axis)
  if (!is.null(cutoff)) {
    E <- which(M >= cutoff & row(M) != col(M), arr.ind = TRUE)
    if (nrow(E) > 0) points(E[, 2], q + 1 - E[, 1], pch = 16, col = dot_col, cex = 1.1)
  }
  invisible(NULL)
}

#' Within-cluster sum of squares of the covariates
#'
#' Section 5 of Niu et al. (2024) compares the PxG and graph-only partitions by how
#' homogeneous the covariates are within clusters, with a "total sum of squares" of 509029.4
#' for PxG and 645584.9 for the graph-only model (576.5 and 731.1 per tumour). This function
#' returns three quantities:
#' \describe{
#'   \item{\code{ss}}{the total within-cluster sum of squares \eqn{\sum_k SS_k}, with
#'     \eqn{SS_k} the sum of squared deviations of the covariates from their column means in
#'     cluster \eqn{k};}
#'   \item{\code{weighted_ss}}{\eqn{\sum_k n_k \, SS_k}: the script accumulates each cluster's
#'     sum of squares inside a loop over the cluster's members, which multiplies it by the
#'     cluster size;}
#'   \item{\code{script_value}}{the number the script actually prints. Its expression
#'     \code{sum((RNA[idx, 1:3] - m)^2)} subtracts the length-3 mean vector \code{m} from an
#'     \code{n_k x 3} matrix, and R recycles \code{m} down the columns, so most entries have a
#'     different column's mean subtracted. \code{script_value} evaluates exactly that
#'     expression and reproduces 509029.4 and 645584.9 for the authors' partitions.}
#' }
#' The three agree on the direction of the comparison (PxG clusters are more homogeneous in the
#' covariates than graph-only clusters) but not on the scale; \code{ss} is the quantity the
#' paper describes.
#'
#' @param X Covariate matrix, \code{n x p}.
#' @param labels Cluster labels of length \code{n}.
#' @return A list with \code{ss}, \code{weighted_ss}, \code{script_value}, \code{per_cluster}
#'   (a data frame with cluster, size, SS) and \code{n}.
#' @export
within_cluster_ss <- function(X, labels) {
  stopifnot(nrow(X) == length(labels))
  ks <- sort(unique(labels))
  per <- do.call(rbind, lapply(ks, function(k) {
    idx <- which(labels == k)
    Xi <- X[idx, , drop = FALSE]
    m <- colMeans(Xi)
    ss <- sum(sweep(Xi, 2, m)^2)
    script_k <- 0
    for (i in idx) script_k <- script_k + sum((Xi - m)^2)   # the script's expression, recycling included
    data.frame(cluster = k, size = length(idx), ss = ss, script_value = script_k)
  }))
  list(ss = sum(per$ss), weighted_ss = sum(per$size * per$ss), script_value = sum(per$script_value),
       per_cluster = per, n = nrow(X))
}

#' Posterior partial correlations from a G-Wishart PxG fit
#'
#' \code{\link{PxG_GW}} and \code{\link{GraphOnly_GW}} return \code{Osum}, the posterior mean
#' precision matrix of every sample. The authors' \code{Example1.R} and \code{Example2.R}
#' convert these to partial correlations \eqn{-\Omega_{ij}/\sqrt{\Omega_{ii}\Omega_{jj}}}
#' without the minus sign (they plot \eqn{\Omega_{ij}/\sqrt{\Omega_{ii}\Omega_{jj}}} against
#' the true value computed the same way). This helper returns the same quantity for one pair
#' of nodes, for all samples, so the paper's figures can be reproduced.
#'
#' @param fit A list returned by \code{\link{PxG_GW}} or \code{\link{GraphOnly_GW}}.
#' @param i,j Node indices.
#' @param signed If \code{TRUE}, return the conventional partial correlation (with the minus
#'   sign). Default \code{FALSE} follows the scripts.
#' @return A numeric vector of length \code{n}.
#' @export
partial_correlation_path <- function(fit, i, j, signed = FALSE) {
  O <- fit$Osum
  r <- O[i, j, ] / sqrt(O[i, i, ] * O[j, j, ])
  if (signed) -r else r
}
