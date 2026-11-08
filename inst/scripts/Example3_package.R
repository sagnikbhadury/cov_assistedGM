### Section 4.3 "Example 3: piecewise constant precision matrix" of Niu et al. (2024), run with
### the covassistedGM package. Same models, arguments and seeds as the authors' Example3.R.

library(covassistedGM)
sim <- simulate_example3(seed = 2022)
Y <- sim$Y; X <- sim$X; p <- 10; q <- 50; n <- 500

## the three pseudo-likelihood models (a few minutes each)
set.seed(2022)
A <- PxG_PL(YData = Y, XData = X, spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1,
            mu0 = rep(0, p), sigsq0 = 1, b1 = 1, b2 = 0.1, M = 1, N = 3,
            burnin = 500, nmcmc = 1000, z0 = sample(1:3, n, TRUE, c(1/3, 1/3, 1/3)), B = 1)
A$DIC  # 104164.7
set.seed(2022)
B <- CovOnly_PL(YData = Y, XData = X, spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1,
                mu0 = rep(0, p), sigsq0 = 1, b1 = 1, b2 = 0.1, M = 1, N = 3,
                burnin = 500, nmcmc = 1000, z0 = sample(1:3, n, TRUE, c(1/3, 1/3, 1/3)))
B$DIC  # 155453.3
set.seed(2022)
C <- GraphOnly_PL(YData = Y, XData = X, spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1,
                  mu0 = rep(0, p), sigsq0 = 1, b1 = 1, b2 = 0.1, M = 1, N = 3,
                  burnin = 500, nmcmc = 1000, z0 = sample(1:3, n, TRUE, c(1/3, 1/3, 1/3)), B = 1)
C$DIC  # 106138.6

## Figure 7: clustering results in the first two covariates
cols <- c("#1170AA", "#EF6F6A", "#009E73")
op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
plot(X[, 1], X[, 2], col = cols[sim$z], pch = 16, xlab = "x1", ylab = "x2", main = "(A) truth")
plot(X[, 1], X[, 2], col = cols[B$Label], pch = 16, xlab = "x1", ylab = "x2", main = "(B) covariate-only")
plot(X[, 1], X[, 2], col = cols[C$Label], pch = 16, xlab = "x1", ylab = "x2", main = "(C) graph-only")
par(op)

## Figure 8: per-cluster graphs under PxG (Graph_NS with the script's settings)
graphs <- cluster_graphs(Y, A$Label, spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.05,
                         burnin = 100, nmcmc = 500, cutoff = 0.5, seed = 2022)
op <- par(mfrow = c(1, 2))
for (k in seq_along(graphs)) plot_edge_heatmap(graphs[[k]]$prob, cutoff = 0.5, main = names(graphs)[k], cex_axis = 0.4)
par(op)

## Figure 9A: BM13
set.seed(2022)
bm <- LSGraph(Y = Y, X0 = X, bG = 3.0, DG = diag(50), k = 0, which = 1, m = 1, M = 1, Kmat = diag(2),
              c = 500, alpha = 0.5, beta = 0.5, gamma0 = rep(0, 10), Adj0 = diag(50),
              burnin = 1000, nmc = 2000, verbose = FALSE)
plot_edge_heatmap(bm$Adj.mean, cutoff = 0.5, main = "(A) BM13", cex_axis = 0.4)

## Figure 9B: mSSL (GitHub package, optional)
if (requireNamespace("mSSL", quietly = TRUE)) {
  set.seed(2022)
  ms <- fit_mSSL(X, Y)
  plot_edge_heatmap(ms$abs_partial_cor, cutoff = NULL, main = "(B) mSSL", cex_axis = 0.4)
}

## Figure 9C: JRNS, initialised from lasso + graphical lasso as in the script
init <- jrns_initial_values(X, Y)          # needs glmnet and glasso
fit_JRNS <- JRNS(B_ini = init$B, Omega_ini = init$Omega, Y = Y, X = X, nmc = 1000, burnin = 2000)
freq <- jrns_edge_frequencies(fit_JRNS)
plot_edge_heatmap(freq, cutoff = 0.5, main = "(C) JRNS", cex_axis = 0.4)
