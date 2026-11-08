### Section 4.2 "Example 2: linear precision matrix with constant graph" of Niu et al. (2024),
### run with the covassistedGM package. Same models, arguments and seeds as the authors' Example2.R.

library(covassistedGM)
sim <- simulate_example2(seed = 2020)
Y <- sim$Y; X <- sim$X; x <- sim$x

## PxG (about 10 minutes)
set.seed(2020)
A <- PxG_GW(Y = Y, X = X, b = 3, D = diag(5), alphaG = 0.2,
            mu0 = 0, sigsq0 = 1, b1 = 10, b2 = 0.1, M = 1, N = 15,
            burnin = 1000, nmcmc = 2000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 3)

## Figure 5: partial correlations along x, MSE per edge (paper: 0.0077, 0.0078, 0.0078, 0.0125)
truerho <- x / 1.4
op <- par(mfrow = c(1, 4), mar = c(4.5, 4.5, 2, 1))
for (e in list(c(1, 2), c(2, 3), c(3, 4), c(4, 5))) {
  est <- partial_correlation_path(A, e[1], e[2]); o <- order(x)
  plot(x[o], truerho[o], type = "l", xlab = "x", ylab = bquote(rho[.(paste0(e[1], e[2]))]),
       main = sprintf("MSE %.4f", mean((est - truerho)^2)))
  lines(lowess(x, est, f = 0.4), col = "red")
}
par(op)

## Figure 6A: BM13 (LSGraph), as in the script
set.seed(2022)
bm <- LSGraph(Y = Y, X0 = X, bG = 3.0, DG = diag(5), k = 0, which = 1, m = 1, M = 1, Kmat = diag(2),
              c = 300, alpha = 0.5, beta = 0.5, gamma0 = 0, Adj0 = diag(5), burnin = 1000, nmc = 2000,
              verbose = FALSE)
plot_edge_heatmap(bm$Adj.mean, cutoff = 0.5, main = "(A) BM13")

## Figure 6B: mSSL, if the GitHub package is installed
## remotes::install_github("skdeshpande91/multivariate_SSL", subdir = "mSSL")
if (requireNamespace("mSSL", quietly = TRUE)) {
  set.seed(2022)
  ms <- fit_mSSL(X, Y)
  plot_edge_heatmap(ms$abs_partial_cor, cutoff = NULL, main = "(B) mSSL |partial correlation|")
}
