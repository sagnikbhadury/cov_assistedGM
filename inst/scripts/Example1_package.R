### Section 4.1 "Example 1: piecewise linear precision matrix" of Niu et al. (2024), run with the
### covassistedGM package. Same models, arguments and seeds as the authors' Example1.R; the data
### block is simulate_example1(), the k-lasso loop is klasso(). The three samplers take a few
### minutes each; k-lasso about 100 minutes (load the authors' klasso.RData instead, see below).

library(covassistedGM)
sim <- simulate_example1(seed = 2020)
Y <- sim$Y; X <- sim$X; x <- sim$x

set.seed(2020)
A <- PxG_GW(Y = Y, X = X, b = 3, D = diag(3), alphaG = 0.1,
            mu0 = 0, sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5,
            burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 1)
A$DIC  # 2193.304

set.seed(2020)
B <- CovOnly_GW(Y = Y, X = X, b = 3, D = diag(3), alphaG = 0.1,
                mu0 = 0, sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5,
                burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 1)
B$DIC  # 2584.573

set.seed(2020)
C <- GraphOnly_GW(Y = Y, X = X, b = 3, D = diag(3), alphaG = 0.1,
                  mu0 = 0, sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5,
                  burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 1)
C$DIC  # 3189.946

## k-lasso (about 100 minutes). To skip, load the authors' result:
## load(system.file("original_supplement/RData/klasso.RData", package = "covassistedGM"))
set.seed(2020)
Omega_est <- klasso(Y, X, h_grid = seq(0.1, 1, length.out = 20), l_grid = seq(0.1, 1, length.out = 20))

## Figures 3 and 4 in base graphics (the authors use ggplot2/ggpubr; see vignette("pxg-simulations"))
xseq <- seq(-1, 1, by = 0.01)
interval <- function(v) if (v < -1/3) 1 else if (v > 1/3) 3 else 2
op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 2, 1))
for (e in list(c(1, 2), c(1, 3), c(2, 3))) {
  plot(xseq, sapply(xseq, function(v) sim$G0(interval(v), v)[e[1], e[2]]), type = "l", ylim = c(0, 1),
       xlab = "x", ylab = bquote(p[.(paste0(e[1], e[2]))]))
  lines(lowess(x, A$Gsum[e[1], e[2], ], f = 0.3), col = "red")
}
for (e in list(c(1, 2), c(1, 3), c(2, 3))) {
  tr <- sapply(seq_along(x), function(s) sim$Omega(interval(x[s]), x[s])[e[1], e[2]] / 1.2)
  pxg <- partial_correlation_path(A, e[1], e[2])
  kl <- Omega_est[e[1], e[2], ] / sqrt(Omega_est[e[1], e[1], ] * Omega_est[e[2], e[2], ])
  plot(x, tr, type = "l", ylim = c(-0.2, 1), xlab = "x", ylab = bquote(rho[.(paste0(e[1], e[2]))]),
       main = sprintf("MSE PxG %.4f, k-lasso %.4f", mean((pxg - tr)^2), mean((kl - tr)^2)))
  lines(lowess(x, pxg, f = 0.3), col = "red"); lines(lowess(x, kl, f = 0.3), col = "blue")
}
par(op)
