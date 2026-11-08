test_that("LSGraph (BM13) runs on Example 2 and returns edge frequencies", {
  sim <- simulate_example2(2020)
  set.seed(2022)
  fit <- LSGraph(Y = sim$Y, X0 = sim$X, bG = 3, DG = diag(5), k = 0, which = 1, m = 1, M = 1,
                 Kmat = diag(2), c = 300, alpha = 0.5, beta = 0.5, gamma0 = 0, Adj0 = diag(5),
                 burnin = 20, nmc = 40, verbose = FALSE)
  expect_equal(dim(fit$Adj.mean), c(5L, 5L))
  expect_true(all(fit$Adj.mean >= 0 & fit$Adj.mean <= 1))
  expect_equal(fit$accept + fit$reject + fit$nondecomp, 40)
  expect_length(fit$gamma.mean, 1)
})

test_that("MCMCGraph, HIWsim and meanOmegaHIW run", {
  sim <- simulate_example2(2020)
  set.seed(1)
  g <- MCMCGraph(Y = sim$Y, bG = 3, DG = diag(5), Adj0 = diag(5), burnin = 10, nmc = 20)
  expect_equal(dim(g$Adj.est), c(5L, 5L))
  Adj <- diag(5); Adj[1, 2] <- Adj[2, 1] <- 1; Adj[2, 3] <- Adj[3, 2] <- 1
  h <- HIWsim(Adj, 3, diag(5))
  expect_true(isSymmetric(unname(h$Sigma), tol = 1e-8))
  expect_true(all(eigen(h$Sigma)$values > 0))
  m <- meanOmegaHIW(Adj, 3, diag(5))
  expect_equal(dim(m), c(5L, 5L))
})

test_that("JRNS samplers run and the edge-frequency helper works", {
  sim <- simulate_example3(2022)
  Y <- sim$Y[1:100, 1:6]; X <- sim$X[1:100, 1:3]
  set.seed(3)
  fit <- JRNS(B_ini = matrix(0, 3, 6), Omega_ini = diag(6), Y = Y, X = X, nmc = 20, burnin = 10)
  expect_equal(dim(fit$Omega), c(6L, 6L, 21L))
  expect_equal(dim(fit$Bhat), c(3L, 6L))
  f <- jrns_edge_frequencies(fit)
  expect_true(isSymmetric(f)); expect_equal(diag(f), rep(0, 6))
  sb <- StepB(B_ini = matrix(0, 3, 6), sigma2_ini = rep(1, 6), Y = Y, X = X, nmc = 10, burnin = 5)
  expect_equal(dim(sb$Bhat), c(3L, 6L))
  so <- StepOmega(S = crossprod(Y) / 100, temp = diag(6), nmc = 10, burnin = 5, n = 100)
  expect_equal(dim(so$Omegahat), c(6L, 6L))
})

test_that("jrns_initial_values reproduces the script's initialisation shape", {
  skip_if_not_installed("glmnet"); skip_if_not_installed("glasso")
  sim <- simulate_example3(2022)
  init <- jrns_initial_values(sim$X[1:150, ], sim$Y[1:150, 1:8])
  expect_equal(dim(init$B), c(10L, 8L))
  expect_equal(dim(init$Omega), c(8L, 8L))
})

test_that("k-lasso helpers run on a few samples", {
  skip_if_not_installed("huge")
  sim <- simulate_example1(2020)
  idx <- seq(1, 300, by = 30)
  K <- klasso(sim$Y[idx, ], sim$X[idx, , drop = FALSE], h_grid = c(0.5, 1), l_grid = c(0.3, 0.8), verbose = FALSE)
  expect_equal(dim(K), c(3L, 3L, length(idx)))
  km <- kernel_mc(sim$Y[idx, ], sim$X[idx, , drop = FALSE], h = 0.5)
  expect_equal(dim(km$Sig), c(length(idx), 3L, 3L))
  expect_equal(ekernel(c(0, 0), c(0, 0), 1), dnorm(0))
})

test_that("fit_mSSL fails cleanly when the GitHub package is absent", {
  skip_if(requireNamespace("mSSL", quietly = TRUE))
  expect_error(fit_mSSL(matrix(rnorm(20), 10), matrix(rnorm(30), 10)), "mSSL")
})
