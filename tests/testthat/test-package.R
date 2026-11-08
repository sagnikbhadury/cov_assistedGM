test_that("breast cancer data have the dimensions stated in the paper", {
  data(brca_rna); data(brca_rppa)
  expect_equal(dim(brca_rna), c(883L, 3L))
  expect_equal(dim(brca_rppa), c(883L, 15L))
  expect_equal(colnames(brca_rna), c("ERBB2", "ESR1", "PGR"))
  expect_equal(unname(round(colMeans(brca_rppa), 8)), rep(0, 15))
  expect_equal(unname(round(apply(brca_rppa, 2, sd), 8)), rep(1, 15))
})

test_that("simulators reproduce the scripts' first draws", {
  s1 <- simulate_example1(2020)
  expect_equal(dim(s1$Y), c(300L, 3L))
  expect_true(all(diff(s1$x) >= 0))
  s2 <- simulate_example2(2020)
  expect_equal(dim(s2$Y), c(300L, 5L))
  expect_equal(sum(s2$G0) - 5, 8)
})

test_that("PxG_PL runs, returns labels in 1:N and the Gsum array", {
  data(brca_rna); data(brca_rppa)
  set.seed(1)
  fit <- PxG_PL(YData = brca_rppa[1:120, ], XData = brca_rna[1:120, ],
                spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
                mu0 = rep(0, 3), sigsq0 = 1, b1 = 30, b2 = 0.1,
                M = 10, N = 3, burnin = 5, nmcmc = 20, z0 = sample(1:3, 120, TRUE), B = 1,
                verbose = FALSE)
  expect_length(fit$Label, 120)
  expect_true(all(fit$Label %in% 1:3))
  expect_equal(dim(fit$Gsum), c(14L, 15L, 120L))
  expect_true(all(fit$Gsum >= 0 & fit$Gsum <= 20))
  expect_true(is.finite(fit$DIC))
  P <- cluster_edge_probabilities(fit, nmcmc = 20, node_names = colnames(brca_rppa))
  expect_true(all(sapply(P, function(M) isSymmetric(M) && all(M >= 0 & M <= 1))))
})

test_that("comparison PL models and Graph_NS run", {
  data(brca_rna); data(brca_rppa)
  set.seed(2)
  co <- CovOnly_PL(YData = brca_rppa[1:80, ], XData = brca_rna[1:80, ],
                   spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
                   mu0 = rep(0, 3), sigsq0 = 1, b1 = 30, b2 = 0.1,
                   M = 10, N = 3, burnin = 5, nmcmc = 10, z0 = sample(1:3, 80, TRUE), verbose = FALSE)
  expect_true(is.finite(co$DIC))
  go <- GraphOnly_PL(YData = brca_rppa[1:80, ], XData = brca_rna[1:80, ],
                     spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
                     mu0 = rep(0, 3), sigsq0 = 1, b1 = 30, b2 = 0.1,
                     M = 10, N = 3, burnin = 5, nmcmc = 10, z0 = sample(1:3, 80, TRUE), B = 1, verbose = FALSE)
  expect_true(is.finite(go$DIC))
  G <- Graph_NS(YData = brca_rppa[1:80, ], spike = 0.1, slab = 5, a1 = 1, a2 = 1, alphaG = 0.5,
                burnin = 5, nmcmc = 20, verbose = FALSE)
  expect_equal(dim(G), c(15L, 15L))
  S <- symmetrize_graph(G)
  expect_true(isSymmetric(S))
  expect_equal(diag(S), rep(1, 15))
})

test_that("PxG_GW runs on the Example 1 data", {
  skip_if_not_installed("BDgraph"); skip_if_not_installed("invgamma")
  sim <- simulate_example1(2020)
  set.seed(2020)
  fit <- PxG_GW(Y = sim$Y, X = sim$X, b = 3, D = diag(3), alphaG = 0.1,
                mu0 = 0, sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5,
                burnin = 2, nmcmc = 5, z0 = sample(1:5, 300, TRUE), B = 1, verbose = FALSE)
  expect_length(fit$Label, 300)
  expect_equal(dim(fit$Gsum), c(3L, 3L, 300L))
  expect_true(all(fit$Gsum >= 0 & fit$Gsum <= 1))
  expect_length(partial_correlation_path(fit, 1, 2), 300)
})

test_that("within_cluster_ss reproduces the script's weighted quantity", {
  X <- matrix(rnorm(40), 20, 2); lab <- rep(1:2, each = 10)
  w <- within_cluster_ss(X, lab)
  m1 <- colMeans(X[1:10, ]); m2 <- colMeans(X[11:20, ])
  ss1 <- sum(sweep(X[1:10, ], 2, m1)^2); ss2 <- sum(sweep(X[11:20, ], 2, m2)^2)
  expect_equal(w$ss, ss1 + ss2)
  expect_equal(w$weighted_ss, 10 * ss1 + 10 * ss2)
})
