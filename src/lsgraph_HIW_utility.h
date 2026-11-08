// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file HIW_utility.h of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see src/README_sources.md.
#ifndef HIW_UTILITY_h
#define HIW_UTILITY_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>

arma::mat wishart_InvA_rnd(double df, arma::mat S);

arma::uvec arma_union(arma::uvec a, arma::uvec b);

arma::uvec arma_setdiff(arma::uvec x, arma::uvec y);

arma::mat rMNorm(arma::vec m, arma::mat V);

Rcpp::List HIWsim(const arma::umat& Adj, const double bG, const arma::mat& DG);

arma::mat meanOmegaHIW(const arma::umat& Adj, const double bG, const arma::mat& DG);

#endif

