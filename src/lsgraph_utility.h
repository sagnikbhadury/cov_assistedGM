// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file utility.h of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
#ifndef GUARD_UTILITY_h
#define GUARD_UTILITY_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>
#include <vector>
using std::vector;


double log_multi_gamma(unsigned int p, double n);


double log_pmargin(arma::mat Sgamma, unsigned int df, 
	double bG, arma::mat DG, unsigned int n);


double log_margin(arma::mat Sgamma, vector<arma::uvec> C, vector<arma::uvec> S, 
	double bG, arma::mat DG, unsigned int p_gamma, unsigned int n, unsigned int k, double c);


unsigned int DiscUnif(unsigned int L);
unsigned int armaDisc(unsigned int L);


arma::mat makesplinematrix(arma::mat X, unsigned int k, unsigned int which, double m, double M, arma::mat Kmat);
arma::uvec getsplineind(arma::uvec ind, unsigned int p, unsigned int k);

Rcpp::List LSGraph( 
	const arma::mat& Y, const arma::mat& X0, const double bG, const arma::mat& DG, 
	const unsigned int k, const unsigned int which, const double m, const double M, const arma::mat& Kmat,
	const double c, const double alpha, const double beta,
	const arma::uvec gamma0, const arma::umat& Adj0,
	const unsigned int burnin, const unsigned int nmc, const bool verbose);

Rcpp::List MCMCGraph(
	const arma::mat& Y, const double bG, const arma::mat& DG, 
	const arma::umat& Adj0,
 	const unsigned int burnin, const unsigned int nmc);

#endif

