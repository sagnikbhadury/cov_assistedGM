// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file graphfun.h of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
#ifndef GUARD_GRAPHFUN_h
#define GUARD_GRAPHFUN_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>


#include <vector>
using std::vector;
// not need to #include "lsgraph_graph_utility.h" here

bool mcard(arma::umat& Adj);
void makedecompgraph(const arma::umat& Adj0, vector< arma::uvec >& Cliques, vector< arma::uvec >& Separators);

#endif


