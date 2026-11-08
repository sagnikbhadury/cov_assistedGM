// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file graph_utility.h of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see src/README_sources.md.
#ifndef GUARD_GRAPH_UTILITY_h
#define GUARD_GRAPH_UTILITY_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>

#include <vector>
using std::vector;

arma::uvec append(arma::uvec vp, unsigned int last);
arma::uvec erase(arma::uvec vp, unsigned int index);


// armadillo vector has a nested typedef 'value_type'
template< typename ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE vintersect( ARMA_VECTOR_TYPE first, ARMA_VECTOR_TYPE second ){

  first = sort(first);
  second = sort(second);
  vector< typename ARMA_VECTOR_TYPE::value_type > output;
  set_intersection( first.begin(), first.end(), second.begin(), second.end(), back_inserter(output) );

  ARMA_VECTOR_TYPE result = arma::conv_to< ARMA_VECTOR_TYPE >::from(output);
  return result;
}
// reference[http://www.cplusplus.com/forum/general/70687/]

#endif


// sort() function may be reloaded so not need for #include <algorithm>
// same thing to set_intersection() and back_inserter() in <vector>

