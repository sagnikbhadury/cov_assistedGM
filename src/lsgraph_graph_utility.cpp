// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file graph_utility.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo> // not necessary, already in "graph_utility.h"
#include "lsgraph_graph_utility.h"


arma::uvec append(arma::uvec vp, unsigned int last){
  const unsigned int L = vp.n_elem;
  arma::uvec comb(L+1);
  for(unsigned int i=0; i<L; i++) { comb(i) = vp(i); }
  comb(L) = last;
  
  return comb;
}


arma::uvec erase(arma::uvec vp, unsigned int index){
  const unsigned int L = vp.n_elem;
  arma::uvec after(L-1);
  unsigned int j = 0;
  for(unsigned int i=0; i<L; i++){
    if( i != index){
      after(j) = vp(i);
      j += 1;
    }
  }
  
  return after;
}

