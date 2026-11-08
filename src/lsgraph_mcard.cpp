// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file mcard.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>

#include "lsgraph_graph_utility.h"
#include "lsgraph_graphfun.h"

// maximum cardinality search of a graph given adjacency matrix Adj on p nodes 
// here we make sure to set the diagonal entries of Adj to 1's as per our convention

bool mcard(arma::umat& Adj){
  
  const unsigned int p = Adj.n_rows;
  // make sure the diagonal entries of Adj is 1
  // Adj should only contains 0 and 1
  Adj = Adj - diagmat(diagvec(Adj)) + arma::eye< arma::umat >(p, p);
  bool isdecomp = 1;
  
  arma::uvec Nnodes(1); // 'Nnodes' stores common indices, starts from 1
  Nnodes(0) = 1;        // 'Nnodes' contains nodes which have been check already
  
  unsigned int nnum = 1;
  
  arma::uvec UNnodes(p-1); // 'UNnodes' stores common indices, starts from 1
                           // 'UNnodes' contains nodes which haven't been check yet
  for(unsigned int i=0; i<(p-1); i++) { UNnodes(i) = i + 2; } // initialize 2:p

  arma::uvec nn, maxj, non0, nej;
  arma::urowvec rowj;
  unsigned int j;

  while( isdecomp && (nnum<p) ){
    // cout<<"nnum = "<<nnum<<endl; /***/
    nn = sum( Adj.submat(UNnodes-1, Nnodes-1), 1 ); //indices start from 0 for Adj
    // cout<<"nn : "<<endl<<nn<<endl; /***/

    j = nn.index_max(); // needs armadillo 7.8 or can use next two lines instead
    // maxj = find(nn==nn.max());
    // j = maxj(0); // index j is C++ matrix index, starts from 0

    Nnodes = append(Nnodes, UNnodes(j)); // 'UNnodes' takes C++ matrix index
    // cout<<"Nnodes : "<<endl<<Nnodes<<endl; /***/

    UNnodes = erase(UNnodes, j); // 'UNnodes' takes C++ matrix index
    // cout<<"UNnodes : "<<endl<<UNnodes<<endl; /***/

    j = Nnodes(Nnodes.n_elem-1); // j is the last index in 'Nnodes', now j is common index, starts from 1
    // cout<<"j = "<<j<<endl; /***/

    rowj = Adj.row(j-1); // Adj using index starts from 0
    non0 = find(rowj!=0)+1; // find out which one in 'rowj' is non-zero, but restore index to common index, starts from 1

    nej = vintersect(Nnodes, non0); // now, 'Nnodes' and 'non0' both have common indices, starts from 1
    // cout<<"nej : "<<endl<<nej<<endl;  /***/
    isdecomp = all(all(Adj(nej-1,nej-1))); // Adj using index starts from 0

    nnum += 1;
  }

  return isdecomp;
}



