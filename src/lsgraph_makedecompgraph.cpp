// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file makedecompgraph.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see src/README_sources.md.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>

#include <vector>
using std::vector;
#include "lsgraph_graph_utility.h"
#include "lsgraph_graphfun.h"


void makedecompgraph(const arma::umat& Adj0, vector< arma::uvec >& Cliques, vector< arma::uvec >& Separators){
    Cliques.clear();
    Separators.clear();
    Cliques.shrink_to_fit();    // require c++11
    Separators.shrink_to_fit(); // require c++11

    arma::umat Adj = Adj0;
    const unsigned int p = Adj.n_rows;

    arma::uvec order(p);
    for(unsigned int j=0; j<p; j++) { order(j) = j + 1; } // order = 1:p

    unsigned int i = 1;
    unsigned int b;
    arma::urowvec nn;

    while(i<p){
        nn = sum( Adj.submat(0, i, i-1, p-1), 0 ); // Adj(1:i, i+1:p) -> Adj(0:i-1, i:p-1)
        b = nn.index_max(); // 'b' is C++ index, starts from 0, b+1 is the real index
        //cout<<b+1<<endl; /***/
        order.swap_rows(i, b+i); // you want to swap i+1 and (b+1)+i -> here do i and b+i
        //cout<<order<<endl; /***/
        i += 1;
        Adj = Adj0;
        Adj = Adj.submat(order-1, order-1); // 'order' is 1:p, but we need 0:p-1 so use 'order - 1'
        //cout<<Adj<<endl; /***/
    }

    unsigned int numberofcliques = 1;

    arma::uvec elementC(1);
    elementC(0) = 0;

    Cliques.push_back(elementC);
    Separators.push_back(elementC);
    Separators.push_back(elementC);

    elementC(0) = 1;
    Cliques.push_back(elementC);

    i = 2;
    arma::uvec ind(1), indsub, temp;

    while(i<=p){
        ind(0) = i;
        if( accu(Adj.submat(ind-1, Cliques[numberofcliques]-1)) == Cliques[numberofcliques].n_elem ){ // index - 1
            Cliques[numberofcliques] = append(Cliques[numberofcliques], i);
        }
        else{
            numberofcliques += 1;
            indsub = find( Adj.submat(i-1, 0, i-1, i-1)==1 ) + 1; // + 1 to be real index
            temp = find(indsub == i);
            if(temp.is_empty() == 1) { indsub = append(indsub, i); }
            Cliques.push_back(indsub);
        }
        i += 1;
    }

    for(unsigned int k=1; k<=numberofcliques; k++){
        Cliques[k] = sort( order(Cliques[k]-1) ); // combine the following two commands
        // Cliques[k] = order(Cliques[k]-1); // index - 1      
        // Cliques[k] = sort(Cliques[k]);
    }

    arma::uvec UN = Cliques[1], intersect, UNtemp;
    unsigned int length;

    for(unsigned l=2; l<=numberofcliques; l++){
        intersect = sort( vintersect(UN, Cliques[l]) );  // combine the following two commands
        // intersect = vintersect(UN, Cliques[l]);
        // intersect = sort(intersect);
        Separators.push_back(intersect);

        UNtemp = UN;
        length = (Cliques[l]).n_elem;
        for(unsigned int m=0; m<length; m++){

            temp = find(UNtemp == Cliques[l](m));
            if(temp.is_empty() == 1) { UN = append( UN, Cliques[l](m) ); }
        }
    }

}
