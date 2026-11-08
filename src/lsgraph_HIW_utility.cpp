// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file HIW_utility.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>
// must use armadillo > 7.7
// must use c++11
// must link to LAPACK and BLAS
#include "lsgraph_HIW_utility.h"

/*
 Generates 1 draw from a Wishart distribution - allows for singular Wishart too, 
 in cases that the distn parameter S is rank deficient of rank r.

 K is W_p(df,A) with sum-of-squares parameter S = A^{-1} and d.o.f. df
 Dimension p is implicit

 Usual nonsingular case: r=p<df; now allow for r<p  with r integral
 Note that  E(K)=df.S^{-1}  in this notation 

 Noticing pdf is p(K) = cons. |K|^((df-p-1)/2) exp(-trace(K S)/2) 
 in case usual modification 

 Returns matrix W of dimension p by p of rank r=rank(S)

 EXAMPLE: reference posterior in a normal model N(0,Sigma) with precision 
    	  matrix Omega = Sigma^{-1}
	      Random sample of size n has sample var matrix V=S/n with S=\sum_i x_ix_i'
    	  Ref posterior for precision mx Omega is W_p(n,A) with A=S^{-1}
    	  e.g., K = wishart_InvA_rnd(n,n*V); draw 1000 samples

 Useful for looking at posteriors for correlations and eigenvalues
 and also for posterior on graph - looking for elements of Omega near 0

 To deal with numerical problems when p is larger, or high collinearity
 and/or rank deficiency, use line 37 - 43
*/


arma::mat wishart_InvA_rnd(double df, arma::mat S){

	unsigned int p = S.n_cols;
	arma::mat K(p, p, arma::fill::zeros);

	arma::mat P, D, V; // we don't use V
	arma::vec d;
	svd(P, d, V, S);
	D = diagmat(1.0/sqrt(d));
	arma::uvec i = find(d>max(d)*1e-9);
	unsigned int r = i.n_elem;
	P = P.cols(i)* D.cols(i);

/* OR replace above 7 lines with the 2 lines below
	arma::mat P = chol(inv(S)).t();
	unsigned int r = p;
*/

	arma::mat U(p, r, arma::fill::zeros);
	arma::rowvec temp;
	//arma::vec rgamma0(1);
  double rgamma0;
	
	// Now generate W_p(df,A) using Bartlett decomposition
	for(unsigned int j=0; j<r; j++){
		// temp = arma::randn<arma::rowvec>(r-j);
		// rgamma0 = arma::randg<arma::vec>(1, arma::distr_param((df-j)/2.0, 2.0));
		// temp(0) = sqrt(rgamma0(0));
		
	  temp.zeros(r-j);
	  for(unsigned int h=0; h<(r-j); h++){ temp(h) = norm_rand(); }
		rgamma0 =	R::rgamma((df-j)/2.0, 2.0);
		temp(0) = sqrt(rgamma0);

		U.submat(j, j, j, r-1) = temp;
	}

	U = U * P.t();
	K = U.t() * U;

	return K;
}


arma::uvec arma_union(arma::uvec a, arma::uvec b){
	unsigned int La = a.n_elem;
	unsigned int Lb = b.n_elem;

	arma::uvec c(La+Lb, arma::fill::zeros);
	for(unsigned int i=0; i<La; i++) c(i) = a(i);
	for(unsigned int j=La; j<(La+Lb); j++) c(j) = b(j-La);

	c = unique(c);
	c = sort(c);

	return c;
}


// http://stackoverflow.com/questions/29724083/trying-to-write-a-setdiff-function-using-rcpparmadillo-gives-compilation-error
arma::uvec arma_setdiff(arma::uvec x, arma::uvec y){
    x = unique(x); // unique() also sorts the vector
    y = unique(y);
    arma::uvec z = x;
    arma::uvec q1;

    for (unsigned int j = 0; j < y.n_elem; j++){
        q1 = find(z == y(j));
        if (!q1.empty()) {
            z.shed_row(q1(0));
        }
    }

    return z;
}


/* RNG RNG RNG RNG RNG RNG RNG RNG RNG RNG RNG RNG RNG*/
arma::mat rMNorm(arma::vec m, arma::mat V){
	unsigned int p = m.n_elem;
  
	// arma::mat rnorm = m + (chol(V)).t() * arma::randn<arma::vec>(p);
	
	arma::vec normvec(p, arma::fill::zeros);
	for(unsigned int i=0; i<p; i++){
	  normvec(i) = norm_rand();
	}
	arma::mat rnorm = m + (chol(V)).t() * normvec;
	
	return rnorm;
}


