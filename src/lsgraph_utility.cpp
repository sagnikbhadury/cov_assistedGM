// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file utility.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see src/README_sources.md.
// #define MATHLIB_STANDALONE // put it on top
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>
#include <cmath> // for log(), M_PI, lgamma()
#include <vector>
using std::vector;

#include "lsgraph_utility.h"

double log_multi_gamma(unsigned int p, double n){
	double f;
	f = 0.25*p*(p-1)*log(M_PI);

	for(unsigned int i=1; i<=p; i++){
		f += lgamma(n + 0.5 - 0.5*i); // pay attention to type conversion
	}

	return f;
}


double log_pmargin(arma::mat Sgamma, unsigned int df, 
	double bG, arma::mat DG, unsigned int n){

	double f;
	f = 0.5 * n * df * log(2) + \
		0.5 * (bG + df -1) * log(det(DG)) - \
		0.5 * (bG + n + df -1) * log(det(DG + Sgamma)) + \
		log_multi_gamma(df, 0.5*(bG + n + df -1)) - \
		log_multi_gamma(df, 0.5*(bG + df -1));

	return f;
}


double log_margin(arma::mat Sgamma, vector<arma::uvec> C, vector<arma::uvec> S, 
	double bG, arma::mat DG, unsigned int p_gamma, unsigned int n, unsigned int k, double c){

	unsigned int q = Sgamma.n_rows;
	unsigned int numberofcliques = static_cast<unsigned int>( C.size() ) - 1;

	double f;
	f = -0.5 * n * q * log(2*M_PI) - \
		0.5 * p_gamma * (k+1) * q * log(c+1) + \
		log_pmargin(Sgamma(C[1]-1,C[1]-1), (C[1]).n_elem, bG, DG(C[1]-1,C[1]-1), n);

	if(numberofcliques>1){
		for(unsigned int i=2; i<=numberofcliques; i++){
			f = f + log_pmargin(Sgamma(C[i]-1,C[i]-1), (C[i]).n_elem, bG, DG(C[i]-1,C[i]-1), n) - \
			 log_pmargin(Sgamma(S[i]-1,S[i]-1), (S[i]).n_elem, bG, DG(S[i]-1,S[i]-1), n);
		}
	}

	return f;
}


unsigned int DiscUnif(unsigned int L){

	double sum = 1.0/L;
	double u = unif_rand();

	unsigned int i = 1;
	while(sum<u){
		i += 1;
		sum += 1.0/L;
	}
	return i;
}
// set_seed(23454,123);


unsigned int armaDisc(unsigned int L){

	double sum = 1.0/L;
	arma::vec u = arma::randu<arma::vec>(1);

	unsigned int i = 1;
	while(sum<u(0)){
		i += 1;
		sum += 1.0/L;
	}
	return i;
}


// k is number of knots
arma::mat makesplinematrix(arma::mat X, unsigned int k, unsigned int which, double m, double M, arma::mat Kmat){

	if(k == 0) { return X; }

	unsigned int n = X.n_rows;
	unsigned int p = X.n_cols;
	arma::mat Xspline(n, (k+1)*p, arma::fill::zeros);
	Xspline.cols(0,p-1) = X;

	if(which == 0){
		double h = (M-m)/(k+1);
		arma::vec theta = m + h*arma::linspace<arma::vec>(1, k, k);
		arma::mat Xtemp(n,p);
		for(unsigned int i=0; i<k; i++){
			Xtemp = X - theta(i);
			Xtemp(find(Xtemp<0)).zeros();
			Xspline.cols((i+1)*p, (i+2)*p-1) = Xtemp;
		}
	}

	else if(which == 1){
		m = X.min();
		M = X.max();
		double h = (M-m)/(k+1);
		arma::vec theta = m + h*arma::linspace<arma::vec>(1, k, k);
		arma::mat Xtemp(n,p);
		for(unsigned int i=0; i<k; i++){
			Xtemp = X - theta(i);
			Xtemp(find(Xtemp<0)).zeros();
			Xspline.cols((i+1)*p, (i+2)*p-1) = Xtemp;
		}
	}

	else if(which == 2){
		double h;
		arma::mat THETA(k, p);
		arma::mat Xtemp(n, p);

		for(unsigned int i=0; i<p; i++){
			m = min(X.col(i));
			M = max(X.col(i));
			h = (M - m)/(k+1);
			THETA.col(i) = m + h*arma::linspace<arma::vec>(1, k, k);
		}

		for(unsigned int i=0; i<k; i++){
			Xtemp = X;
			for(unsigned int j=0; j<p; j++){
				Xtemp.col(j) = Xtemp.col(j) - THETA(i, j);
			}
			Xtemp(find(Xtemp<0)).zeros();
			Xspline.cols((i+1)*p, (i+2)*p-1) = Xtemp;
		}
	}

	else if(which == 3){
		if(Kmat.n_rows != p){
			Rcpp::stop("error in Kmat");
		}
		if(Kmat.n_cols != k){
			Rcpp::stop("error in Kmat");
		} 

		arma::mat Xtemp(n, p, arma::fill::zeros);
		for(unsigned int i=0; i<k; i++){
			Xtemp = X;
			for(unsigned int j=0; j<n; j++){
				for(unsigned int l=0; l<p; l++){
					Xtemp(j, l) = Xtemp(j, l) - Kmat(l, i);
				}
			}
			Xtemp(find(Xtemp<0)).zeros();
			Xspline.cols((i+1)*p, (i+2)*p-1) = Xtemp;
		}
	}

	return Xspline;
}


arma::uvec getsplineind(arma::uvec ind, unsigned int p, unsigned int k){
	if(k == 0) { return ind; }

	unsigned int L = ind.n_elem;
	arma::uvec splineind((k+1)*L, arma::fill::zeros);

	for(unsigned int i=0; i<L; i++){
		for(unsigned int j=0; j<k+1; j++){
			splineind(i + j*L) = ind(i) + j*p;
		} 
	}

	return splineind;
}

