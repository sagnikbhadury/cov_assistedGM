// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file HIWsim.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>
#include <vector>
using std::vector;

#include "lsgraph_graphfun.h" // makedecompgraph()
#include "lsgraph_HIW_utility.h"


// Sample the HIW_G(bG,DG) distribution on a graph G with adjacency matrix Adj
//
// [[Rcpp::export]]
Rcpp::List HIWsim(const arma::umat& Adj, const double bG, const arma::mat& DG){

	GetRNGstate();

  arma::mat Sigma, Omega;
  
	unsigned int p = Adj.n_cols;
	Sigma.zeros(p, p);
	Omega.zeros(p, p);

	vector<arma::uvec> Cliques, Separators;
	makedecompgraph(Adj, Cliques, Separators);
	unsigned int numberofcliques = static_cast<unsigned int>( Cliques.size() ) - 1;


	/* ************************************************** */
	// Create some working arrays that are computed only once
	arma::mat C1 = inv(DG(Cliques[1]-1, Cliques[1]-1)/bG); // seems didn't use C1
	arma::uvec c1, UN, cid, sid, dif;
	c1 = Cliques[1];
	UN = c1;

	vector<arma::mat> DSi(numberofcliques+1), DRS(numberofcliques+1), mU(numberofcliques+1);
	arma::mat a0(1, 1, arma::fill::zeros);
	DSi[0] = a0; DSi[1] = a0;
	DRS[0] = a0; DSi[1] = a0;
	mU[0]  = a0; mU[1]  = a0;

	unsigned int sizedif;

	for(unsigned int i=2; i<=numberofcliques; i++){
		sid = Separators[i];
		DSi[i] = inv(DG(sid-1, sid-1));
		cid = Cliques[i];
		dif = arma_setdiff(cid, UN);
		UN = arma_union(cid, UN);
		sizedif = dif.n_elem;
		DRS[i] = DG(dif-1, dif-1) - DG(dif-1, sid-1) * DSi[i] * DG(sid-1, dif-1);
		DRS[i] = ( DRS[i] + (DRS[i]).t() )/2.0;
		mU[i] = DG(dif-1, sid-1) * DSi[i];
	}
	/* ************************************************** */


	// MC Sampling
	UN = c1;
	arma::mat Sigmaj(p, p, arma::fill::zeros);
	Sigmaj(c1-1, c1-1) = inv( wishart_InvA_rnd( bG+(Cliques[1]).n_elem-1, DG(Cliques[1]-1, Cliques[1]-1) ) );
	
	arma::mat SigRS, Ui;

	for(unsigned int i=2; i<=numberofcliques; i++){
		dif = arma_setdiff(Cliques[i], UN);
		UN = arma_union(Cliques[i], UN);
		sizedif = dif.n_elem;
		sid = Separators[i];
		
		SigRS = inv( wishart_InvA_rnd( bG+(Cliques[i]).n_elem-1, DRS[i] ) );
		Ui = rMNorm( vectorise((mU[i]).t(), 0), kron(SigRS, DSi[i]) );
		/* reshape may be dangerous !!! PAY ATTANTION !!! */
		Ui.reshape(Ui.n_elem/sizedif, sizedif);
		Sigmaj(dif-1, sid-1) = Ui.t() * Sigmaj(sid-1, sid-1);
		Sigmaj(sid-1, dif-1) = (Sigmaj(dif-1, sid-1)).t();
		Sigmaj(dif-1, dif-1) = SigRS + Sigmaj(dif-1, sid-1) * inv(Sigmaj(sid-1, sid-1)) * Sigmaj(sid-1, dif-1);
	}


	// Next, completion operation for sampled variance matrix
	arma::uvec H = c1;
	arma::uvec h;
	for(unsigned int i=2; i<=numberofcliques; i++){
		dif = arma_setdiff(Cliques[i], H);
		sid = Separators[i];
		h = arma_setdiff(H, sid);
		Sigmaj(dif-1, h-1) = Sigmaj(dif-1, sid-1) * inv(Sigmaj(sid-1, sid-1)) * Sigmaj(sid-1, h-1);
		Sigmaj(h-1, dif-1) = (Sigmaj(dif-1, h-1)).t();
		H = arma_union(H, Cliques[i]);
	}

	Sigma = Sigmaj;


	// Next, computing the corresponding sampled precision matrix
	arma::cube Caux(p, p, numberofcliques, arma::fill::zeros);
	arma::cube Saux = Caux;
	cid = Cliques[1];
	// Caux(cid-1, cid-1, 0) = inv(Sigmaj(cid-1, cid-1));
	(Caux.slice(0)).submat(cid-1, cid-1) = inv(Sigmaj(cid-1, cid-1));

	for(unsigned int i=2; i<=numberofcliques; i++){
		cid = Cliques[i];
		// Caux(cid-1, cid-1, i-1) = inv(Sigmaj(cid-1, cid-1));
		(Caux.slice(i-1)).submat(cid-1, cid-1) = inv(Sigmaj(cid-1, cid-1));
		sid = Separators[i];
		// Saux(sid-1, sid-1, i-1) = inv(Sigmaj(sid-1, sid-1));
		(Saux.slice(i-1)).submat(sid-1, sid-1) = inv(Sigmaj(sid-1, sid-1));
	}

	Omega = sum(Caux, 2) - sum(Saux, 2);

	PutRNGstate();
	
	return Rcpp::List::create(Rcpp::Named("Sigma") = Sigma,
                            Rcpp::Named("Omega") = Omega);

}


// posterior mean of Omega on decomposable graph G when its inverse Sigma~HIW(bG,DG)
// input:  G - the usual cell array whose 2 cells are the comps and seps
//
// [[Rcpp::export]]
arma::mat meanOmegaHIW(const arma::umat& Adj, const double bG, const arma::mat& DG){

	unsigned int p = Adj.n_cols;
	arma::mat Omega(p, p, arma::fill::zeros);

	vector<arma::uvec> Cliques, Separators;
	makedecompgraph(Adj, Cliques, Separators);
	unsigned int k = static_cast<unsigned int>( Cliques.size() ) - 1;

	arma::uvec j = Cliques[1];
	Omega(j-1,j-1) = (bG+(Cliques[1]).n_elem-1) * inv(DG(j-1, j-1)); // first component 

	for(unsigned int c=2; c<=k; c++){
		j = Separators[c]; // visit separator c
		Omega(j-1, j-1) = Omega(j-1, j-1) - (bG+(Separators[c]).n_elem-1) * inv(DG(j-1, j-1));
		j = Cliques[c]; // visit component c
		Omega(j-1, j-1) = Omega(j-1, j-1) + (bG+(Cliques[c]).n_elem-1) * inv(DG(j-1, j-1));
	}

	return Omega;
}


