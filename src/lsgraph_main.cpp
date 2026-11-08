// BM13 competitor (Bhadra & Mallick 2013; Niu, Pati & Mallick 2020): file LSGraph.cpp of the authors' LSGraph 0.1
// package, shipped in the JASA supplement of Niu et al. (2024). Packaging changes only, see inst/README_sources.md.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo>
#include <cmath> // log()
#include <vector>
using std::vector;

#include "lsgraph_graphfun.h"
#include "lsgraph_utility.h"


// [[Rcpp::export]]
Rcpp::List LSGraph( 
	const arma::mat& Y, const arma::mat& X0, const double bG, const arma::mat& DG, 
	const unsigned int k, const unsigned int which, const double m, const double M, const arma::mat& Kmat,
	const double c, const double alpha, const double beta,
	const arma::uvec gamma0, const arma::umat& Adj0,
	const unsigned int burnin, const unsigned int nmc, const bool verbose = true){

	// arma::arma_rng::set_seed(armaseed);
	GetRNGstate();

	arma::uvec gamma = gamma0;
	arma::umat Adj = Adj0;
	arma::mat X;

	unsigned int n = Y.n_rows;
	unsigned int p = X0.n_cols;
	unsigned int q = Y.n_cols;
	arma::mat In = arma::eye<arma::mat>(n,n);

	if(k == 0){
		X = X0;
	}
	else{
		X = makesplinematrix(X0, k, which, m, M, Kmat);
	}

	arma::uvec nonzeroind = find(gamma > 0);
	unsigned int nvar = nonzeroind.n_elem;
	arma::mat Xgamma, Px, Sgamma, Xgamma_star, Px_star, Sgamma_star;
	arma::uvec nonzeroind_spline;

	if(nvar == 0){
		Sgamma = Y.t() * Y;
	}
	else{

		if(k == 0){
			Xgamma = X.cols(nonzeroind);
		}
		else{
			nonzeroind_spline = getsplineind(nonzeroind, p, k);
			Xgamma = X.cols(nonzeroind_spline);
		}

		Px = Xgamma * inv(Xgamma.t() * Xgamma) * Xgamma.t();
		Sgamma = Y.t() * (In - c/(c+1)*Px) * Y;
	}


	vector<arma::uvec> Cliques, Separators, Cliques_star, Separators_star;
	makedecompgraph(Adj, Cliques, Separators);

	arma::umat indmat(q,q);
	for(unsigned int i = 0; i<q*q; i++) indmat(i)=i;
	indmat = trimatu(indmat);
	indmat.diag().zeros();
	arma::uvec upperind = find(indmat);
	unsigned int Lupperind = upperind.n_elem;

	double delta = 2.0/(double)(q-1);
	arma::umat gamma_mcmc(p, nmc, arma::fill::zeros);
	arma::ucube Adj_mcmc(q, q, nmc, arma::fill::zeros);

	
	int accept = 0;
	int reject = 0;
	int nondecomp = 0;

	arma::uvec gamma_map(p, arma::fill::zeros);
	arma::vec gamma_mean(p, arma::fill::zeros);
	arma::umat Adj_map(q, q, arma::fill::zeros);
	arma::mat Adj_mean(q, q, arma::fill::zeros); 
	double MAPlike = 0.0;
	arma::vec post(nmc, arma::fill::zeros);

	/****************************************************/
	arma::uvec gamma_star, nonzeroind_star;
	unsigned int changeind, nvar_star, changeind_g, indg;
	arma::vec u(1);
	double ratio, logmargin, logmargin_star, r, prop, r_g, logmarging, logmarging_star;
	bool sign;
	// double deter;
	arma::umat Adj_star;

	for(unsigned int iter=1; iter<=(burnin+nmc); iter++){

		/*********** M-H for gamma ***********/
		gamma_star = gamma;
		// changeind = armaDisc(p) - 1;
		changeind = DiscUnif(p) - 1;

		if( gamma(changeind) == 0 ){
			if(unif_rand()<alpha){
				gamma_star(changeind) = 1;
				ratio = (double)(nvar + 1)/(double)(p - nvar) * (1 - alpha)/alpha;
				sign = 1;
			}
			else{
				sign = 0;
			}
		}

		else{
			if(unif_rand()<(1-alpha)){
				gamma_star(changeind) = 0;
				ratio = (double)(p - nvar + 1)/(double)(nvar) * alpha/(1 - alpha);
				sign = 1;
			}
			else{
				sign = 0;
			}
		}


		if(sign == 1){
			nonzeroind_star = find(gamma_star > 0);
			nvar_star = nonzeroind_star.n_elem;
			if(nvar_star > 0){

				if(k == 0){
					Xgamma_star = X.cols(nonzeroind_star);
				}
				else{
					nonzeroind_spline = getsplineind(nonzeroind_star, p, k);
					Xgamma_star = X.cols(nonzeroind_spline);
				}

				Px_star = Xgamma_star * inv(Xgamma_star.t() * Xgamma_star) * Xgamma_star.t();
				Sgamma_star = Y.t() * (In - c/(c+1)*Px_star) * Y;

				// deter = rcond(Px_star);
				// if(deter<1e-18) { cout<<"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING"<<endl; )
			}
			else{
				Sgamma_star = Y.t() * Y;
			}

			logmargin = log_margin(Sgamma, Cliques, Separators, bG, DG, nvar, n, k, c);
			logmargin_star = log_margin(Sgamma_star, Cliques, Separators, bG, DG, nvar_star, n, k, c);
			r = logmargin_star - logmargin + log(ratio);

			if(log(unif_rand())<r){
				gamma = gamma_star;
				Sgamma = Sgamma_star;
				nvar = nvar_star;
				if(iter > burnin){
					gamma_mcmc.col(iter-burnin-1) = gamma;
				}
				if(verbose) Rcpp::Rcout<<"accept gamma"<<std::endl;
				if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
			}
			else{
				if(iter > burnin){
					gamma_mcmc.col(iter-burnin-1) = gamma;
				}
				if(verbose) Rcpp::Rcout<<"reject gamma"<<std::endl;
				if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
			}
		}

		else{
			if(iter > burnin){
				gamma_mcmc.col(iter-burnin-1) = gamma;
			}
			if(verbose) Rcpp::Rcout<<"same gamma"<<std::endl;
			if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
		}


		/*********** M-H and Stochastic Search for graphs ***********/
		Adj_star = Adj;
		// changeind_g = armaDisc(Lupperind) - 1;
		changeind_g = DiscUnif(Lupperind) - 1;
		indg = upperind(changeind_g);
		
		if(Adj(indg) == 0){
			if(unif_rand()<beta){
				Adj_star(indg) = 1;
				Adj_star = Adj_star.t();
				Adj_star(indg) = 1;
				prop = delta/(1-delta)*(1-beta)/beta;
				sign = 1;
			}
			else{
				sign = 0;
			}
		}

		else{
			if(unif_rand()<(1-beta)){
				Adj_star(indg) = 0;
				Adj_star = Adj_star.t();
				Adj_star(indg) = 0;
				prop = (1-delta)/delta*beta/(1-beta);
				sign = 1;				
			}
			else{
				sign = 0;
			}
		}

		if(sign == 1){
			if(mcard(Adj_star) == 1){
				makedecompgraph(Adj_star, Cliques_star, Separators_star);
				logmarging = log_margin(Sgamma, Cliques, Separators, bG, DG, nvar, n, k, c);
				logmarging_star = log_margin(Sgamma, Cliques_star, Separators_star, bG, DG, nvar, n, k, c);
				r_g = logmarging_star - logmarging + log(prop);

				if(log(unif_rand())<r_g){
					Adj = Adj_star;
					Cliques = Cliques_star;
					Separators = Separators_star;
					if(iter>burnin){
						Adj_mcmc.slice(iter-burnin-1) = Adj;
						post(iter-burnin-1) = logmarging_star;
						accept += 1;
					}
					if(verbose) Rcpp::Rcout<<"accept G"<<std::endl;
					if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
				}
				else{
					if(iter>burnin){
						Adj_mcmc.slice(iter-burnin-1) = Adj;
						post(iter-burnin-1) = logmarging;
						reject = reject + 1;
					}
					if(verbose) Rcpp::Rcout<<"reject G"<<std::endl;
					if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
				}
			}

			else{
				if(iter>burnin){
					Adj_mcmc.slice(iter-burnin-1) = Adj;
					post(iter-burnin-1) = logmarging;
					nondecomp += 1;
				}
				if(verbose) Rcpp::Rcout<<"reject nondecomp"<<std::endl;
				if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
			}
		}


		else{
			if(iter>burnin){
				Adj_mcmc.slice(iter-burnin-1) = Adj;
				post(iter-burnin-1) = logmarging;
				nondecomp += 1;
			}
			if(verbose) Rcpp::Rcout<<"reject nondecomp"<<std::endl;
			if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<std::endl;
		}
	}

	unsigned int mm = index_max(post);
	gamma_map = gamma_mcmc.col(mm);
	Adj_map   = Adj_mcmc.slice(mm);
	MAPlike = post(mm);
	gamma_mean = arma::conv_to<arma::vec>::from(sum(gamma_mcmc, 1))*(1.0/nmc);
	// Adj_mean   = arma::conv_to<arma::mat>::from(sum(Adj_mcmc, 2))*(1.0/nmc);
	arma::umat Temp_Adj(q, q, arma::fill::zeros);
	for(unsigned int i=0; i<nmc; i++){
		Temp_Adj = Temp_Adj + Adj_mcmc.slice(i);
	}
	Adj_mean = arma::conv_to<arma::mat>::from(Temp_Adj)*(1.0/nmc);

	PutRNGstate();
	
	return Rcpp::List::create(Rcpp::Named("accept") = accept,
                            Rcpp::Named("reject") = reject,
                            Rcpp::Named("nondecomp") = nondecomp,
                            Rcpp::Named("gamma.map") = gamma_map,
                            Rcpp::Named("gamma.mean") = gamma_mean,
                            Rcpp::Named("Adj.map") = Adj_map,
                            Rcpp::Named("Adj.mean") = Adj_mean,
                            Rcpp::Named("MAP.like") = MAPlike);

}



// [[Rcpp::export]]
Rcpp::List MCMCGraph(
	const arma::mat& Y, const double bG, const arma::mat& DG, 
	const arma::umat& Adj0,
 	const unsigned int burnin, const unsigned int nmc){

	GetRNGstate();

	unsigned int n = Y.n_rows;
	unsigned int q = Y.n_cols;
	arma::umat Adj = Adj0;

	arma::ucube Adj_mcmc(q, q, nmc, arma::fill::zeros);
	arma::mat Sgamma = Y.t() * Y;

	vector<arma::uvec> Cliques, Separators, Cliques_star, Separators_star;
	makedecompgraph(Adj, Cliques, Separators);

	arma::umat indmat(q,q);
	for(unsigned int i = 0; i<q*q; i++) indmat(i)=i;
	indmat = trimatu(indmat);
	indmat.diag().zeros();
	arma::uvec upperind = find(indmat);
	unsigned int Lupperind = upperind.n_elem;

	double delta = 2.0/(double)(q-1);

	unsigned int changeind_g, indg;
	double prop, r_g, logmarging, logmarging_star;
	logmarging = 0;
	arma::umat Adj_star;
	arma::vec post(nmc, arma::fill::zeros);

	for(unsigned int iter=1; iter<=(burnin+nmc); iter++){

		Adj_star = Adj;
		// changeind_g = armaDisc(Lupperind) - 1;
		changeind_g = DiscUnif(Lupperind) - 1;
		indg = upperind(changeind_g);
		
		if(Adj(indg) == 0){
			Adj_star(indg) = 1;
			Adj_star = Adj_star.t();
			Adj_star(indg) = 1;
			prop = delta/(1-delta);
		}
		else{
			Adj_star(indg) = 0;
			Adj_star = Adj_star.t();
			Adj_star(indg) = 0;
			prop = (1-delta)/delta;
		}

		if(mcard(Adj_star) == 1){
			makedecompgraph(Adj_star, Cliques_star, Separators_star);
			logmarging = log_margin(Sgamma, Cliques, Separators, bG, DG, 0, n, 0, 0.0);
			logmarging_star = log_margin(Sgamma, Cliques_star, Separators_star, bG, DG, 0, n, 0, 0.0);
			r_g = logmarging_star - logmarging + log(prop);

			if(log(unif_rand())<r_g){
				Adj = Adj_star;
				Cliques = Cliques_star;
				Separators = Separators_star;
				if(iter>burnin){
					// Adj_mcmc = Adj_mcmc + Adj;
					Adj_mcmc.slice(iter-burnin-1) = Adj;
					post(iter-burnin-1) = logmarging_star;
				}
				// cout<<"accept G"<<endl;
				// cout<<"iteration = "<<iter<<endl;
			}
			else{
				if(iter>burnin){
					// Adj_mcmc = Adj_mcmc + Adj;
					Adj_mcmc.slice(iter-burnin-1) = Adj;
					post(iter-burnin-1) = logmarging;
				}
				// cout<<"reject G"<<endl;
				// cout<<"iteration = "<<iter<<endl;
			}
		}

		else{
			if(iter>burnin){
				// Adj_mcmc = Adj_mcmc + Adj;
				Adj_mcmc.slice(iter-burnin-1) = Adj;
				post(iter-burnin-1) = logmarging;
			}
			// cout<<"reject nondecomp"<<endl;
			// cout<<"iteration = "<<iter<<endl;
		}
	}

	unsigned int m = index_max(post);
	arma::umat Adj_MAP = Adj_mcmc.slice(m);
	double MAPlike = post(m);

	PutRNGstate();

	return Rcpp::List::create(Rcpp::Named("Adj.est") = Adj_MAP,
                            Rcpp::Named("MAP.like") = MAPlike);
}




