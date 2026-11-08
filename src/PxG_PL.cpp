/**************************************************
This C++ script contains all functions for the proposed PxG model 
with pseudo-likelihood specification in section 3.2.2
It is used for simulation studies in section 4.3 & real data analysis in section 5
To reproduce all results in section 4.3 run R script "Example3.R"
To reproduce all results in section 5 run R script "RealDataAnalysis.R"
 **************************************************/

/**************************************************
There are three main functions in this C++ script: PxG_PL, CovOnly_PL, GraphOnly_PL
"PxG_PL": the proposed PxG model
"CovOnly_PL": covariate-only model where graph is assumed to be constant among samples
"GraphOnly_PL": graph-only model where graph varies with samples but does not depend on covariates
For all three functions, pseudo-likelihood model is used for likelihood and graphical prior 
 **************************************************/

/* list of arguments for all three functions */
// YData: response data or graphical data (n by q)
// XData: covariate data (n by p)
// spike: "spike" variance parameter of the spike-and-slab prior, see section 3.2.2
// slab: "slab" variance parameter of the spike-and-slab prior, see section 3.2.2
// a1: shape parameter of inverse gamma prior for the variance parameter of the spike-and-slab prior, see section 3.2.2
// a2: rate parameter of inverse gamma prior for the variance parameter of the spike-and-slab prior, see section 3.2.2
// alphaG: prior edge inclusion probability
// mu0: mean of normal priors in the auxiliary model, see section 3.1
// sigsq0: variance of priors in the auxiliary model, see section 3.1
// b1: shape parameter of inverse gamma prior in the auxiliary model, see section 3.1
// b2: rate parameter of inverse gamma prior in the auxiliary model, see section 3.1
// M: scaling parameter in the stick-breaking process
// N: maximum number of clusters allowed
// burnin: number of burn-in samples
// nmcmc: number of MCMC samples
// z0: initial cluster allocation
// B: number of graph updates in each iteration

/* list of values for all three functions */
// Label: posterior estimates of cluster labels using the least squares clustering method (Dahl, 2006) 
// DIC: DIC value
// LogTotLike: log total likelihood under the MCMC sample (the one in "Label") 
//             selected by the least squares clustering method (Dahl, 2006) 
// Var: variance of log total likelihood under all MCMC samples of cluster labels
// Gsum: 3-dimensional array (q-1 by q by n), graph selection results for each sample using pseudo-likelihood



// Package version of the authors' PxG-PL.cpp (Niu, Ni, Pati & Mallick, JASA 2024).
// Changes relative to the supplement are limited to packaging:
//   * includes reordered so that RcppArmadillo.h comes first (R.h / Rmath.h are
//     pulled in by Rcpp) and the C++11 plugin removed (current Armadillo needs at least
//     C++14, which R's default standard provides);
//   * std::cout replaced by Rcpp::Rcout, switched on by a trailing `verbose`
//     argument (default TRUE, matching the original per-iteration printing);
//   * every numerical statement is unchanged.
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <vector>
using std::vector;

/* test invgamma function */
// [[Rcpp::export]]
double testinvgamma(double b1, double b2){
  double a = 0.0;
  for(unsigned int i=0; i<1000; i++){
    a = a + 1/R::rgamma(b1, 1/b2);
  }
  return a/1000;
}



/* 
  update covariate related parameters in all N clusters
  x1, ..., xn ~ N(xMu, xSigSq * I)
          xMu ~ N(m0,  sigsq0 * xSigSq * I)
      xSigSq ~ IG(b1, b2)
  update: xSigSq | x1, ..., xn ~ IG(postb1, postb2)
             xMu | xSigSq ~ Normal(postMu, postSigSq)
*/
void UpdateCovParam(arma::uvec z, unsigned int N,
                    arma::mat X, arma::mat& xMu, arma::vec& xSigSq, 
                    arma::vec mu0, double sigsq0, double b1, double b2){ // X here is p by n
  // unsigned int n = X.n_cols;
  unsigned int p = X.n_rows; // number of covariates
  
  unsigned int ni; // number of samples in each cluster
  arma::uvec indi; // sample indices in each cluster
  arma::mat Xi;    // all samples in the ith cluster
  arma::vec postMu, munew;
  double postSigSq, postb1, postb2, sigsqnew, TSS;
  
  for(unsigned int i=0; i<N; i++){
    indi = find(z==(i+1));
    ni = indi.n_elem;
    Xi = X.cols(indi);
    
    //update covariate related parameters
    if(ni==0){
      // draw from Normal-IG prior
      xSigSq(i) = 1/R::rgamma(b1, 1/b2); // ATTENTION R::rgamma(shape, <SCALE>)
      xMu.col(i) = mu0 + sqrt(sigsq0*xSigSq(i)) * arma::randn<arma::vec>(p);
    }
    else{
      // draw variance 'xSigSq': xSigSq given all data in this cluster (not conditional)
      postb1 = b1 + 0.5*ni;
      munew = sum(Xi, 1) + mu0/sigsq0;
      sigsqnew = (ni*sigsq0+1)/sigsq0;
      TSS = arma::accu(Xi % Xi) + arma::dot(mu0, mu0)/sigsq0 - arma::dot(munew, munew)/sigsqnew;
      postb2 = b2 + 0.5*TSS;
      xSigSq(i) = 1/R::rgamma(postb1, 1/postb2); // ATTENTION R::rgamma(shape, <SCALE>)
      
      // draw mean 'xMu'
      postMu = munew/sigsqnew;
      postSigSq = xSigSq(i)/sigsqnew;
      xMu.col(i) = postMu + sqrt(postSigSq) * arma::randn<arma::vec>(p);
    }
  }
}



/* covariate log likelihood -- multi var normal */
// [[Rcpp::export]]
double CovlogLike(arma::vec x, arma::vec mu, double sigsq){
  unsigned int p = x.n_elem;
  double tss = arma::dot(x-mu, x-mu);
  
  return -0.5*p*log(2*M_PI*sigsq) - tss/(2*sigsq);
}



/* remove the one dimenion for each sample 
   data augmentation function to reshape data for neighborhood selection regression 
*/
// [[Rcpp::export]]
arma::cube shrunk(arma::mat Y){ // input Y is q by n 
  unsigned int n = Y.n_cols;
  unsigned int q = Y.n_rows;
  arma::cube Yshrunk(q-1, q, n, arma::fill::zeros);
  arma::mat Yslice(q-1, q, arma::fill::zeros);
  arma::mat y(q, 1, arma::fill::zeros);
  
  for(unsigned int i=0; i<n; i++){
    for(unsigned int j=0; j<q; j++){
      y = Y.col(i);
      y.shed_row(j);
      Yslice.col(j) = y;
    }
    Yshrunk.slice(i) = Yslice;
  }
  
  return Yshrunk;
}



/* graph log pseudo likelihood -- multi linear reg */
// [[Rcpp::export]]
double logPseudoLike(arma::vec y, arma::mat yshrunk, arma::mat Beta, arma::vec vsigsq){
  unsigned int q = y.n_elem;
  double sumreg = 0.0;
  double mu;
  
  for(unsigned int i=0; i<q; i++){
    mu = arma::dot(yshrunk.col(i), Beta.col(i));
    sumreg = sumreg - 0.5*log(2*M_PI*vsigsq(i)) - (y(i)-mu)*(y(i)-mu)/(2*vsigsq(i));
  }
  
  return sumreg;
}



/* draw a random binary vector of length q-1 according to alphaG */
// [[Rcpp::export]]
arma::uvec randg(unsigned int q, double alphaG){
  arma::uvec v(q-1, arma::fill::zeros);
  for(unsigned int i=0; i<q-1; i++){
    if(unif_rand()<alphaG) { v(i) = 1; }
  }
  return v;
}



/* draw a random vector with spike and slab variance of length q-1 */
// [[Rcpp::export]]
arma::vec varvec(arma::uvec g, unsigned int q, double spike, double slab){
  arma::vec v(q-1, arma::fill::zeros);
  for(unsigned int i=0; i<q-1; i++){
    if(g(i)==1) { v(i) = sqrt(slab)*R::rnorm(0,1); }
    else { v(i) = sqrt(spike)*R::rnorm(0,1); }
  }
  return v;
}



/* draw an inverse gamma random vector with length q */
// [[Rcpp::export]]
arma::vec gammavec(unsigned int q, double a, double b){
  arma::vec v(q, arma::fill::zeros);
  for(unsigned int i=0; i<q; i++) { v(i) = 1/R::rgamma(a, 1/b); } // ATTENTION R::rgamma(shape, <SCALE>)
  return v;
}



/* update graph related parameters in all N clusters */
void UpdateGraphParam(arma::uvec z, unsigned int N,
                      arma::mat Y, arma::cube Yshrunk,
                      vector<arma::umat>& ListG, vector<arma::mat>& ListB, arma::mat& ySigSq,
                      double spike, double slab, double a1, double a2, double alphaG, unsigned int B){ // Y here is q by n
  // unsigned int n = Y.n_cols;
  unsigned int q = Y.n_rows; // number of nodes 
  
  unsigned int ni; // number of sample in each cluster
  arma::uvec indi; // sample indices in each cluster
  arma::vec yi;    // samples in the ith cluster
  arma::cube cubeY; // neighborhood selection data format
  arma::mat subY, Yi, Ai, invAi, temp; // covariance matrix in the ith cluster
  arma::vec postMu1, postMu2;
  double posta1, posta2;
  
  for(unsigned int i=0; i<N; i++){
    indi = find(z==(i+1));
    ni = indi.n_elem;
    
    // update graph related parameters
    if(ni==0){
      // draw from prior
      for(unsigned int j=0; j<q; j++){
        ListG[i].col(j) = randg(q, alphaG);
        ListB[i].col(j) = varvec(ListG[i].col(j), q, spike, slab);
      }
      ySigSq.col(i) = gammavec(q, a1, a2);
    }
    else{
      subY = Y.cols(indi);
      subY = subY.t();
      cubeY = Yshrunk.slices(indi);
      
      for(unsigned int j=0; j<q; j++){
        yi = subY.col(j);  // ni by 1
        Yi = cubeY.col(j); // q by ni
        
        // update ListG & ListB & ySigSq
        /* updating scheme
           gamma_1 | beta_1, ySigSq
            ySigSq | beta_all
          beta_all | gamma_all, ySigSq
           gamma_2 | beta_2, ySigSq
            ySigSq | beta_all
          beta_all | gamma_all, ySigSq
          ... ...              */
        arma::uvec gamma(q-1, arma::fill::zeros);
        arma::vec   beta(q-1, arma::fill::zeros);
        double SigSq = 1.0;

        for(unsigned int round=1; round<=B; round++){

        for(unsigned int k=0; k<q-1; k++){
          /* update gamma_i given beta_i
             P(gamma_i=0 | beta_i) = P(beta_i | gamma_i=0) * P(gamma_i=0)
                                   = N(beta_i; 0, "SigSq" * spike) * (1-alphaG)
             P(gamma_i=1 | beta_i) = P(beta_i | gamma_i=1) * P(gamma_i=1)
                                   = N(beta_i; 0, "SigSq" * slab) * alphaG           
          */
          double g0, g1, prob;
          g0 = (1-alphaG) * exp(-0.5*beta(k)*beta(k)/spike/SigSq) / sqrt(spike*SigSq);
          g1 =     alphaG * exp(-0.5*beta(k)*beta(k)/slab/SigSq)  / sqrt(slab*SigSq);
          prob = g1/(g0+g1);
          if(unif_rand()<prob) { gamma(k) = 1; }
          else { gamma(k) = 0; }
          
          // given gamma construct covariance matrix for beta prior
          arma::uvec ind0 = find(gamma == 0);
          arma::uvec ind1 = find(gamma == 1);
          arma::vec vspike(ind0.n_elem, arma::fill::ones);
          arma::vec  vslab(ind1.n_elem, arma::fill::ones);
          vspike = vspike / spike; // " / " or " * "
          vslab = vslab / slab;    // " / " or " * "
          arma::vec invsigsq(q-1, arma::fill::zeros);
          invsigsq.elem(ind0) = vspike;
          invsigsq.elem(ind1) = vslab;
          arma::mat invSigSq = diagmat(invsigsq);
          
          /* update ySigSq given beta_all
             yi ~ N(Yi*beta_all, SigSq * I)
             beta_all | gamma_all, SigSq ~ N(0, "SigSq" * Sigma)
             SigSq ~ IG(a1, a2)
             SigSq | yi, Yi, beta_all ~ IG(posta1, posta2)
          */
          posta1 = a1 + 0.5*ni + 0.5*(q-1);
          postMu2 = Yi.t() * beta;
          temp = beta.t()*invSigSq*beta;
          posta2 = a2 + 0.5*arma::dot(yi-postMu2, yi-postMu2) + 0.5*temp(0,0);
          SigSq = 1/R::rgamma(posta1, 1/posta2);
          
          /* update beta_all given gamma_all and SigSq
             yi ~ N(Yi*beta_all, SigSq * I)
             beta_all | gamma_all ~ N(0, "SigSq" * Sigma)
             beta_all | yi, Yi, gamma_all, SigSq ~ N(postMu1, *)
          */
          Ai = Yi*Yi.t() + invSigSq;
          invAi = inv(Ai);
          postMu1 = invAi*Yi*yi;
          beta = postMu1 + sqrt(SigSq) * chol(invAi, "lower") * arma::randn<arma::vec>(q-1);
        }

        }
        ListG[i].col(j) = gamma;
        ySigSq(j,i) = SigSq;
        ListB[i].col(j) = beta;
      }
    }
  }
}



/* get membership matrix */
// [[Rcpp::export]]
arma::umat getMemberMat(arma::uvec z){
  unsigned int n = z.n_elem;
  arma::uvec zs = unique(z);
  arma::umat Membership(n, n, arma::fill::zeros);
  
  for(unsigned int i:zs){
    arma::umat indi = find(z==i);
    indi = sort(indi);
    unsigned int ni = indi.n_elem;
    arma::umat One1(1, ni-1, arma::fill::ones); 
    for(unsigned int j=0; j<ni; j++){
      arma::umat indT = indi;
      indT.shed_row(j);
      arma::uvec indr(1);
      indr(0) = indi(j);
      arma::uvec indc = indT;
      Membership.submat(indr, indc) = One1;
    }
  }
  
  return Membership;
}



/* discrete sampler with respect to a probability vector */
// [[Rcpp::export]]
unsigned int DiscSample(arma::vec prob){
  double u = unif_rand();
  unsigned int i = 0;
  double Sum = prob(0);
  
  while(Sum<u){
    i = i + 1;
    Sum = Sum + prob(i);
  }
  return i;
}



/* calculate the total log likelihood w.r.t. the current partition */
double getlogTotLikeX(arma::uvec z, arma::mat X, // X here is p by n
                      arma::mat xMu, arma::vec xSigSq){
  unsigned int p = X.n_rows; // number of covariates
  arma::uvec zs = unique(z);
  arma::uvec indi; // indices of samples in the ith cluster
  unsigned int ni; // number of samples in the ith cluster
  arma::mat Xi;    // samples in the ith cluster
  double TSS;
  double logTotX = 0.0;
  
  for(unsigned int i:zs){ // i from 1 to N
    indi = find(z==i);
    ni = indi.n_elem;
    Xi = X.cols(indi);
    TSS = arma::accu(Xi % Xi) - 2*arma::dot(sum(Xi, 1), xMu.col(i-1)) + ni*arma::dot(xMu.col(i-1), xMu.col(i-1));
    logTotX = logTotX - 0.5*p*ni*log(2*M_PI*xSigSq(i-1)) - 0.5*TSS/xSigSq(i-1);
  }
  
  return logTotX;
}



/* calculate the total GRAPH log likelihood w.r.t. the current partition */
double getlogTotLikeY(arma::uvec z, arma::mat Y, arma::cube Yshrunk, // Y here is q by n
                      vector<arma::mat> ListB, arma::mat ySigSq){
  unsigned int q = Y.n_rows; // number of nodes 
  arma::uvec zs = unique(z);
  arma::uvec indi; // indices of samples in the ith cluster
  unsigned int ni; // number of samples in the ith cluster
 
  arma::vec yi;    // samples in the ith cluster
  arma::cube cubeY; // neighborhood selection data format
  arma::mat subY, Yi; // covariance matrix in the ith cluster
  arma::vec Mu;
  double logTotY = 0.0;
  
  for(unsigned int i:zs){
    indi = find(z==i);
    ni = indi.n_elem;
    
    subY = Y.cols(indi);
    subY = subY.t();
    cubeY = Yshrunk.slices(indi);
    
    for(unsigned int j=0; j<q; j++){
      yi = subY.col(j);  // ni by 1
      Yi = cubeY.col(j); // q by ni
      Mu = Yi.t()*ListB[i-1].col(j);
      logTotY = logTotY - 0.5*ni*log(2*M_PI*ySigSq(j,i-1)) - 0.5*arma::dot(yi-Mu, yi-Mu)/ySigSq(j,i-1);
    }
  }
  
  return logTotY;
}



/* calculate the total log likelihood w.r.t. the current partition */
double getlogTotLike(arma::uvec z, arma::mat X, arma::mat Y, arma::cube Yshrunk,
                     arma::mat xMu, arma::vec xSigSq,
                     vector<arma::mat> ListB, arma::mat ySigSq){
  unsigned int p = X.n_rows; // number of covariates
  unsigned int q = Y.n_rows; // number of nodes 
  arma::uvec zs = unique(z);
  arma::uvec indi; // indices of samples in the ith cluster
  unsigned int ni; // number of samples in the ith cluster
  
  arma::mat Xi;    // samples in the ith cluster
  double TSS;
  arma::vec yi;    // samples in the ith cluster
  arma::cube cubeY; // neighborhood selection data format
  arma::mat subY, Yi; // covariance matrix in the ith cluster
  arma::vec Mu;

  double logTotY = 0.0;
  double logTotX = 0.0;
  for(unsigned int i:zs){
    indi = find(z==i);
    ni = indi.n_elem;
    
    Xi = X.cols(indi);
    TSS = arma::accu(Xi % Xi) - 2*arma::dot(sum(Xi, 1), xMu.col(i-1)) + ni*arma::dot(xMu.col(i-1), xMu.col(i-1));
    logTotX = logTotX - 0.5*p*ni*log(2*M_PI*xSigSq(i-1)) - 0.5*TSS/xSigSq(i-1);
 
    subY = Y.cols(indi);
    subY = subY.t();
    cubeY = Yshrunk.slices(indi);
    for(unsigned int j=0; j<q; j++){
      yi = subY.col(j);  // ni by 1
      Yi = cubeY.col(j); // q by ni
      Mu = Yi.t()*ListB[i-1].col(j);
      logTotY = logTotY - 0.5*ni*log(2*M_PI*ySigSq(j,i-1)) - 0.5*arma::dot(yi-Mu, yi-Mu)/ySigSq(j,i-1);
    }
  }
  
  return logTotY + logTotX;
}



/* run MCMC chain to get log total likelihood for estimating mean and variance for covariates only as one cluster */
// [[Rcpp::export]]
arma::vec getlogTotLikeCov(arma::mat XData, 
                        arma::vec mu0, double sigsq0, double b1, double b2,
                        unsigned int burnin, unsigned int nmcmc){
  unsigned int n = XData.n_rows;
  unsigned int p = XData.n_cols;
  arma::mat X = XData.t();
  arma::vec logTotLike(nmcmc, arma::fill::zeros);
  
  arma::vec postMu, munew, PostMean;
  double postSigSq, postb1, postb2, sigsqnew, TSS, PostVar, RSS;
  double SV = 0.0; arma::vec SM(p, arma::fill::zeros);
  
  for(unsigned int iter=1; iter<=(burnin+nmcmc); iter++){
    // draw variance given all data in this cluster (not conditional)
    postb1 = b1 + 0.5*n;
    munew = sum(X, 1) + mu0/sigsq0;
    sigsqnew = (n*sigsq0+1)/sigsq0;
    TSS = arma::accu(X % X) + arma::dot(mu0, mu0)/sigsq0 - arma::dot(munew, munew)/sigsqnew;
    postb2 = b2 + 0.5*TSS;
    PostVar = 1/R::rgamma(postb1, 1/postb2); // ATTENTION R::rgamma(shape, <SCALE>)
    
    // draw mean given variance
    postMu = munew/sigsqnew;
    postSigSq = PostVar/sigsqnew;
    PostMean = postMu + sqrt(postSigSq) * arma::randn<arma::vec>(p);
    
    if(iter>burnin){
      SV = SV + PostVar; SM = SM + PostMean;
      RSS = arma::accu(X % X) - 2*arma::dot(sum(X, 1), PostMean) + n*arma::dot(PostMean, PostMean);
      logTotLike(iter-burnin-1) = -0.5*p*n*log(2*M_PI*PostVar) - 0.5*TSS/PostVar;
    }
  }
  // cout<<"posterior mean of mu: "<<SM/nmcmc<<endl;
  // cout<<"posterior mean of sigsq: "<<SV/nmcmc<<endl;
  
  return logTotLike;
}



/* run MCMC chain using neighbor selection to get log total likelihood for Y as one cluster */
arma::vec getlogTotLikeG(arma::mat YData,
                         double spike, double slab, double a1, double a2, double alphaG, 
                         unsigned int burnin, unsigned int nmcmc, bool verbose){
  unsigned int n = YData.n_rows;
  unsigned int q = YData.n_cols;
  arma::mat Y = YData.t();
  arma::cube Yshrunk = shrunk(Y);
  Y = Y.t();
  arma::vec logTotLike(nmcmc, arma::fill::zeros);
  
  arma::vec yi, Mu;
  arma::mat Yi, Ai, invAi, temp;
  arma::vec postMu1, postMu2;
  double posta1, posta2;
  
  for(unsigned int iter=1; iter<=(burnin+nmcmc); iter++){
    if(verbose) Rcpp::Rcout<<"iteration i = "<<iter<<std::endl;
    
    for(unsigned int j=0; j<q; j++){
      yi = Y.col(j);
      Yi = Yshrunk.col(j);
      
      arma::uvec gamma(q-1, arma::fill::zeros);
      arma::vec   beta(q-1, arma::fill::zeros);
      double SigSq = 1.0;
      for(unsigned int k=0; k<q-1; k++){
        double g0, g1, prob;
        g0 = (1-alphaG) * exp(-0.5*beta(k)*beta(k)/spike/SigSq) / sqrt(spike*SigSq);
        g1 =     alphaG * exp(-0.5*beta(k)*beta(k)/slab/SigSq)  / sqrt(slab*SigSq);
        prob = g1/(g0+g1);
        if(unif_rand()<prob) { gamma(k)=1; }
        
        // given gamma construct covariance matrix for beta prior
        arma::uvec ind0 = find(gamma == 0);
        arma::uvec ind1 = find(gamma == 1);
        arma::vec vspike(ind0.n_elem, arma::fill::ones);
        arma::vec  vslab(ind1.n_elem, arma::fill::ones);
        vspike = vspike / spike; // " / " or " * "
        vslab = vslab / slab;    // " / " or " * "
        arma::vec invsigsq(q-1, arma::fill::zeros);
        invsigsq.elem(ind0) = vspike;
        invsigsq.elem(ind1) = vslab;
        arma::mat invSigSq = diagmat(invsigsq);

        posta1 = a1 + 0.5*n + 0.5*(q-1);
        postMu2 = Yi.t() * beta;
        temp = beta.t()*invSigSq*beta;
        posta2 = a2 + 0.5*arma::dot(yi-postMu2, yi-postMu2) + 0.5*temp(0,0);
        SigSq = 1/R::rgamma(posta1, 1/posta2);
        
        Ai = Yi*Yi.t() + invSigSq;
        invAi = inv(Ai);
        postMu1 = invAi*Yi*yi;
        beta = postMu1 + sqrt(SigSq) * chol(invAi, "lower") * arma::randn<arma::vec>(q-1);
      }
      
      if(iter>burnin){
        Mu = Yi.t()*beta;
        logTotLike(iter-burnin-1) = logTotLike(iter-burnin-1) - 0.5*n*log(2*M_PI*SigSq) - 0.5*arma::dot(yi-Mu, yi-Mu)/SigSq;
      }
    }
  }

  return logTotLike;
}



/* blocked Gibbs sampler for the covariate-only model */
// [[Rcpp::export]]
Rcpp::List CovOnly_PL(arma::mat YData, arma::mat XData, 
                      double spike, double slab, double a1, double a2, double alphaG, 
                      arma::vec mu0, double sigsq0, double b1, double b2, 
                      double M, unsigned int N, // N is the maximum number of clusters
                      unsigned int burnin, unsigned int nmcmc, arma::uvec z0, bool verbose = true){ // z0 is the initial labeling
  unsigned int n = XData.n_rows;
  unsigned int p = XData.n_cols;
  arma::mat X = XData.t(); // using transpose for easy vector allocation
  
  arma::umat Z(nmcmc, n, arma::fill::zeros); // recording all MCMC labeling samples
  arma::ucube Membership(n, n, nmcmc, arma::fill::zeros); // membership matrix
  arma::mat SUM(n, n, arma::fill::zeros);
  
  arma::uvec z = z0; // initial cluster allocation
  
  arma::vec V(N, arma::fill::zeros); // stick-breaking weights
  V(N-1) = 1.0;
  arma::vec prob(N, arma::fill::zeros); // cluster probability
  arma::uvec nc(N, arma::fill::zeros);  // number of samples in each cluster
  
  arma::mat xMu(p, N, arma::fill::zeros); // cluster mean
  arma::vec xSigSq(N, arma::fill::ones);  // cluster variance
  
  arma::uvec zi;
  arma::vec weights(N, arma::fill::zeros); // posterior prob for each cluster
  double maxweight, mindist;
  
  arma::vec logTotLikeX(nmcmc, arma::fill::zeros); // log total likelihood for covariates
  arma::vec logTotLikeY(nmcmc, arma::fill::zeros); // log total likelihood for graphs 
  
  for(unsigned int iter=1; iter<=(burnin+nmcmc); iter++){
    if(verbose) Rcpp::Rcout<<"iteration i = "<<iter<<std::endl;
    
    // first counting number of samples in each current cluster
    for(unsigned int i=1; i<=N; i++) {
      zi = find(z==i);
      nc(i-1) = zi.n_elem;
    }
    
    for(unsigned int i=0; i<(N-1); i++){
      arma::uvec ind = arma::linspace<arma::uvec>(i+1, N-1, N-i-1);
      V(i) = R::rbeta(1+nc(i), M+sum(nc(ind)));
    }
    
    arma::vec Vp = 1 - V;
    prob(0) = V(0);
    for(unsigned int i=1; i<N; i++){
      arma::uvec ind = arma::linspace<arma::uvec>(0, i-1, i);
      prob(i) = V(i)*prod(Vp(ind));
    }
    
    for(unsigned int i=0; i<n; i++){
      for(unsigned int j=0; j<N; j++){
        weights(j) = log(prob(j)) + CovlogLike(X.col(i), xMu.col(j), xSigSq(j));
      }
      
      maxweight = max(weights);
      weights = weights - maxweight;
      weights = exp(weights);
      weights = weights/sum(weights);
      z(i) = DiscSample(weights) + 1;
    }
    // update all cluster-related parameters
    UpdateCovParam(z, N, X, xMu, xSigSq, mu0, sigsq0, b1, b2);
    arma::uvec aa = find(z==1);
    // cout<<aa.n_elem<<endl;
    // cout<<xMu<<endl; cout<<xSigSq<<endl;
    
    if(iter>burnin){
      Z.row(iter-burnin-1) = z.t();
      Membership.slice(iter-burnin-1) = getMemberMat(z);
      SUM = SUM + Membership.slice(iter-burnin-1);
      
      /* calculate the total log likelihood w.r.t. the current partition */
      logTotLikeX(iter-burnin-1) = getlogTotLikeX(z, X, xMu, xSigSq);
    }
  }
  
  if(verbose) Rcpp::Rcout<<"Calculating Membership Matrix..."<<std::endl;
  SUM = SUM/nmcmc;
  arma::vec dist(nmcmc, arma::fill::zeros);
  for(unsigned int i=0; i<nmcmc; i++){
    dist(i) = arma::accu( (Membership.slice(i) - SUM) % (Membership.slice(i) - SUM) );
  }
  mindist = min(dist);
  arma::uvec minid = find(dist==mindist);
  arma::uvec result = (Z.row(minid(0))).t();

  if(verbose) Rcpp::Rcout<<"Calculating DIC..."<<std::endl;
  logTotLikeY = getlogTotLikeG(YData, spike, slab, a1, a2, alphaG, burnin, nmcmc, verbose);
  double maxG = max(logTotLikeY);
  arma::uvec maxGid = find(logTotLikeY==maxG);
  
  if(verbose) Rcpp::Rcout<<"DIC = "<< -2*logTotLikeX(minid(0)) + var(logTotLikeX) - 2*logTotLikeY(maxGid(0)) + var(logTotLikeY) <<std::endl;
  if(verbose) Rcpp::Rcout<<"Log Total Likelihood is "<< logTotLikeX(minid(0)) + logTotLikeY(maxGid(0)) <<std::endl;
  if(verbose) Rcpp::Rcout<<"Variance of Log Tot Like is "<< var(logTotLikeX) + var(logTotLikeY) <<std::endl;

  return Rcpp::List::create(Rcpp::Named("Label") = result,
                            Rcpp::Named("DIC") = -2*logTotLikeX(minid(0)) + var(logTotLikeX) - 2*logTotLikeY(maxGid(0)) + var(logTotLikeY),
                            Rcpp::Named("logTotLike") = logTotLikeX(minid(0)) + logTotLikeY(maxGid(0)),
                            Rcpp::Named("Var") = var(logTotLikeX) + var(logTotLikeY));
}



/* blocked Gibbs sampler for the graph-only model */
// [[Rcpp::export]]
Rcpp::List GraphOnly_PL(arma::mat YData, arma::mat XData, 
                        double spike, double slab, double a1, double a2, double alphaG, 
                        arma::vec mu0, double sigsq0, double b1, double b2, 
                        double M, unsigned int N, // N is the maximum number of clusters
                        unsigned int burnin, unsigned int nmcmc, arma::uvec z0, unsigned int B, bool verbose = true){ // z0 is the initial labeling
  unsigned int n = YData.n_rows;
  unsigned int q = YData.n_cols;
  arma::mat Y = YData.t(); // using transpose for easy vector allocation
  arma::cube Yshrunk = shrunk(Y);
  
  arma::umat Z(nmcmc, n, arma::fill::zeros); // recording all MCMC labeling samples
  arma::ucube Membership(n, n, nmcmc, arma::fill::zeros); // membership matrix
  arma::mat SUM(n, n, arma::fill::zeros);
  
  arma::uvec z = z0; // initial cluster allocation
  
  arma::vec V(N, arma::fill::zeros); // stick-breaking weights
  V(N-1) = 1.0;
  arma::vec prob(N, arma::fill::zeros); // cluster probability
  arma::uvec nc(N, arma::fill::zeros);  // number of samples in each cluster
  
  vector<arma::umat> ListG;
  vector<arma::mat> ListB;
  arma::umat G0(q-1, q, arma::fill::zeros);
  arma::mat  B0(q-1, q, arma::fill::zeros);
  for(unsigned int i=0; i<N; i++) { ListG.push_back(G0); ListB.push_back(B0); }
  arma::mat ySigSq(q, N, arma::fill::ones);
  
  arma::uvec zi;
  arma::vec weights(N, arma::fill::zeros); // posterior prob for each cluster
  double maxweight, mindist;
  
  arma::vec logTotLikeY(nmcmc, arma::fill::zeros); // log total likelihood for graphs
  arma::vec logTotLikeX(nmcmc, arma::fill::zeros); // log total likelihood for covariates
  
  for(unsigned int iter=1; iter<=(burnin+nmcmc); iter++){
    if(verbose) Rcpp::Rcout<<"iteration i = "<<iter<<std::endl;
    
    // first counting number of samples in each current cluster
    for(unsigned int i=1; i<=N; i++) {
      zi = find(z==i);
      nc(i-1) = zi.n_elem;
    }
    
    for(unsigned int i=0; i<(N-1); i++){
      arma::uvec ind = arma::linspace<arma::uvec>(i+1, N-1, N-i-1);
      V(i) = R::rbeta(1+nc(i), M+sum(nc(ind)));
    }
    
    arma::vec Vp = 1 - V;
    prob(0) = V(0);
    for(unsigned int i=1; i<N; i++){
      arma::uvec ind = arma::linspace<arma::uvec>(0, i-1, i);
      prob(i) = V(i)*prod(Vp(ind));
    }
    
    for(unsigned int i=0; i<n; i++){
      for(unsigned int j=0; j<N; j++){
        weights(j) = log(prob(j)) + logPseudoLike(Y.col(i), Yshrunk.slice(i), ListB[j], ySigSq.col(j));
      }
      
      maxweight = max(weights);
      weights = weights - maxweight;
      weights = exp(weights);
      weights = weights/sum(weights);
      z(i) = DiscSample(weights) + 1;
    }
    // update all cluster-related parameters
    UpdateGraphParam(z, N, Y, Yshrunk, ListG, ListB, ySigSq, spike, slab, a1, a2, alphaG, B);
    
    if(iter>burnin){
      Z.row(iter-burnin-1) = z.t();
      Membership.slice(iter-burnin-1) = getMemberMat(z);
      SUM = SUM + Membership.slice(iter-burnin-1);
      
      /* calculate the total log likelihood w.r.t. the current partition */
      logTotLikeY(iter-burnin-1) = getlogTotLikeY(z, Y, Yshrunk, ListB, ySigSq);
    }
  }
  
  if(verbose) Rcpp::Rcout<<"Calculating Membership Matrix..."<<std::endl;
  SUM = SUM/nmcmc;
  arma::vec dist(nmcmc, arma::fill::zeros);
  for(unsigned int i=0; i<nmcmc; i++){
    dist(i) = arma::accu( (Membership.slice(i) - SUM) % (Membership.slice(i) - SUM) );
  }
  mindist = min(dist);
  arma::uvec minid = find(dist==mindist);
  arma::uvec result = (Z.row(minid(0))).t();
  
  if(verbose) Rcpp::Rcout<<"Calculating DIC..."<<std::endl;
  /*
  unsigned int p = XData.n_cols;
  arma::mat X = XData.t(); // X is p by n
  arma::vec MuHat = arma::mean(X, 1);
  double RSS = arma::accu(X % X) - 2*arma::dot(sum(X, 1), MuHat) + n*arma::dot(MuHat, MuHat);
  logTotLikeX = -0.5*p*n*log(2*M_PI*RSS/(p*n)) - 0.5*p*n;
  */
  logTotLikeX = getlogTotLikeCov(XData, mu0, sigsq0, b1, b2, burnin, nmcmc);
  double maxX = max(logTotLikeX);
  arma::uvec maxXid = find(logTotLikeX==maxX);
  
  if(verbose) Rcpp::Rcout<<"DIC = "<< -2*logTotLikeY(minid(0)) + var(logTotLikeY) - 2*logTotLikeX(maxXid(0)) + var(logTotLikeX) <<std::endl;
  if(verbose) Rcpp::Rcout<<"Log Total Likelihood is "<< logTotLikeY(minid(0)) + logTotLikeX(maxXid(0)) <<std::endl;
  if(verbose) Rcpp::Rcout<<"Variance of Log Tot Like is "<< var(logTotLikeY) + var(logTotLikeX) <<std::endl;
  
  return Rcpp::List::create(Rcpp::Named("Label") = result,
                            Rcpp::Named("DIC") = -2*logTotLikeY(minid(0)) + var(logTotLikeY) - 2*logTotLikeX(maxXid(0)) + var(logTotLikeX),
                            Rcpp::Named("logTotLike") = logTotLikeY(minid(0)) + logTotLikeX(maxXid(0)),
                            Rcpp::Named("Var") = var(logTotLikeY) + var(logTotLikeX));
}



/* blocked Gibbs sampler for the PxG model */
// [[Rcpp::export]]
Rcpp::List PxG_PL(arma::mat YData, arma::mat XData, 
                  double spike, double slab, double a1, double a2, double alphaG, 
                  arma::vec mu0, double sigsq0, double b1, double b2, 
                  double M, unsigned int N, // N is the maximum number of clusters
                  unsigned int burnin, unsigned int nmcmc, arma::uvec z0, unsigned int B, bool verbose = true){ // z0 is the initial labeling
  unsigned int n = YData.n_rows;
  unsigned int q = YData.n_cols;
  unsigned int p = XData.n_cols;
  arma::mat X = XData.t(); // using transpose for easy vector allocation
  arma::mat Y = YData.t(); // using transpose for easy vector allocation
  arma::cube Yshrunk = shrunk(Y);
  
  arma::umat Z(nmcmc, n, arma::fill::zeros); // recording all MCMC labeling samples
  arma::ucube Membership(n, n, nmcmc, arma::fill::zeros); // membership matrix
  arma::mat SUM(n, n, arma::fill::zeros);
  
  arma::uvec z = z0; // initial cluster allocation
  
  arma::vec V(N, arma::fill::zeros); // stick-breaking weights
  V(N-1) = 1.0;
  arma::vec prob(N, arma::fill::zeros); // cluster probability
  arma::uvec nc(N, arma::fill::zeros);  // number of samples in each cluster
  
  vector<arma::umat> ListG;
  vector<arma::mat> ListB;
  arma::umat G0(q-1, q, arma::fill::zeros);
  arma::mat  B0(q-1, q, arma::fill::zeros);
  for(unsigned int i=0; i<N; i++) { ListG.push_back(G0); ListB.push_back(B0); }
  arma::mat ySigSq(q, N, arma::fill::ones);
  
  arma::mat xMu(p, N, arma::fill::zeros); // cluster covariates mean
  arma::vec xSigSq(N, arma::fill::ones);  // cluster covariates variance
  
  arma::uvec zi;
  arma::vec weights(N, arma::fill::zeros); // posterior prob for each cluster
  double maxweight, mindist;
  
  arma::vec logTotLike(nmcmc, arma::fill::zeros); // record log total likelihood
  arma::ucube Gx(q-1, q, n, arma::fill::zeros);
  
  for(unsigned int iter=1; iter<=(burnin+nmcmc); iter++){
    if(verbose) Rcpp::Rcout<<"iteration i = "<<iter<<std::endl;
    
    // first counting number of samples in each current cluster
    for(unsigned int i=1; i<=N; i++) {
      zi = find(z==i);
      nc(i-1) = zi.n_elem;
    }
    
    for(unsigned int i=0; i<(N-1); i++){
      arma::uvec ind = arma::linspace<arma::uvec>(i+1, N-1, N-i-1);
      V(i) = R::rbeta(1+nc(i), M+sum(nc(ind)));
    }
    
    arma::vec Vp = 1 - V;
    prob(0) = V(0);
    for(unsigned int i=1; i<N; i++){
      arma::uvec ind = arma::linspace<arma::uvec>(0, i-1, i);
      prob(i) = V(i)*prod(Vp(ind));
    }
    
    for(unsigned int i=0; i<n; i++){
      for(unsigned int j=0; j<N; j++){
        weights(j) = log(prob(j)) + logPseudoLike(Y.col(i), Yshrunk.slice(i), ListB[j], ySigSq.col(j));
        weights(j) = weights(j) + CovlogLike(X.col(i), xMu.col(j), xSigSq(j));
      }
      
      maxweight = max(weights);
      weights = weights - maxweight;
      weights = exp(weights);
      weights = weights/sum(weights);
      z(i) = DiscSample(weights) + 1;
    }
    // update all cluster-related parameters of Y
    UpdateGraphParam(z, N, Y, Yshrunk, ListG, ListB, ySigSq, spike, slab, a1, a2, alphaG, B);
    // update all cluster-related parameters of X
    UpdateCovParam(z, N, X, xMu, xSigSq, mu0, sigsq0, b1, b2);
    // cout<<ySigSq<<endl;
    
    if(iter>burnin){
      Z.row(iter-burnin-1) = z.t();
      Membership.slice(iter-burnin-1) = getMemberMat(z);
      SUM = SUM + Membership.slice(iter-burnin-1);

      /* calculate the total log likelihood w.r.t. the current partition */
      logTotLike(iter-burnin-1) = getlogTotLike(z, X, Y, Yshrunk, xMu, xSigSq, ListB, ySigSq);

      for(unsigned int s=0; s<n; s++){ Gx.slice(s) = Gx.slice(s) + ListG[z(s)-1]; }
    }
  }
  
  if(verbose) Rcpp::Rcout<<"Calculating Membership Matrix..."<<std::endl;
  SUM = SUM/nmcmc;
  arma::vec dist(nmcmc, arma::fill::zeros);
  for(unsigned int i=0; i<nmcmc; i++){
    dist(i) = arma::accu( (Membership.slice(i) - SUM) % (Membership.slice(i) - SUM) );
  }
  mindist = min(dist);
  arma::uvec minid = find(dist==mindist);
  arma::uvec result = (Z.row(minid(0))).t();
  
  if(verbose) Rcpp::Rcout<<"Calculating DIC..."<<std::endl;
  if(verbose) Rcpp::Rcout<<"DIC = "<< -2*logTotLike(minid(0)) + var(logTotLike) <<std::endl;
  if(verbose) Rcpp::Rcout<<"Log Total Likelihood is "<< logTotLike(minid(0)) <<std::endl;
  if(verbose) Rcpp::Rcout<<"Variance of Log Tot Like is "<< var(logTotLike) <<std::endl;
  
  return Rcpp::List::create(Rcpp::Named("Label") = result,
                            Rcpp::Named("DIC") = -2*logTotLike(minid(0)) + var(logTotLike),
                            Rcpp::Named("logTotLike") = logTotLike(minid(0)),
                            Rcpp::Named("Var") = var(logTotLike),
                            Rcpp::Named("Gsum") = Gx);
}



/* estimating graph via neighbor selection / pseudo likelihood */
// [[Rcpp::export]]
arma::mat Graph_NS(arma::mat YData, double spike, double slab, 
                   double a1, double a2, double alphaG,
                   unsigned int burnin, unsigned int nmcmc, bool verbose = true){
  
  unsigned int n = YData.n_rows;
  unsigned int q = YData.n_cols;
  arma::mat Y = YData.t(); // using transpose for easy vector allocation
  arma::cube Yshrunk = shrunk(Y);
  
  double posta1, posta2;
  arma::vec postMu1, postMu2;
  arma::mat yi, Yi, Ai, invAi, temp;
  
  arma::mat Est(q, q, arma::fill::zeros);
  arma::mat Gamma(q-1, q, arma::fill::zeros);

  for(unsigned int i=0; i<q; i++){
    yi = YData.col(i);
    Yi = Yshrunk.col(i);
    
    arma::uvec gamma(q-1, arma::fill::zeros);
    arma::vec  beta(q-1, arma::fill::ones);
    double SigSq = 1.0;
    
    for(unsigned int iter=1; iter<(burnin+nmcmc); iter++){
      if(verbose) Rcpp::Rcout<<"iteration = "<<iter<<"  node = "<<i+1<<std::endl;
      
      for(unsigned int j=0; j<q-1; j++){ // update q-1 indicators
        double g0, g1, prob;
        g0 = (1-alphaG) * exp(-0.5*beta(j)*beta(j)/spike/SigSq) / sqrt(spike*SigSq);
        g1 =     alphaG * exp(-0.5*beta(j)*beta(j)/slab/SigSq)  / sqrt(slab*SigSq);
        prob = g1/(g0+g1);
        if(unif_rand()<prob) { gamma(j) = 1; }
        else { gamma(j) = 0; }
        
        arma::uvec ind0 = find(gamma == 0);
        arma::uvec ind1 = find(gamma == 1);
        arma::vec vspike(ind0.n_elem, arma::fill::ones);
        arma::vec  vslab(ind1.n_elem, arma::fill::ones);
        vspike = vspike / spike; // " / " or " * "
        vslab = vslab / slab;    // " / " or " * "
        arma::vec invsigsq(q-1, arma::fill::zeros);
        invsigsq.elem(ind0) = vspike;
        invsigsq.elem(ind1) = vslab;
        arma::mat invSigSq = diagmat(invsigsq);
        
        posta1 = a1 + 0.5*n + 0.5*(q-1);
        postMu2 = Yi.t() * beta;
        temp = beta.t()*invSigSq*beta;
        posta2 = a2 + 0.5*arma::dot(yi-postMu2, yi-postMu2) + 0.5*temp(0,0);
        SigSq = 1/R::rgamma(posta1, 1/posta2);
        
        Ai = Yi*Yi.t() + invSigSq;
        invAi = inv(Ai);
        postMu1 = invAi*Yi*yi;
        beta = postMu1 + sqrt(SigSq) * chol(invAi, "lower") * arma::randn<arma::vec>(q-1);
      }
      if(iter > burnin){ Gamma.col(i) = Gamma.col(i) + gamma; }
      
    }
  }
  
  for(unsigned int i=0; i<q; i++){
    for(unsigned int j=0; j<q; j++){
      if(i<j) { Est(i,j) = Gamma(j-1,i)/nmcmc; }
      if(i>j) { Est(i,j) = Gamma(j,i)/nmcmc; }
    }
  }
  
  return Est;
}


