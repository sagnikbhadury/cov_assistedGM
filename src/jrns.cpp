// JRNS competitor (Samanta, Khare & Michailidis 2022): the three files JRNS.cpp, StepB.cpp and
// StepOmega.cpp distributed in JRNS.zip of the JASA supplement of Niu et al. (2024), concatenated
// into one translation unit because RcppArmadilloExtensions/sample.h defines non-inline functions
// and may be included only once per shared library. Exported signatures are namespace-qualified
// (arma::mat, Rcpp::List) so that Rcpp's generated wrappers compile. Packaging changes only, see
// inst/README_sources.md.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <ctime>

using namespace Rcpp;
using namespace arma;
using namespace std;

// ---------------------------------------------------------------- JRNS.cpp
// [[Rcpp::export]]
Rcpp::List JRNS(arma::mat B_ini, arma::mat Omega_ini, arma::mat Y, arma::mat X, int nmc, int burnin, double r1 = 1e-4, double s1 = 1e-8,double alpha1 = 1,double beta1 = 1,double alpha2 = 1, double beta2 = 1,int posdef = 0, int hyp = 0){
	
	
	int n = Y.n_rows;
	int p = B_ini.n_rows;
	int q = B_ini.n_cols;
	mat XtX = X.t()*X;
	mat XtY = X.t()*Y;
	mat B = zeros(p,q);
	B = B_ini ;
	mat M2 = B.t()*X.t()*X;
	mat Omega = zeros(q,q);
	Omega = Omega_ini;
	mat Omega2 = Omega*Omega;
	mat Bhat = zeros(p,q);
	mat Omegahat = zeros(q,q);
	
	double q1 = 1-(1/(double)p);
	vec resid =zeros(p);
	
	double b ;
	double c1;
	double c2 ;
	double a = 0;
    double d2 = 0;
    double lambda = 0;
	double eta = 0;
    double q2 = 1-(1/(double)q);
	double new1;
	double v,rho;
	
	int components = 1;
	int rand = 1;
	double bd = 0;
	double nb, nomega;
	vec P_B = zeros(2);
	vec var_mod = zeros(q);
	cube Bcount = zeros(p,q,2);
	mat Bsum = zeros(p,q);
	cube count2 = zeros(q,q,2);
	mat sum2 = zeros(q,q);
	cube out1 = zeros(p,q,nmc+1);
	cube out2 = zeros(q,q,nmc+1);
	cube B_avg = zeros(p,q,nmc+burnin+1);
    vec myvec = zeros(2);
	vec mhvec = zeros(2);
    vec P = zeros(2);
	vec probdiag = zeros(2);
	mat I = eye(q,q);
	
    clock_t time_start, time_end;
    time_start = clock();
    
    for(int t=0; t<2;t++){
		myvec(t) = t;
	}
    
	for(int it = 0; it < (nmc + burnin+1); it++){
		
		
		mat M1 = XtY*Omega2;
		
		if(hyp == 1){
		uvec idx1 = find(B != 0);
	    nb = idx1.size();	
		q1 = R::rbeta((double)p * (double)q - nb + alpha1, nb + beta1);	
		}			
		
			
		for(int r=0; r<p; r++){
		
		  for(int s=0; s<q; s++){	
		  
		    b = XtX(r,r)*Omega2(s,s);
		  	c2 = M1(r,s) - sum(Omega2.col(s) % M2.col(r)) + B(r,s)*XtX(r,r)*Omega2(s,s);
				
			
		  	if((B(r,s) == 0)){
				eta = R::rgamma(r1, s1);  
			   }
			   else{
				eta = R::rgamma(r1 + 0.5, 1/(0.5*B(r,s)*B(r,s) + s1));   
				
			}
			
            if(eta == 0){
                P_B(1) =0;
			}
			else{
            P_B(1) = sqrt(eta/(b + eta)) * exp(c2*c2/(2*(b + eta))) *(1-q1)/q1;
			}			
            P_B(0) = 1;
			
		    c1 = b + eta;
			if(P_B.has_inf()==true){
			   	
                int components = 1;
                new1 = R::rnorm(c2/c1, sqrt(1/c1));
                  
                if((it > burnin)){
                  	
				Bcount(r,s,components) += 1;
                Bsum(r,s) +=  new1;
                
                }
				
            } else {
            	
            	
            	P_B = P_B/sum(P_B);
            	components = (Rcpp::RcppArmadillo::sample(myvec,1,false, P_B) )(0);
            	if(components == 1){
			        new1 = R::rnorm(c2/c1, sqrt(1/c1));
    		        if((it > burnin))
    		         Bsum(r,s) +=  new1;
				 
				} else {
          		 
				 new1 = 0;
          		}
          		    
            	if((it > burnin))
            	Bcount(r,s,components) += 1;
            	
            	
                }
		  	M2.row(s) = M2.row(s) + (new1-B(r,s)) * XtX.row(r);
			B(r,s) = new1;
			if(it > 0){
				
				B_avg(r,s,it) =  (it*B_avg(r,s,it-1) + new1)/(it+1);
			}else{
				B_avg(r,s,it) = new1;
				
			}
		    }
    	}
		
    	if(it >= burnin)
		  out1.slice(it - burnin) = B;
	    mat E = Y-X*B;
        mat S =  E.t()*E;
		
		if(hyp == 1){
		uvec idx2 = find(Omega != 0);
	    nomega = 0.5*(idx2.size() - (double)q);
		q2 = R::rbeta(0.5 * (double)q * ((double)q -1) - nomega + alpha2, nomega + beta2);
		}
		
		
		for(int j=0; j<(q-1); j++){
            for(int k=(j+1); k<q; k++){
            	
               a = S(j,j) + S(k,k);
               d2 = - sum(Omega.col(j)%S.col(k) + Omega.col(k)%S.col(j)) + Omega(j,k)*a;
               
			   
			   
			   if((Omega(j,k) == 0)){
				lambda = R::rgamma(r1, s1);   
				
			   }
			   else{
				lambda = R::rgamma(r1 + 0.5, 1/(0.5*Omega(j,k)*Omega(j,k) + s1));   
				
			   }
               
			   if(lambda == 0) {
				   P(1) = 0;
			   }
			   else {
				   P(1) = sqrt(lambda/(a+lambda)) * exp(d2*d2/(2*(a+lambda))) *(1-q2)/q2; 
			   }
               P(0) = 1;
			   
		
               if(P.has_inf()==true){
                 int rand = 1;
                 Omega(j,k) = Omega(k,j) = R::rnorm(d2/(a+lambda), 1/sqrt(a + lambda));
				 
				 
                 if((it > burnin)){
				 
                     count2(j,k,rand) += 1;
                     sum2(j,k) +=  Omega(j,k);
					 
					 
                  }
                  
                  
                }  else {
          
          
                 P = P/sum(P);
          
                 rand = (Rcpp::RcppArmadillo::sample(myvec,1,false, P))(0);
          
          
                 if(rand == 1){
                 	
                     Omega(j,k) = Omega(k,j) = R::rnorm(d2/(a+lambda), 1/sqrt(a + lambda));
					 
					  
					 if(it > burnin){
					 sum2(j,k) +=  Omega(j,k);
					 }
                 } else {
                 	
                     Omega(j,k) = Omega(k,j) = 0;
					 
				    }
     
     
                  if((it > burnin))
                     count2(j,k,rand) += 1; 
            
                }
                
        
           }
		   
		   
           bd = sum(Omega.col(j)%S.col(j)) - Omega(j,j)*S(j,j) + (R::rgamma(r1 + 1, 1/(abs(Omega(j,j)) + s1)));
           
		   double mode = (sqrt(bd*bd + 4*S(j,j)*(double)n) - bd)/(2*S(j,j));
		   v = R::rnorm(mode, 0.001);
		   rho = std::min(exp(n * (log(v) - log(Omega(j,j))) - 0.5 * S(j,j) * ((v*v) - (Omega(j,j)*Omega(j,j))) - bd * (v - Omega(j,j))),1.0);
		   
		   
		   mhvec(0) = v;
		   mhvec(1) = Omega(j,j);
		   probdiag(0) = rho;
		   probdiag(1) = 1-rho;
		   
		   Omega(j,j) = (Rcpp::RcppArmadillo::sample(mhvec,1,false,probdiag))(0);
		   if(it > burnin)
		    sum2(j,j) += Omega(j,j);
		   
		   
		   
        }
        bd = sum(Omega.col(q-1)%S.col(q-1)) - Omega(q-1,q-1)*S(q-1,q-1) + (R::rgamma(r1 + 1, 1/(abs(Omega(q-1, q-1)) + s1)));
        double mode = (sqrt(bd*bd + 4*S(q-1,q-1)*(double)n) - bd)/(2*S(q-1,q-1));
		v = R::rnorm(mode, 0.001 );
		rho = std::min(exp(n * (log(v) - log(Omega(q-1,q-1))) - 0.5 * S(q-1,q-1) * ((v*v) - (Omega(q-1,q-1)*Omega(q-1,q-1))) - bd * (v - Omega(q-1,q-1))),1.0);
		
		
		mhvec(0) = v;
		mhvec(1) = Omega(q-1,q-1);
		probdiag(0) = rho;
		probdiag(1) = 1-rho;
		
		Omega(q-1,q-1) = (Rcpp::RcppArmadillo::sample(mhvec,1,false,probdiag))(0);
		if(it > burnin)
          sum2(q-1,q-1) += Omega(q-1,q-1);
		
		Omega2 = Omega*Omega;
		
		if( posdef == 1){
            double eigmin = min(eig_sym(Omega));
			
	        if(eigmin < 0){
	            
                Rcout << "Omega is not pd in iteration : " << it <<"\n";	   
	            mat Omeganew = Omega - (eigmin - 0.01) * I;
    	
		        if(it > burnin){
		            for(int j=0; j<q; j++){
		                sum2(j,j) = sum2(j,j) - Omega(j,j) + Omeganew(j,j);
	                }
		        }
	            Omega = Omeganew;
	        }
	    }
		
	if(it - floor(it/1000)*1000 == 0){
      vec iteration = ones(1)*it;
      iteration.print();
    }	
    if(it >=burnin)
    out2.slice(it - burnin) = Omega; 	
            	         	
    }
	
	for (int mi = 0 ; mi < p ; mi++){
       for(int mj = 0 ; mj < q ; mj++){
       	
          if(Bcount(mi,mj,0)> (0.5*nmc)){
        
             Bhat(mi,mj) = 0;
        
           }else{
              Bhat(mi,mj) = Bsum(mi,mj)/Bcount(mi,mj,1);
        
           }
        }
    }
	double dnmc = (double)nmc;
  
   for(int ni = 0; ni < (q-1); ni++){
    for(int nj = (ni+1); nj < q; nj++){
    	
      if(count2(ni,nj,0) > (0.5*nmc)){
      	
        Omegahat(ni,nj) = 0;
        Omegahat(nj,ni) = 0;
      }else{
        Omegahat(ni,nj)  = sum2(ni,nj)/count2(ni,nj,1);
        Omegahat(nj,ni) = Omegahat(ni,nj) ;
	  }
    }
    Omegahat(ni,ni) = sum2(ni,ni)/dnmc;
  }
    Omegahat(q-1,q-1) = sum2(q-1,q-1)/dnmc; 
    
	time_end = clock();
    double time_taken;
    time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
    return(List::create(Named("Bcount") = Bcount,Named("Bsum") = Bsum, Named("B") = out1, Named("Bcumavg") = B_avg, Named("Omega") = out2, Named("Omegacount") = count2, Named("Omegasum") = sum2, Named("Bhat") = Bhat, Named("Omegahat") = Omegahat, Named("time") = time_taken));
       
}
	

// ---------------------------------------------------------------- StepB.cpp
// [[Rcpp::export]]
Rcpp::List StepB(arma::mat B_ini, arma::vec sigma2_ini, arma::mat Y, arma::mat X, int nmc, int burnin, double r1 = 1e-4, double s1 = 1e-8,double alpha1 = 1, double beta1 = 1, int hyp = 0){
	
	int n = Y.n_rows;
	int p = B_ini.n_rows;
	int q = B_ini.n_cols;
	mat XtX = X.t()*X;
	mat YtX = Y.t()*X;
	mat B = zeros(p,q);
	mat Bhat = zeros(p,q);
	B = B_ini ;
	double tau;
	double eta;
	double q1 = 1-(1/(double)p);
	double nb ;
	vec resid =zeros(p);
	double beta ;
	double c1 ;
	double b;
	double c2 ;
	
	
	
	int components = 1;
	vec P_B = zeros(2);
	vec var_mod = sigma2_ini;
	cube Bcount = zeros(p,q,2);
    mat CI_lower = zeros(p,q);
	mat CI_upper = zeros(p,q);
	cube CI = zeros(p,q,2);
	mat Bsum = zeros(p,q);
	cube out = zeros(p,q,nmc+1);
	vec myvec = zeros(2);
	
	clock_t time_start, time_end;
    time_start = clock();
    
    for(int t=0; t<2;t++){
		myvec(t) = t;
	}
    
	for(int it = 0; it < (nmc + burnin+1); it++){
	 	 	
	    for(int k=0; k<q; k++){
			
			uvec idx = find(B.col(k) != 0);
		    nb = idx.size();
	    	resid = Y.col(k) - X*B.col(k);
	    	
	    	if(hyp == 1){
				q1 = R::rbeta((double)p * (double)q - nb + alpha1, nb + beta1);
			}
	    	
           for(int j=0; j<p; j++){
           	
			   b = XtX(j,j) ;
			   c2 = YtX(k,j) - sum(XtX.col(j)%B.col(k)) + B(j,k)*XtX(j,j);
			   
			   if((B(j,k) == 0)){
				eta = R::rgamma(r1, s1);   
				
			   }
			   else{
				eta = R::rgamma(r1 + 0.5, 1/((0.5*B(j,k)*B(j,k)/var_mod(k)) + s1));  
                	
			   }
               
			   if(eta == 0) {
				   P_B(1) = 0;
			   }
			   else {
				   P_B(1) = sqrt(eta/(b + eta)) * exp(c2*c2/(2*(b + eta))) *(1-q1)/q1; 
			   }
			   
			P_B(0) = 1;
            tau = 1/sqrt(eta);
			c1 = b + eta;
			beta =  s1 + 0.5*sum(resid % resid) + sum(B.col(k) % B.col(k))/(2 * tau*tau);
            var_mod(k) = R:: rgamma( r1 + (nb/2) + (n/2), 1/beta);
            var_mod(k) = 1/var_mod(k);
			   
			   if(P_B.has_inf()==true){
			   	
                  int components = 1;
                  B(j,k) = R::rnorm(c2/c1, sqrt(var_mod(k)/c1));
                  
                  if((it > burnin)){
                  	
				     Bcount(j,k,components) += 1;
                     Bsum(j,k) +=  B(j,k);
                     
                     }
                     
                } else {
                	
                	
                	P_B = P_B/sum(P_B);
		   
			        components = (Rcpp::RcppArmadillo::sample(myvec,1,false, P_B) )(0);
			        
			        if(components == 1){
			        	
    		          B(j,k) = R::rnorm(c2/c1, sqrt(var_mod(k)/c1));
    		          
    		          if((it > burnin))
    		          Bsum(j,k) +=  B(j,k);
    		          
          		    } else {
          		    	
            		B(j,k) = 0;
          		    }
          		    
            	    if((it > burnin))
            		Bcount(j,k,components) += 1;
            	}
            	
            	
       		}
       		
       		
   		} 
   		if(it - floor(it/1000)*1000 == 0){
        vec iteration = ones(1)*it;
        iteration.print();
        }
   		
		if(it >= burnin)
		    out.slice(it - burnin) = B	;
		    
		    
	}
	mat Bcount0 = Bcount.slice(0);
	
    
    for (int mi = 0 ; mi < p ; mi++){
       for(int mj = 0 ; mj < q ; mj++){
       	
          if(Bcount(mi,mj,0)> (nmc/2)){
        
             Bhat(mi,mj) = 0;
        
           }else{
              Bhat(mi,mj) = Bsum(mi,mj)/Bcount(mi,mj,1);
        
           }
        }
    }
	
	time_end = clock();
    double time_taken;
    time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
  
	return(List::create(Named("Bcount") = Bcount,Named("Bsum") = Bsum, Named("B") = out, Named("Bhat") = Bhat, Named("time") = time_taken ));
}
	

// ---------------------------------------------------------------- StepOmega.cpp
// [[Rcpp::export]]
Rcpp::List StepOmega(arma::mat S, arma::mat temp, int nmc, int burnin, int n, double r = 1e-4, double s = 1e-8, double alpha2 = 1, double beta2 = 1,int posdef = 0, int hyp = 0){
	
  int p = S.n_cols;
  mat Omega = zeros(p,p);
  Omega = temp;
  mat Omegahat = zeros(p,p);
  cube out = zeros(p,p,nmc+burnin+1); // out is an array of estimated precision matrices in the sampling stage (i.e after burning stage)
  cube OmegaCount = zeros(p,p,2); // OmegaCount is an array of the votes given to zero and non-zero edges
  vec P = zeros(2);
  vec myvec = zeros(2);
  vec mhvec = zeros(2);
  vec probdiag = zeros(2);
  for(int t=0; t<2;t++)
    myvec(t) = t;
  
  double a = 0;
  double b = 0;
  double bd = 0;
  double lambda = 0;
  double q2 = 1-(1/(double)p);
  double v;
  double rho;
  double nomega;
  int rand = 1;
  mat OmegaSum = zeros(p,p);
  mat I = eye(p,p);
    
  clock_t time_start, time_end;
  time_start = clock();
  
  

  for(int it = 0; it < (nmc + burnin+1); it++){
	
    if( hyp == 1){
		uvec idx = find(Omega != 0);
		nomega = 0.5*(idx.size() - (double)p);
	    q2 = R::rbeta(0.5 * (double)p * ((double)p -1) - nomega + alpha2, nomega + beta2);
	}

    for(int j=0; j<(p-1); j++){
      for(int k=(j+1); k<p; k++){
        
        a = S(j,j) + S(k,k);
        b = - sum(Omega.col(j)%S.col(k) + Omega.col(k)%S.col(j)) + Omega(j,k)*a;
		
        if((Omega(j,k) == 0)){
		 lambda = R::rgamma(r, s);
         	 
	    }
		else{
		lambda = R::rgamma(r + 0.5, 1/(0.5*Omega(j,k)*Omega(j,k) + s));   
		
		}
        
        if(lambda == 0)	{
            P(1) = 0;
	    }  
	    else{
        P(1) = sqrt(lambda/(n*a+lambda)) * exp((n*b)*(n*b)/(2*(n*a+lambda))) *(1-q2)/q2; 
	    }
        P(0) = 1;
		
	    if(P.has_inf()==true){
        	
          rand = 1;
          Omega(j,k) = Omega(k,j) = R::rnorm(n*b/(n*a+lambda), 1/sqrt(n*a + lambda));
          
          if((it > burnin)){
          	
          	OmegaCount(j,k,rand) += 1;
			OmegaCount(k,j,rand) =  OmegaCount(j,k,rand);
          	
          	OmegaSum(j,k)  += Omega(j,k);
          	OmegaSum(k,j)  += Omega(k,j);
          	
		    }
            
            
        }  else {
          
          
          P = P/sum(P);
          
          rand = (Rcpp::RcppArmadillo::sample(myvec,1,false, P))(0);
          
          
          if(rand == 1){
            Omega(j,k) = Omega(k,j) = R::rnorm(n*b/(n*a+lambda), 1/sqrt(n*a + lambda));
            
            if((it > burnin)){
            	
            OmegaSum(j,k) +=  Omega(j,k);
            OmegaSum(k,j) += Omega(k,j);
			}
            
    	      
    	       
          } else {
            Omega(j,k) = Omega(k,j) = 0;
          }
     
     
          if((it > burnin))
            OmegaCount(j,k,rand) += 1; 
		    OmegaCount(k,j,rand) =  OmegaCount(j,k,rand);
            
        }

        
      }
      
      bd = sum(Omega.col(j)%S.col(j)) - Omega(j,j)*S(j,j) + (R::rgamma(r + 1, 1/(abs(Omega(j,j)) + s)))/n;
        
	  double mode = (sqrt(bd*bd + 4*S(j,j)) - bd)/(2*S(j,j));
	  v = R::rnorm(mode, 0.001);
	  rho = std::min(exp(n * (log(v) - log(Omega(j,j))) - 0.5 * S(j,j) * ((v*v) - (Omega(j,j)*Omega(j,j))) - bd * (v - Omega(j,j))),1.0);
	  
	  mhvec(0) = v;
	  mhvec(1) = Omega(j,j);
	  probdiag(0) = rho;
	  probdiag(1) = 1-rho;
	  
	  Omega(j,j) = (Rcpp::RcppArmadillo::sample(mhvec,1,false,probdiag))(0);
	
      if(it > burnin)
        OmegaSum(j,j) += Omega(j,j); 
	  
		
    }
    bd = sum(Omega.col(p-1)%S.col(p-1)) - Omega(p-1,p-1)*S(p-1,p-1) + R::rgamma(r + 1, 1/(abs(Omega(p-1, p-1)) + s))/n;
    	
	double mode = (sqrt(bd*bd + 4*S(p-1,p-1)) - bd)/(2*S(p-1,p-1));
	v = R::rnorm(mode, 0.001);
	rho = std::min(exp(n * (log(v) - log(Omega(p-1,p-1))) - 0.5 * S(p-1,p-1) * ((v*v) - (Omega(p-1,p-1)*Omega(p-1,p-1))) - bd * (v - Omega(p-1,p-1))),1.0);
	
	mhvec(0) = v;
	mhvec(1) = Omega(p-1,p-1);
	probdiag(0) = rho;
	probdiag(1) = 1-rho;
		
	Omega(p-1,p-1) = (Rcpp::RcppArmadillo::sample(mhvec,1,false,probdiag))(0);
	
	if(it > burnin)
      OmegaSum(p-1,p-1) += Omega(p-1,p-1);
		
	    if( posdef == 1){
			
            double eigmin = min(eig_sym(Omega));
			
	        if(eigmin < 0){
	            
                Rcout << "Omega is not pd in iteration : " << it <<"\n";	   
	            mat Omeganew = Omega - (eigmin - 0.01) * I;
    	
		        if(it > burnin){
		            for(int j=0; j<p; j++){
		                OmegaSum(j,j) = OmegaSum(j,j) - Omega(j,j) + Omeganew(j,j);
	                }
		        }
	            Omega = Omeganew;
	        }
	    }
	
    
     if(it - floor(it/1000)*1000 == 0){
      vec iteration = ones(1)*it;
      iteration.print();
    }
    if(it >=burnin)
      out.slice(it - burnin) = Omega; 
    
  }
  
  double dnmc = (double)nmc;
  
  for(int ni = 0; ni < (p-1); ni++){
    for(int nj = (ni+1); nj < p; nj++){
    	
      if((OmegaCount(ni,nj,0) > (0.5*nmc))){
      	
        Omegahat(ni,nj) = 0;
        Omegahat(nj,ni) = 0;
      }else{
        Omegahat(ni,nj)  = OmegaSum(ni,nj)/OmegaCount(ni,nj,1);
        Omegahat(nj,ni) = Omegahat(ni,nj) ;
	  }
    }
    Omegahat(ni,ni) = OmegaSum(ni,ni)/dnmc;
  }
  
  Omegahat(p-1,p-1) = OmegaSum(p-1,p-1)/dnmc; 
  
  time_end = clock();
  double time_taken;
  time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
  return(List::create(Named("OmegaCount") = OmegaCount, Named("OmegaSum") = OmegaSum, Named("Omega") = out, Named("Omegahat") = Omegahat, Named("time") = time_taken));
}
