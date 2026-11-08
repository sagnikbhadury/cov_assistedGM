# Gaussian / G-Wishart implementation of the PxG model and its two comparison models.
# This is the authors' PxG-GW.R from the JASA supplement of Niu, Ni, Pati & Mallick (2024).
# Packaging changes only: library()/require() calls removed in favour of namespace-qualified
# calls (BDgraph::gnorm, BDgraph::rgwish, invgamma::rinvgamma), print() progress replaced by
# message() behind a `verbose` argument, roxygen documentation added. No numerical line changed.
# The original file is shipped verbatim in inst/original_supplement/Code/PxG-GW.R.

### This R script contains all functions for the proposed PxG model 
### with Gaussian-G-Wishart specification in section 3.2.1
### and it is used for simulation studies in section 4.1 and 4.2
### to reproduce all results in section 4.1 run R script "Example1.R"
### to reproduce all results in section 4.2 run R script "Example2.R"

### There are three main functions in this R script: PxG_GW, CovOnly_GW, GraphOnly_GW
### "PxG_GW": the proposed PxG model
### "CovOnly_GW": covariate-only model where graph is assumed to be constant among samples
### "GraphOnly_GW": graph-only model where graph varies with samples but does not depend on covariates
### For all three functions, Gaussian-G-Wishart model is used for likelihood and graphical prior 

### list of arguments for all three functions
### Y: response data or graphical data (n by q)
### X: covariate data (n by p)
### b: degree of freedom for G-Wishart distribution
### D: positive definite matrix for G-Wishart distribution (q by q)
### alphaG: prior edge inclusion probability
### mu0: mean of normal priors in the auxiliary model, see section 3.1
### sigsq0: variance of priors in the auxiliary model, see section 3.1
### b1: shape parameter of inverse gamma prior in the auxiliary model, see section 3.1
### b2: rate parameter of inverse gamma prior in the auxiliary model, see section 3.1
### M: scaling parameter in the stick-breaking process
### N: maximum number of clusters allowed
### burnin: number of burn-in samples
### nmcmc: number of MCMC samples
### z0: initial cluster allocation
### B: number of graph updates in each iteration

### list of values for all three functions
### Label: posterior estimates of cluster labels using the least squares clustering method (Dahl, 2006) 
### Gsum: posterior edge inclusion probability of graph for each sample
### Osum: posterior precision matrix for each sample
### avgdist: average least square distance to the mean for all MCMC samples of cluster labels (Dahl, 2006) 
### DIC: DIC value
### LogTotLike: log total likelihood under the MCMC sample (the one in "Label") 
###             selected by the least squares clustering method (Dahl, 2006) 
### VarTotLike: variance of log total likelihood under all MCMC samples of cluster labels





UpdateCovParam = function(z, N, X, mu0, sigsq0, b1, b2){
  p = ncol(X)
  xMu = matrix(0, p, N)
  xSigSq = rep(1, N)
  
  for(i in 1:N){
    indi = which(z==i)
    ni = length(indi)
    Xi = X[indi, , drop = FALSE]
    
    if(ni==0){
      xSigSq[i] = invgamma::rinvgamma(1, shape = b1, rate = b2)
      xMu[, i] = mu0 + sqrt(sigsq0*xSigSq[i]) * rnorm(p)
    }
    else{
      postb1 = b1 + 0.5*ni*p;
      munew = apply(Xi, 2, sum) + mu0/sigsq0
      sigsqnew = sigsq0/(ni*sigsq0+1)
      TSS = sum(Xi^2) + sum(mu0^2)/sigsq0 - sum(munew^2)*sigsqnew
      postb2 = b2 + 0.5*TSS
      xSigSq[i] = invgamma::rinvgamma(1, shape = postb1, rate = postb2)
      
      postMu = munew*sigsqnew
      postSigSq = xSigSq[i]*sigsqnew
      xMu[,i] = postMu + sqrt(postSigSq) * rnorm(p)
    }
  }
  
  return(list(xMu = xMu, xSigSq = xSigSq))
}


CovlogLike = function(x, mu, sigsq){
  p = length(x)
  tss = sum((x-mu)^2)
  return(-0.5*p*log(sigsq)-tss/(2*sigsq))
}


DrawGraph = function(alphaG, q){
  Adj = diag(rep(1,q))
  edges = rbinom(q*(q-1)/2, 1, alphaG)
  Adj[lower.tri(Adj, diag = FALSE)] = edges
  Adj = t(Adj) + Adj - diag(rep(1,q))
  return(Adj)
}


LogMarginG = function(S, Adj, n, b, D){
  logconst1 = BDgraph::gnorm(Adj, b+n, D+S, 100)
  logconst2 = BDgraph::gnorm(Adj, b, D, 100)
  return( logconst1 - logconst2 )
}


UpdateGraphParam0 = function(z, N, Y, b, D, alphaG, B){
  q = ncol(Y)
  ListG = ListOmg = array(diag(q), c(q,q,N))
  upperind = which(upper.tri(diag(q))==TRUE)
  
  for(i in 1:N){
    indi = which(z==i)
    ni = length(indi)
    
    if(ni==0){
      ListG[,,i] = DrawGraph(alphaG, q)
      ListOmg[,,i] = BDgraph::rgwish(n=1, adj=ListG[,,i], b=b, D=D)
    }
    else{
      Yi = Y[indi,,drop = FALSE]
      Si = t(Yi)%*%Yi
      Adj = diag(q)
      for(round in 1:B){
        
        for(j in upperind){
          Adjt = Adj 
          if(Adjt[j] == 0){
            Adjt[j] = 1
            Adjt = t(Adjt)
            Adjt[j] = 1
            prop = alphaG/(1-alphaG)
          }
          else{
            Adjt[j] = 0
            Adjt = t(Adjt)
            Adjt[j] = 0
            prop = (1-alphaG)/alphaG
          }
          r = LogMarginG(Si, Adjt, ni, b, D) - LogMarginG(Si, Adj, ni, b, D) + log(prop)
          if(r>log(runif(1))){ Adj = Adjt }
        }
        
      }
      ListG[,,i] = Adj; ListOmg[,,i] = BDgraph::rgwish(n=1, adj=ListG[,,i], b=b+ni, D=D+Si)
    }
  }
  
  return(list(ListG=ListG, ListOmg = ListOmg))
}


GraphUpdate = function(Y, Adj0, b, D, alphaG, B){
  Adj = Adjt = Adj0
  n = nrow(Y)
  q = ncol(Y)
  S = t(Y)%*%Y
  upperind = which(upper.tri(Adj)==TRUE)
  
  for(rounds in 1:B){
    for(flipind in upperind){
      Adjt = Adj
      
      if(Adjt[flipind] == 0){
        Adjt[flipind] = 1
        Adjt = t(Adjt)
        Adjt[flipind] = 1
        prop = alphaG/(1-alphaG)
      }
      else{
        Adjt[flipind] = 0
        Adjt = t(Adjt)
        Adjt[flipind] = 0
        prop = (1-alphaG)/alphaG
      }
      
      r = LogMarginG(S, Adjt, n, b, D) - LogMarginG(S, Adj, n, b, D) + log(prop)
      if(r>log(runif(1))){ Adj = Adjt }
    }
  }
  
  return(Adj)
}


UpdateGraphParam = function(z, N, Y, b, D, alphaG, B){
  q = ncol(Y)
  ListG = ListOmg = array(diag(q), c(q,q,N))
  
  for(i in 1:N){
    indi = which(z==i)
    ni = length(indi)
    if(ni==0){
      ListG[,,i] = DrawGraph(alphaG, q)
      ListOmg[,,i] = BDgraph::rgwish(1, ListG[,,i], b, D)
    }
    else{
      Yi = Y[indi, , drop = FALSE]
      ListG[,,i] = GraphUpdate(Yi, diag(q), b, D, alphaG, B)
      ListOmg[,,i] = BDgraph::rgwish(n=1, adj=ListG[,,i], b=b+ni, D=D+t(Yi)%*%Yi)
    }
  }
  
  return(list(ListG=ListG, ListOmg=ListOmg))
}


getMemberMat = function(z){
  n = length(z)
  Membership = matrix(0, n, n)
  for(cls in unique(z)){
    ind = which(z==cls)
    Membership[ind, ind] = 1
  }
  diag(Membership) = 0
  return(Membership)
}


getlogTotLikeX = function(X, z, mu0, sigsq0, b1, b2){
  zs = unique(z)
  logTotX = 0
  for(i in zs){
    id = which(z==i)
    lgoTotX = logTotX + logMultiX(X[id, , drop=FALSE], mu0, sigsq0, b1, b2)
  }
  return(logTotX)
}


getLogLikeY = function(Y, Adj0, b, D, alphaG, nmc, burnin, nmcmc){
  Adj = Adj0
  logTotY = rep(0, nmcmc)
  for(i in 1:(burnin+nmcmc)){
    Adj = GraphUpdate(Y, Adj0, b, D, alphaG, 1)
    if(i>burnin){ logTotY[i-burnin] = logMarginLikeY(Y, Adj, b, D, 100) }
  }
  return(logTotY)
}


#' Covariate-only model, Gaussian / G-Wishart specification
#'
#' Comparison model of Niu et al. (2024) in which the partition is driven by the covariates
#' alone and a single graph is shared by all samples. The DIC combines the covariate-model
#' likelihood at the selected partition with a separate single-graph MCMC run on \code{Y}.
#'
#' Original authors' code; see \code{\link{PxG_GW}} for the argument conventions.
#'
#' @param Y Response (graphical) data, an \code{n x q} numeric matrix.
#' @param X Covariate data, an \code{n x p} numeric matrix.
#' @param b Degrees of freedom of the G-Wishart prior.
#' @param D Positive definite \code{q x q} scale matrix of the G-Wishart prior.
#' @param alphaG Prior edge inclusion probability.
#' @param mu0 Prior mean (length \code{p}) of the cluster means in the covariate model (paper Section 3.1).
#' @param sigsq0 Prior variance scale of the cluster means in the covariate model.
#' @param b1,b2 Shape and rate of the inverse-gamma prior on the covariate variance.
#' @param M Concentration parameter of the truncated stick-breaking (Dirichlet process) prior.
#' @param N Maximum number of clusters (truncation level).
#' @param burnin Number of burn-in iterations.
#' @param nmcmc Number of MCMC iterations kept after burn-in.
#' @param z0 Initial cluster allocation, an integer vector of length \code{n} with values in \code{1:N}.
#' @param B Number of graph (edge-flip sweep) updates per iteration.
#' @param verbose Logical; print per-iteration progress as in the original scripts. Default \code{TRUE}.
#' @return A list with \code{Label}, \code{avgdist}, \code{DIC}, \code{LogTotLike} and \code{VarTotLike}
#'   as in \code{\link{PxG_GW}} (no per-sample graphs are returned).
#' @references Niu, Y., Ni, Y., Pati, D. and Mallick, B. K. (2024). Covariate-Assisted Bayesian
#'   Graph Learning for Heterogeneous Data. \emph{Journal of the American Statistical Association},
#'   119(547), 1985-1999. \doi{10.1080/01621459.2023.2233744}
#'
#'   Dahl, D. B. (2006). Model-based clustering for expression data via a Dirichlet process mixture model.
#'   In \emph{Bayesian Inference for Gene Expression and Proteomics}, Cambridge University Press.
#' @seealso \code{\link{PxG_GW}}, \code{\link{GraphOnly_GW}}
#' @export
CovOnly_GW = function(Y, X, b, D, alphaG, mu0, sigsq0, b1, b2, M, N, burnin, nmcmc, z0, B, verbose = TRUE){
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  Z = matrix(0, nrow = nmcmc, ncol = n)
  Membership = array(0, c(n,n,nmcmc))
  SUM = matrix(0, n, n)
  z = z0
  
  V = rep(0, N); V[N] = 1
  prob = rep(0, N); nc = rep(0, N)
  weights = rep(0, N)
  
  xMu = matrix(0, p, N)
  xSigSq = rep(1, N)

  vlogTotLikeX = rep(0, nmcmc)
  
  for(iter in 1:(burnin+nmcmc)){
    if(verbose) message("iteration i = ", iter)
    
    for(i in 1:N){ zi = which(z==i); nc[i] = length(zi) }
    for(i in 1:(N-1)){ V[i] = rbeta(1, 1+nc[i], M+sum(nc[(i+1):N])) }
    Vp = 1 - V; prob[1] = V[1]
    for(i in 2:N){ prob[i] = V[i]*prod(Vp[1:(i-1)]) }

    for(i in 1:n){
      for(j in 1:N){ 
        weights[j] = log(prob[j])
        weights[j] = weights[j] + CovlogLike(X[i,,drop=FALSE], xMu[,j], xSigSq[j])
      }

      maxweight = max(weights)
      weights = weights - maxweight
      weights = exp(weights)
      weights = weights/sum(weights)
      z[i] = sample(1:N, 1, prob = weights)
    }

    if(verbose) message("unique z = ", length(unique(z)))
    Cov = UpdateCovParam(z, N, X, mu0, sigsq0, b1, b2)
    xMu = Cov$xMu; xSigSq = Cov$xSigSq
    
    if(iter>burnin){
      Z[iter-burnin, ] = z
      Membership[,, iter-burnin] = getMemberMat(as.vector(z))
      SUM = SUM + Membership[,, iter-burnin]
      
      vlogTotLikeX[iter-burnin] = getlogTotLikeX(X, z, mu0, sigsq0, b1, b2)
    }
  }
  
  if(verbose) message("Calculating Membership Matrix...")
  SUM = SUM/nmcmc
  dist = rep(0, nmcmc)
  for(i in 1:nmcmc){ dist[i] = sum((Membership[,,i]-SUM)^2) }
  minZ = which.min(dist)
  result = Z[minZ,]

  if(verbose) message("Calculating DIC...")
  vlogTotLikeY = getLogLikeY(Y, diag(q), b, D, alphaG, 100, burnin, nmcmc)
  LogTotLike = vlogTotLikeX[minZ] + max(vlogTotLikeY)
  VarTotLike = var(vlogTotLikeX) + var(vlogTotLikeY)
  DIC = -2*vlogTotLikeX[minZ] + var(vlogTotLikeX) - 2*max(vlogTotLikeY) + var(vlogTotLikeY)
  
  return(list(Label = result, avgdist = SUM,
              DIC = DIC, LogTotLike = LogTotLike, VarTotLike = VarTotLike))
}


getlogTotLikeY = function(Y, ListG, z, b, D, nmc){
  zs = unique(z)
  logTotY =  0
  for(i in zs){
    id = which(z==i)
    logTotY = logTotY + logMarginLikeY(Y[id, , drop=FALSE], ListG[ , , i], b, D, nmc)
  }
  return(logTotY)
}


#' Graph-only model, Gaussian / G-Wishart specification
#'
#' Comparison model of Niu et al. (2024) in which the partition is driven by the graphical
#' data alone (a Dirichlet process mixture of Gaussian graphical models); the covariates enter
#' the DIC only through a single-cluster marginal likelihood.
#'
#' Original authors' code; see \code{\link{PxG_GW}} for the argument conventions.
#'
#' @param Y Response (graphical) data, an \code{n x q} numeric matrix.
#' @param X Covariate data, an \code{n x p} numeric matrix.
#' @param b Degrees of freedom of the G-Wishart prior.
#' @param D Positive definite \code{q x q} scale matrix of the G-Wishart prior.
#' @param alphaG Prior edge inclusion probability.
#' @param mu0 Prior mean (length \code{p}) of the cluster means in the covariate model (paper Section 3.1).
#' @param sigsq0 Prior variance scale of the cluster means in the covariate model.
#' @param b1,b2 Shape and rate of the inverse-gamma prior on the covariate variance.
#' @param M Concentration parameter of the truncated stick-breaking (Dirichlet process) prior.
#' @param N Maximum number of clusters (truncation level).
#' @param burnin Number of burn-in iterations.
#' @param nmcmc Number of MCMC iterations kept after burn-in.
#' @param z0 Initial cluster allocation, an integer vector of length \code{n} with values in \code{1:N}.
#' @param B Number of graph (edge-flip sweep) updates per iteration.
#' @param verbose Logical; print per-iteration progress as in the original scripts. Default \code{TRUE}.
#' @return A list with components
#' \describe{
#'   \item{Label}{posterior cluster labels chosen by Dahl's (2006) least-squares criterion}
#'   \item{Gsum}{\code{q x q x n} array of posterior edge inclusion probabilities, one slice per sample}
#'   \item{Osum}{\code{q x q x n} array of posterior mean precision matrices, one slice per sample}
#'   \item{avgdist}{\code{n x n} posterior co-clustering (pairwise allocation) matrix}
#'   \item{DIC}{deviance information criterion}
#'   \item{LogTotLike}{log total likelihood at the selected MCMC sample}
#'   \item{VarTotLike}{variance of the log total likelihood over the kept samples}
#' }
#' @references Niu, Y., Ni, Y., Pati, D. and Mallick, B. K. (2024). Covariate-Assisted Bayesian
#'   Graph Learning for Heterogeneous Data. \emph{Journal of the American Statistical Association},
#'   119(547), 1985-1999. \doi{10.1080/01621459.2023.2233744}
#'
#'   Dahl, D. B. (2006). Model-based clustering for expression data via a Dirichlet process mixture model.
#'   In \emph{Bayesian Inference for Gene Expression and Proteomics}, Cambridge University Press.
#' @seealso \code{\link{PxG_GW}}, \code{\link{CovOnly_GW}}
#' @export
GraphOnly_GW = function(Y, X, b, D, alphaG, mu0, sigsq0, b1, b2, M, N, burnin, nmcmc, z0, B, verbose = TRUE){
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  Z = matrix(0, nrow = nmcmc, ncol = n)
  Membership = array(0, c(n,n,nmcmc))
  SUM = matrix(0, n, n)
  z = z0
  
  V = rep(0, N); V[N] = 1
  prob = rep(0, N); nc = rep(0, N)
  weights = rep(0, N)
  
  ListG = array(diag(q), c(q, q, N))
  ListOmg = array(diag(q), c(q, q, N))
  postG = array(0, c(q,q,n))
  postOmg = array(0, c(q,q,n))

  vlogTotLikeY = rep(0, nmcmc)
  
  for(iter in 1:(burnin+nmcmc)){
    if(verbose) message("iteration i = ", iter)
    
    for(i in 1:N){ zi = which(z==i); nc[i] = length(zi) }
    for(i in 1:(N-1)){ V[i] = rbeta(1, 1+nc[i], M+sum(nc[(i+1):N])) }
    Vp = 1 - V; prob[1] = V[1]
    for(i in 2:N){ prob[i] = V[i]*prod(Vp[1:(i-1)]) }
    
    for(i in 1:n){
      for(j in 1:N){ 
        weights[j] = log(prob[j]) + 0.5*log(det(ListOmg[,,j])) - 0.5*Y[i,, drop = FALSE]%*%ListOmg[,,j]%*%t(Y[i,,drop = FALSE]) 
        weights[j] = weights[j]
      }
      
      maxweight = max(weights)
      weights = weights - maxweight
      weights = exp(weights)
      weights = weights/sum(weights)
      z[i] = sample(1:N, 1, prob = weights)
    }

    if(verbose) message("unique z = ", length(unique(z)))
    G = UpdateGraphParam(z, N, Y, b, D, alphaG, B)
    ListG = G$ListG; ListOmg = G$ListOmg
    
    if(iter>burnin){
      Z[iter-burnin, ] = z
      Membership[,, iter-burnin] = getMemberMat(as.vector(z))
      SUM = SUM + Membership[,, iter-burnin]

      vlogTotLikeY[iter-burnin] = getlogTotLikeY(Y, ListG, z, b, D, 100)
      
      for(k in 1:n){
        postG[,, k] = postG[,, k] + ListG[,, z[k]]
        postOmg[,, k] = postOmg[,, k] + ListOmg[,, z[k]]
      }
      
    }
  }
  
  if(verbose) message("Calculating Membership Matrix...")
  SUM = SUM/nmcmc
  dist = rep(0, nmcmc)
  for(i in 1:nmcmc){ dist[i] = sum((Membership[,,i]-SUM)^2) }
  minZ = which.min(dist)
  result = Z[minZ,]

  if(verbose) message("Calculating DIC...")
  LogTotLikeX = logMultiX(X, mu0, sigsq0, b1, b2)
  LogTotLike = vlogTotLikeY[minZ] + LogTotLikeX
  VarTotLike = var(vlogTotLikeY)
  DIC = -2*vlogTotLikeY[minZ] + var(vlogTotLikeY) - 2*LogTotLikeX
  
  return(list(Label = result, Gsum = postG/nmcmc, Osum = postOmg/nmcmc, avgdist = SUM,
              DIC = DIC, LogTotLike = LogTotLike, VarTotLike = VarTotLike))
}


logMarginLikeY = function(Y, Adj, b, D, nmc){
  n = nrow(Y)
  q = ncol(Y)
  S = t(Y)%*%Y
  logconst1 = BDgraph::gnorm(Adj, b+n, D+S, nmc)
  logconst2 = BDgraph::gnorm(Adj, b, D, nmc)
  return( -0.5*q*n*log(2*pi) + logconst1 - logconst2 )
}


logMultiX = function(X, mu0, sigsq0, b1, b2){
  n = nrow(X)
  p = ncol(X)
  if(n==1){
    AX = sum((X-mu0)*(X-mu0))*sigsq0/(1+sigsq0)
  }
  else{
    a = mu0/sigsq0 + as.vector(apply(X,2,sum))
    AX = sum(X*X) + sum(mu0*mu0)/sigsq0 - sigsq0/(1+n*sigsq0)*sum(a*a)
  }
  temp1 = -0.5*p*n*log(2*pi) - 0.5*log(det(diag(p)+n*sigsq0))
  temp2 = lgamma(b1+0.5*p*n) - (b1+0.5*p*n)*log(b2+0.5*AX)

  return(temp1 + temp2 + b1*log(b2) - lgamma(b1))
}


getlogTotLike = function(Y, X, ListG, z, b, D, nmc, mu0, sigsq0, b1, b2){ 
  zs = unique(z)
  logTotY = logTotX = 0
  for(i in zs){
    id = which(z==i)
    logTotY = logTotY + logMarginLikeY(Y[id, , drop=FALSE], ListG[ , , i], b, D, nmc)
    lgoTotX = logTotX + logMultiX(X[id, , drop=FALSE], mu0, sigsq0, b1, b2)
  }
  return(logTotY + logTotX)
}


#' PxG model, Gaussian / G-Wishart specification
#'
#' Blocked Gibbs sampler for the covariate-assisted partition-and-graph (PxG) model of
#' Niu et al. (2024), Section 3.2.1: a truncated Dirichlet process partitions the samples,
#' each cluster carries its own Gaussian graphical model with a G-Wishart prior on the precision
#' matrix, and the covariates enter through a Gaussian auxiliary model so that the partition is
#' informed by both the responses and the covariates. G-Wishart normalising constants are
#' approximated with \code{BDgraph::gnorm}, so this specification is intended for small \code{q}
#' (the paper uses \code{q = 3} and \code{q = 5}); use \code{\link{PxG_PL}} for larger graphs.
#'
#' The function body is the authors' code from the JASA supplement (\code{PxG-GW.R}), unchanged
#' apart from namespace-qualified calls and the \code{verbose} switch. See the package README,
#' section "Notes on the original code", for remarks on this specification.
#'
#' @param Y Response (graphical) data, an \code{n x q} numeric matrix.
#' @param X Covariate data, an \code{n x p} numeric matrix.
#' @param b Degrees of freedom of the G-Wishart prior.
#' @param D Positive definite \code{q x q} scale matrix of the G-Wishart prior.
#' @param alphaG Prior edge inclusion probability.
#' @param mu0 Prior mean (length \code{p}) of the cluster means in the covariate model (paper Section 3.1).
#' @param sigsq0 Prior variance scale of the cluster means in the covariate model.
#' @param b1,b2 Shape and rate of the inverse-gamma prior on the covariate variance.
#' @param M Concentration parameter of the truncated stick-breaking (Dirichlet process) prior.
#' @param N Maximum number of clusters (truncation level).
#' @param burnin Number of burn-in iterations.
#' @param nmcmc Number of MCMC iterations kept after burn-in.
#' @param z0 Initial cluster allocation, an integer vector of length \code{n} with values in \code{1:N}.
#' @param B Number of graph (edge-flip sweep) updates per iteration.
#' @param verbose Logical; print per-iteration progress as in the original scripts. Default \code{TRUE}.
#' @return A list with components
#' \describe{
#'   \item{Label}{posterior cluster labels chosen by Dahl's (2006) least-squares criterion}
#'   \item{Gsum}{\code{q x q x n} array of posterior edge inclusion probabilities, one slice per sample}
#'   \item{Osum}{\code{q x q x n} array of posterior mean precision matrices, one slice per sample}
#'   \item{avgdist}{\code{n x n} posterior co-clustering (pairwise allocation) matrix}
#'   \item{DIC}{deviance information criterion}
#'   \item{LogTotLike}{log total likelihood at the selected MCMC sample}
#'   \item{VarTotLike}{variance of the log total likelihood over the kept samples}
#' }
#' @references Niu, Y., Ni, Y., Pati, D. and Mallick, B. K. (2024). Covariate-Assisted Bayesian
#'   Graph Learning for Heterogeneous Data. \emph{Journal of the American Statistical Association},
#'   119(547), 1985-1999. \doi{10.1080/01621459.2023.2233744}
#'
#'   Dahl, D. B. (2006). Model-based clustering for expression data via a Dirichlet process mixture model.
#'   In \emph{Bayesian Inference for Gene Expression and Proteomics}, Cambridge University Press.
#' @seealso \code{\link{CovOnly_GW}}, \code{\link{GraphOnly_GW}}, \code{\link{PxG_PL}}
#' @examples
#' \donttest{
#' ## a short run on the 3-node piecewise example (paper Section 4.1 uses burnin 500, nmcmc 1000)
#' sim <- simulate_example1(seed = 2020)
#' set.seed(2020)
#' fit <- PxG_GW(Y = sim$Y, X = sim$X, b = 3, D = diag(3), alphaG = 0.1,
#'               mu0 = 0, sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5,
#'               burnin = 20, nmcmc = 50, z0 = sample(1:5, 300, TRUE), B = 1, verbose = FALSE)
#' table(fit$Label)
#' }
#' @export
PxG_GW = function(Y, X, b, D, alphaG, mu0, sigsq0, b1, b2, M, N, burnin, nmcmc, z0, B, verbose = TRUE){
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  Z = matrix(0, nrow = nmcmc, ncol = n)
  Membership = array(0, c(n,n,nmcmc))
  SUM = matrix(0, n, n)
  z = z0
  
  V = rep(0, N); V[N] = 1
  prob = rep(0, N); nc = rep(0, N)
  weights = rep(0, N)
  
  xMu = matrix(0, p, N)
  xSigSq = rep(1, N)
  ListG = array(diag(q), c(q, q, N))
  ListOmg = array(diag(q), c(q, q, N))
  postG = array(0, c(q,q,n))
  postOmg = array(0, c(q,q,n))

  vlogTotLike = rep(0, nmcmc)
  
  for(iter in 1:(burnin+nmcmc)){
    if(verbose) message("iteration i = ", iter)
    
    for(i in 1:N){ zi = which(z==i); nc[i] = length(zi) }
    for(i in 1:(N-1)){ V[i] = rbeta(1, 1+nc[i], M+sum(nc[(i+1):N])) }
    Vp = 1 - V; prob[1] = V[1]
    for(i in 2:N){ prob[i] = V[i]*prod(Vp[1:(i-1)]) }
    
    for(i in 1:n){
      for(j in 1:N){ 
        weights[j] = log(prob[j]) + 0.5*log(det(ListOmg[,,j])) - 0.5*Y[i,, drop = FALSE]%*%ListOmg[,,j]%*%t(Y[i,,drop = FALSE]) 
        weights[j] = weights[j] + CovlogLike(X[i,,drop=FALSE], xMu[,j], xSigSq[j])
      }
      
      maxweight = max(weights)
      weights = weights - maxweight
      weights = exp(weights)
      weights = weights/sum(weights)
      z[i] = sample(1:N, 1, prob = weights)
    }

    if(verbose) message("unique z = ", length(unique(z)))
    G = UpdateGraphParam(z, N, Y, b, D, alphaG, B)
    ListG = G$ListG; ListOmg = G$ListOmg
    Cov = UpdateCovParam(z, N, X, mu0, sigsq0, b1, b2)
    xMu = Cov$xMu; xSigSq = Cov$xSigSq
    
    if(iter>burnin){
      Z[iter-burnin, ] = z
      Membership[,, iter-burnin] = getMemberMat(as.vector(z))
      SUM = SUM + Membership[,, iter-burnin]

      vlogTotLike[iter-burnin] = getlogTotLike(Y, X, ListG, z, b, D, 100, mu0, sigsq0, b1, b2)
      
      for(k in 1:n){
        postG[,, k] = postG[,, k] + ListG[,, z[k]]
        postOmg[,, k] = postOmg[,, k] + ListOmg[,, z[k]]
      }
      
    }
  }
  
  if(verbose) message("Calculating Membership Matrix...")
  SUM = SUM/nmcmc
  dist = rep(0, nmcmc)
  for(i in 1:nmcmc){ dist[i] = sum((Membership[,,i]-SUM)^2) }
  minZ = which.min(dist)
  result = Z[minZ,]

  if(verbose) message("Calculating DIC...")
  LogTotLike = vlogTotLike[minZ]
  VarTotLike = var(vlogTotLike)
  DIC = -2*vlogTotLike[minZ] + var(vlogTotLike)

  return(list(Label = result, Gsum = postG/nmcmc, Osum = postOmg/nmcmc, avgdist = SUM,
              DIC = DIC, LogTotLike = LogTotLike, VarTotLike = VarTotLike))
}


#' Posterior edge inclusion probabilities for a single Gaussian graphical model
#'
#' Single-cluster G-Wishart graph sampler (edge-flip Metropolis moves with
#' \code{BDgraph::gnorm} marginal likelihoods) used by the authors as a helper.
#'
#' @param Y Data matrix, \code{n x q}.
#' @param Adj0 Initial \code{q x q} adjacency matrix (ones on the diagonal).
#' @param b,D G-Wishart degrees of freedom and scale matrix.
#' @param alphaG Prior edge inclusion probability.
#' @param burnin,nmcmc Burn-in and kept iterations.
#' @return A \code{q x q} matrix of posterior edge inclusion frequencies.
#' @export
GraphEst = function(Y, Adj0, b, D, alphaG, burnin, nmcmc){
  Adj = Adjt = Adj0
  n = nrow(Y)
  q = ncol(Y)
  S = t(Y)%*%Y
  upperind = which(upper.tri(Adj)==TRUE)
  SUM = matrix(0, q, q)
  
  for(iter in 1:(burnin+nmcmc)){
    for(flipind in upperind){
      Adjt = Adj
      
      if(Adjt[flipind] == 0){
        Adjt[flipind] = 1
        Adjt = t(Adjt)
        Adjt[flipind] = 1
        prop = alphaG/(1-alphaG)
      }
      else{
        Adjt[flipind] = 0
        Adjt = t(Adjt)
        Adjt[flipind] = 0
        prop = (1-alphaG)/alphaG
      }
      
      r = LogMarginG(S, Adjt, n, b, D) - LogMarginG(S, Adj, n, b, D) + log(prop)
      if(r>log(runif(1))){ Adj = Adjt }
    }
    
    if(iter>burnin){ SUM = SUM + Adj }
  }
  SUM = SUM/nmcmc
  return(SUM)
}

