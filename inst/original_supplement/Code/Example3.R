### This R script reproduces section 4.3 "Example 3: piecewise constant precision matrix"
### including Figure 7: posterior clustering results for all three models
### Figure 8: posterior edge inclusion probability for graphs in both clusters under PxG
### Figure 9: posterior results regarding graph estimates under BM13, mSSL and JRNS

rm(list = ls())

## loading libraries 
library(BDgraph)
library(invgamma)

## generating simulation data from the 50-dimensional GGM in section 4.3
set.seed(2022)
p = 10
q = 50
n = 500
X = matrix(0, n, p)
for(i in 1:(n/2)){
  X[i,] = rnorm(p, 0, 1)
}
for(i in (n/2+1):n){
  X[i,] = rnorm(p, 2, 1)
}

Y = matrix(0, n, q)
Omega1 = matrix(0, q, q)
Omega1[upper.tri(Omega1, diag = FALSE)] = rbinom(q*(q-1)/2, 1, 0.0085)
Omega1 = Omega1 + t(Omega1); diag(Omega1) = 1; G1 = Omega1
Omega2 = matrix(0, q, q)
Omega2[upper.tri(Omega2, diag = FALSE)] = rbinom(q*(q-1)/2, 1, 0.0085)
Omega2 = Omega2 + t(Omega2); diag(Omega2) = 1; G2 = Omega2
Omega1 = rgwish(1,Omega1,3,diag(q))
Omega2 = rgwish(1,Omega2,3,diag(q))
Sigma1 = solve(Omega1)
Sigma2 = solve(Omega2)
for(i in 1:(n/2)){ Y[i,] = t(t(chol(Sigma1))%*%rnorm(q,0,1)) }
for(i in (n/2+1):n){ Y[i,] = t(t(chol(Sigma2))%*%rnorm(q,0,1)) }



## loading libraries and compiling all C++ functions
library(Rcpp)
library(RcppArmadillo)
sourceCpp("PxG-PL.cpp")

## fitting simulation data with the proposed PxG model
set.seed(2022)
A = PxG_PL(YData = Y, XData = X, 
           spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1, 
           mu0 = rep(0, p), sigsq0 = 1, b1 = 1, b2 = 0.1, 
           M = 1, N = 3, burnin = 500, nmcmc = 1000, z0 = sample(1:3,n,TRUE, c(1/3,1/3,1/3)), B = 1)
A$DIC # 104164.7
save(A, file = "Example3A.RData")

## fitting simulation data with the covariate-only model
set.seed(2022)
B = CovOnly_PL(YData = Y, XData = X, 
               spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1, 
               mu0 = rep(0, p), sigsq0 = 1, b1 = 1, b2 = 0.1, 
               M = 1, N = 3, burnin = 500, nmcmc = 1000, z0 = sample(1:3,n,TRUE, c(1/3,1/3,1/3)))
B$DIC # 155453.3
save(B, file = "Example3B.RData")

## fitting simulation data with the graph-only model
set.seed(2022)
C = GraphOnly_PL(YData = Y, XData = X, 
                 spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1, 
                 mu0 = rep(0, p), sigsq0 = 1, b1 = 1, b2 = 0.1, 
                 M = 1, N = 3, burnin = 500, nmcmc = 1000, z0 = sample(1:3,n,TRUE, c(1/3,1/3,1/3)), B = 1)
C$DIC # 106138.6
save(C, file = "Example3C.RData")



## loading libraries for plot
library(ggplot2)
library(ggpubr)

## plot Figure 7
subX = as.data.frame(X[,1:2])
subX[501,] = c(NA, NA)
subX$Cluster = c(rep("Cluster 1",n/2), rep("Cluster 2",n/2), rep("Cluster 3", 1))
cols <- c("#1170AA", "#EF6F6A", "#009E73")
fg1 = ggplot(subX, aes(x = V1, y = V2, color = Cluster,shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=12, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = expression(x[1]), y = expression(x[2]), caption="(A) True cluster allocation")

subX = as.data.frame(X[,1:2])
Bind = rep(0,500)
Bind[which(B$Label==2)] = "Cluster 1"
Bind[which(B$Label==1)] = "Cluster 2"
Bind[which(B$Label==3)] = "Cluster 3"
subX$Cluster = Bind
cols <- c("#1170AA", "#EF6F6A", "#009E73")
fg2 = ggplot(subX, aes(x = V1, y = V2, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=12, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = expression(x[1]), y = expression(x[2]), caption="(B) Covariate-only model")

subX = as.data.frame(X[,1:2])
Cind = rep(0,500)
Cind[which(C$Label==2)] = "Cluster 1"
Cind[which(C$Label==1)] = "Cluster 2"
Cind[which(C$Label==3)] = "Cluster 3"
subX$Cluster = Cind
cols <- c("#1170AA", "#EF6F6A", "#009E73")
fg3 = ggplot(subX, aes(x = V1, y = V2, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=12, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = expression(x[1]), y = expression(x[2]), caption="(C) Graph-only model")

fg0 = ggarrange(fg1, fg2, fg3, ncol = 3, nrow = 1, common.legend = TRUE); fg0



## too save time, uncomment the next line to load results from the PxG model
## load("Example3A.RData")
## calculating graph estimate for each cluster under the PxG model
set.seed(2022)
AG1 = Graph_NS(YData = Y[which(A$Label==2),], spike = 0.01, slab = 5, 
               a1 = 1, a2 = 1, alphaG = 0.05, burnin = 100, nmcmc = 500)

set.seed(2022)
AG2 = Graph_NS(YData = Y[which(A$Label==3),], spike = 0.01, slab = 5, 
               a1 = 1, a2 = 1, alphaG = 0.05, burnin = 100, nmcmc = 500)

save(AG1, AG2, file = "Example3AG.RData")



## making pseudo-likelihood graph estimates symmetric
AG1new = AG1
for(i in 1:50){
  for(j in 1:50){
    AG1new[i,j] = max(AG1[i,j], AG1[j,i])
  }
}
diag(AG1new) = 1

AG2new = AG2
for(i in 1:50){
  for(j in 1:50){
    AG2new[i,j] = max(AG2[i,j], AG2[j,i])
  }
}
diag(AG2new) = 1

xaxis = as.character(1:50)
yaxis = xaxis[50:1]
data = expand.grid(X=xaxis, Y=yaxis)
diag(G1) = 0
diag(G2) = 0

## plot Figure 8
G1T = G1F = matrix(0,50,50)
for(i in 1:50){
  for(j in 1:50){
    if( (G1[i,j]==1)&(AG1new[i,j]>=0.5) ){ G1T[i,j] = 1 }
    if( (G1[i,j]==1)&(AG1new[i,j]<0.5) ){ G1F[i,j] = 1 }
  }
}
gg1 = AG1new[50:1,]
data$freq = as.vector(t(gg1))
data$vT = as.vector(t(G1T[50:1,]))
data$vF = as.vector(t(G1F[50:1,]))
fg1 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "Cluster 1") +
      theme(axis.text.x = element_text(angle = 0, vjust = 0.75, hjust = 0.5)) + 
      theme(axis.text.y = element_text(angle = 0, hjust = 0.75, vjust = 0.5)) + 
      scale_fill_continuous(high = "black", low = "white") + 
      scale_x_discrete("Cluster 1", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
      scale_y_discrete("", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
      theme(axis.text=element_text(size=10, face = "bold"),
            axis.title=element_text(size=12, face="bold"), 
            legend.position="none",
            axis.title.x = element_text(margin = margin(c(7,0,0,0))),
            aspect.ratio=1, plot.margin=unit(c(0,-3,0,-3.65), "cm")) +
      geom_point(col = "green", aes(size=ifelse(vT, "dot", "no_dot"))) +
      geom_point(col = "red", aes(size=ifelse(vF, "dot", "no_dot"))) +
      scale_size_manual(values=c(dot=0.5, no_dot=NA), guide="none") + 
      theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=0.75))

G2T = G2F = matrix(0,50,50)
for(i in 1:50){
  for(j in 1:50){
    if( (G2[i,j]==1)&(AG2new[i,j]>=0.5) ){ G2T[i,j] = 1 }
    if( (G2[i,j]==1)&(AG2new[i,j]<=0.5) ){ G2F[i,j] = 1 }
  }
}
gg2 = AG2new[50:1,]
data$freq = as.vector(t(gg2))
data$vT = as.vector(t(G2T[50:1,]))
data$vF = as.vector(t(G2F[50:1,]))
fg2 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "Cluster 2") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.75, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 0.75, vjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white") + 
  scale_x_discrete("Cluster 2", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
  scale_y_discrete("", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
  theme(axis.text=element_text(size=10, face = "bold"),
        axis.title=element_text(size=12, face="bold"), 
        legend.position="right",
        axis.title.x = element_text(margin = margin(c(7,0,0,0))),
        aspect.ratio=1, plot.margin=unit(c(0,-3,0,-3.65), "cm")) +
  geom_point(col = "green", aes(size=ifelse(vT, "dot", "no_dot"))) +
  geom_point(col = "red", aes(size=ifelse(vF, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=0.5, no_dot=NA), guide="none") +
  theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=0.75))

fg0 = ggarrange(fg1, fg2, ncol = 2, nrow = 1, common.legend = FALSE); fg0



## fitting simulation data with mSSL: The Multivariate Spike-and-Slab LASSO
## uncomment the next two lines to install the R package for mSSL from Github
## library(devtools)
## devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")

## loading libraries
library(SSLASSO)
library(MASS)
library(mSSL)

set.seed(2022)
fit_mSSL_dpe = mSSL_dpe(X,Y)
Omega1 = fit_mSSL_dpe$Omega
OmegaT = Omega1
for(i in 1:50){
  for(j in 1:50){
    OmegaT[i,j] = OmegaT[j,i] = abs(Omega1[i,j]/sqrt(Omega1[i,i]*Omega1[j,j]))
  }
}
GmSSL = matrix(0,50,50)
GmSSL[which(OmegaT>0)] = 1
diag(GmSSL) = 0

## Generating Figure 9 (B) in section 4.3
xaxis = as.character(1:50)
yaxis = xaxis[50:1]
data = expand.grid(X = xaxis, Y = yaxis)

Omega1 = OmegaT[50:1,]
data$value = as.vector(t(Omega1))
data$edge = as.vector(t(GmSSL[50:1,]))
pmSSL = ggplot(data, aes(X, Y, fill = value)) + geom_tile() + labs(y="", x = "mSSL") +
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) + 
      theme(axis.text.y = element_text(angle = 0, hjust = 0.5)) + 
      scale_fill_continuous(high = "black", low = "white", limits = c(0,1)) + 
      scale_x_discrete("(B) mSSL", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
      scale_y_discrete("", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
      theme(axis.text=element_text(size=10, face = "bold"),
            axis.title=element_text(size=12,face="bold"),
            axis.title.x = element_text(margin = margin(c(7,0,0,0))), 
            aspect.ratio=1, plot.margin=unit(c(0,0,0,0), "cm")) +
      geom_point(col = "red", aes(size=ifelse(edge, "dot", "no_dot"))) +
      scale_size_manual(values=c(dot=0.5, no_dot=NA), guide="none") + 
      theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=0.75))



## fitting simulation data with BM13 (Bhadra & Mallick, 2013)
## uncomment the next line to install the R package for BM13 from file "LSGraph_0.1.tar.gz"
## install.packages("LSGraph_0.1.tar.gz", repos = NULL)

## loading libraries
library(LSGraph)
set.seed(2022)
result = LSGraph(Y = Y, X = X, bG = 3.0, DG = diag(50), k = 0,
                 which = 1, m = 1, M = 1, Kmat = diag(2), c = 500, alpha = 0.5, beta = 0.5, 
                 gamma0 = as.vector(rep(0,10)), Adj = diag(50), burnin = 1000, nmc = 2000)
postAdj = result$Adj.mean
GLS = matrix(0,50,50)
GLS[which(postAdj>=0.5)] = 1
diag(GLS) = 0

## Generating Figure 9 (A) in section 4.3
xaxis = as.character(1:50)
yaxis = xaxis[50:1]
data = expand.grid(X = xaxis, Y = yaxis)

postAdj1 = postAdj[50:1,]
data$freq = as.vector(t(postAdj1))
data$edge = as.vector(t(GLS[50:1,]))
pBM13 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "BM13") +
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) + 
      theme(axis.text.y = element_text(angle = 0, hjust = 0.5)) + 
      scale_fill_continuous(high = "black", low = "white") + 
      scale_x_discrete("(A) BM13", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
      scale_y_discrete("", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
      theme(axis.text=element_text(size=10, face = "bold"),
            axis.title=element_text(size=12,face="bold"),
            axis.title.x = element_text(margin = margin(c(7,0,0,0))), 
            aspect.ratio=1, plot.margin=unit(c(0,0,0,0), "cm")) + 
      geom_point(col = "red", aes(size=ifelse(edge, "dot", "no_dot"))) +
      scale_size_manual(values=c(dot=0.5, no_dot=NA), guide="none") + 
      theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=0.75))



## fitting simulation data with JRNS (Samanta et al., 2022)
## three cpp files ("JRNS.cpp", "StepB.cpp", "StepOmega.cpp") can be found in "JRNS.zip"
## loading libraries 
library(gglasso)
library(glasso)
library(glmnet)
library(mvtnorm)
library(rms)
library(dlm)

## compiling JRNS related cpp files
sourceCpp("JRNS.cpp")
sourceCpp("StepB.cpp")
sourceCpp("StepOmega.cpp")

## INITIAL VALUE OF THE PARAMETERS
B_lasso=matrix(0,p,q)
for(i in 1:q){
  cvglm=cv.glmnet(X,Y[,i])
  
  mod_lasso=glmnet(X,Y[,i],alpha=1,family="gaussian",lambda=cvglm$lambda.1se,intercept=F,standardize=F)
  B_lasso[,i]=as.vector(mod_lasso$beta)
}

er=Y-X%*%(B_lasso)
s0=(1/n)*t(er)%*%er
mod_glasso=glasso(s=s0,rho=sqrt(log(p)/n))
sigma_inv_lasso=mod_glasso$wi

fit_JRNS = JRNS(B_ini=B_lasso,Omega_ini=sigma_inv_lasso,Y=Y,X=X,nmc=1000,burnin=2000)

postOmega = fit_JRNS$Omega
S = matrix(0, 50, 50)
for(i in 1:1001){
  temp = postOmega[,,i]
  temp[which(temp!=0)] = 1
  S = S + temp
}
S = S/1001
postAdj = S

GLS = matrix(0,50,50)
GLS[which(postAdj>=0.5)] = 1
diag(GLS) = 0

## Generating Figure 9 (C) in section 4.3
xaxis = as.character(1:50)
yaxis = xaxis[50:1]
data = expand.grid(X = xaxis, Y = yaxis)

postAdj1 = postAdj[50:1,]
data$freq = as.vector(t(postAdj1))
data$edge = as.vector(t(GLS[50:1,]))
pJRNS = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "JRNS") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white") + 
  scale_x_discrete("(A) JRNS", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
  scale_y_discrete("", c("1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50")) +
  theme(axis.text=element_text(size=10, face = "bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.title.x = element_text(margin = margin(c(7,0,0,0))), 
        aspect.ratio=1, plot.margin=unit(c(0,0,0,0), "cm")) + 
  geom_point(col = "red", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=0.5, no_dot=NA), guide="none") + 
  theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=0.75))

pg0 = ggarrange(pBM13, pmSSL, pJRNS, ncol = 3, nrow = 1, common.legend = FALSE); pg0


