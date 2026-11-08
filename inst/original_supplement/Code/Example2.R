### This R script reproduces section 4.2 "Example 2: linear precision matrix with constant graph"
### including Figure 5: posterior estimates of partial correlations under the PxG model
### and Figure 6: posterior estimates of edge inclusion probability (BM13) or partial correlation (mSSL)

rm(list = ls())

## reading functions for the simulation study
## "PxG-GW.R" contains all functions for Gaussian-G-Wishart model (section 3.2.1)
source("PxG-GW.R")

## generating simulation data from the 5-dimensional GGM in section 4.2
## first, generating covariate x
set.seed(2020)
x = runif(300, -0.8, 0.8)

## true precision matrix as a linear function of x
Omega = function(x){
  G0 = diag(5) - diag(5)
  G0[1,2] = G0[2,3] = G0[3,4] = G0[4,5] = 1
  G0 = G0 + t(G0)

  Omega = x*G0
  diag(Omega) = 1.4
  return(Omega)
}

## true graph
G0 = diag(5) - diag(5)
G0[1,2] = G0[2,3] = G0[3,4] = G0[4,5] = 1
G0 = G0 + t(G0)
diag(G0) = 1

## generating response y according to the true precision matrix
y = matrix(0, 300, 5)
for(i in 1:300){
  temp = solve(Omega(x[i]))
  y[i,] = t(t(chol(temp))%*%rnorm(5,0,1)) 
}
X = matrix(x, 300, 1); Y = y



## fitting simulation data with the proposed PxG model
set.seed(2020)
A = PxG_GW(Y = Y, X = X, b = 3, D = diag(5), alphaG = 0.2,
          mu0 = as.vector(rep(0,1)), sigsq0 = 1, b1 = 10, b2 = 0.1, M = 1, N = 15, 
          burnin = 1000, nmcmc = 2000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 3)
save(A, file="Example2.RData")



## generating true edge and partial correlation sequences for plot purposes
xseq = seq(-0.8, 0.8, by = 0.01)
Gseq = Oseq = array(0, c(5,5,length(xseq)))
for(i in 1:length(xseq)){
  Oseq[,,i] = Omega(xseq[i])
}

## true partial correlation sequences
w12 = Oseq[1,2,]/1.4
w23 = Oseq[2,3,]/1.4
w34 = Oseq[3,4,]/1.4
w45 = Oseq[4,5,]/1.4



## load("Example2.RData")
## the following codes are for generating Figure 5 in section 4.2
## loading libraries for plot
library(ggplot2)
library(ggpubr)
## true partial correlation for each sample of covariate x
truerho = x/1.4

a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[1,2,i]/sqrt(A$Osum[1,1,i]*A$Osum[2,2,i])
## MSE for edge 1-2
MSE12 = mean((a-truerho)^2); MSE12 # 0.0077
temp1 = data.frame(x1=xseq, y1=w12)
post = lowess(x, a, f=0.4)
temp2 = data.frame(x2=post$x, y2=post$y)
p1 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red")) +
  xlab("x") + ylab(expression(rho[12])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = -0.35, y = 0.45, label = "MSE 0.0077")

a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[2,3,i]/sqrt(A$Osum[2,2,i]*A$Osum[3,3,i])
## MSE for edge 2-3
MSE23 = mean((a-truerho)^2); MSE23 # 0.0078
temp1 = data.frame(x1=xseq, y1=w23)
post = lowess(x, a, f=0.4)
temp2 = data.frame(x2=post$x, y2=post$y)
p2 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red")) +
  xlab("x") + ylab(expression(rho[23])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = -0.35, y = 0.45, label = "MSE 0.0078")

a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[3,4,i]/sqrt(A$Osum[4,4,i]*A$Osum[3,3,i])
## MSE for edge 3-4
MSE34 = mean((a-truerho)^2); MSE34 # 0.0078
temp1 = data.frame(x1=xseq, y1=w34)
post = lowess(x, a, f=0.4)
temp2 = data.frame(x2=post$x, y2=post$y)
p3 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red")) +
  xlab("x") + ylab(expression(rho[34])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = -0.35, y = 0.45, label = "MSE 0.0078")

a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[4,5,i]/sqrt(A$Osum[4,4,i]*A$Osum[5,5,i])
## MSE for edge 4-5
MSE45 = mean((a-truerho)^2); MSE45 # 0.0125
temp1 = data.frame(x1=xseq, y1=w45)
post = lowess(x, a, f=0.4)
temp2 = data.frame(x2=post$x, y2=post$y)
p4 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red")) +
  xlab("x") + ylab(expression(rho[45])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = -0.35, y = 0.45, label = "MSE 0.0125")

p0 = ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = TRUE); p0



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
for(i in 1:5){
  for(j in 1:5){
    OmegaT[i,j] = OmegaT[j,i] = Omega1[i,j]/sqrt(Omega1[i,i]*Omega1[j,j])
  }
}
## absolute values of partial correlations
OmegaT = abs(OmegaT)

## Generating Figure 6 (B) in section 4.2
xaxis = as.character(1:5)
yaxis = xaxis[5:1]
data = expand.grid(X = xaxis, Y = yaxis)

diag(G0) = 0
OmegaT = OmegaT[5:1,]
data$value = as.vector(t(OmegaT))
data$edge = as.vector(t(G0[5:1,]))
pmSSL = ggplot(data, aes(X, Y, fill = value)) + geom_tile() + labs(y="", x = "(B) mSSL") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white", limits = c(0,1)) + 
  theme(axis.text=element_text(size=10, face = "bold"),
        axis.title=element_text(size=13,face="bold"), aspect.ratio=1,
        plot.margin=unit(c(0,0,0,0), "cm")) +
  geom_point(col = "#EF6F6A", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=7, no_dot=NA), guide="none"); pmSSL



## fitting simulation data with BM13 (Bhadra & Mallick, 2013)
## uncomment the next line to install the R package for BM13 from file "LSGraph_0.1.tar.gz"
## install.packages("LSGraph_0.1.tar.gz", repos = NULL)

## loading libraries
library(LSGraph)
set.seed(2022)
result = LSGraph(Y = Y, X = X, bG = 3.0, DG = diag(5), k = 0,
                 which = 1, m = 1, M = 1, Kmat = diag(2), c = 300, alpha = 0.5, beta = 0.5, 
                 gamma0 = as.vector(rep(0,1)), Adj = diag(5), burnin = 1000, nmc = 2000)
postAdj = result$Adj.mean

## Generating Figure 6 (A) in section 4.2
postAdj = result$Adj.mean
xaxis = as.character(1:5)
yaxis = xaxis[5:1]
data = expand.grid(X = xaxis, Y = yaxis)

postAdj = postAdj[5:1,]
data$prob = as.vector(t(postAdj))
diag(G0) = 0
data$edge = as.vector(t(G0[5:1,]))
pBM13 = ggplot(data, aes(X, Y, fill = prob)) + geom_tile() + labs(y="", x = "(A) BM13") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white", limits=c(0, 1)) + 
  theme(axis.text=element_text(size=10, face = "bold"),
        axis.title=element_text(size=13,face="bold"), aspect.ratio=1,
        plot.margin=unit(c(0,0,0,0), "cm")) +
  scale_size_manual(values=c(dot=7, no_dot=NA), guide="none")+
  geom_point(col = "#EF6F6A", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=7, no_dot=NA), guide="none"); pBM13

## combining 2 figures
fig = ggarrange(pBM13, pmSSL, nrow = 1, ncol = 2, common.legend = FALSE); fig


