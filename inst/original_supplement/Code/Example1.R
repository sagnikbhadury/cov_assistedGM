### This R script reproduces section 4.1 "Example 1: piecewise linear precision matrix"
### including Figure 3: posterior edge inclusion probability
### and Figure 4: posterior estimates of partial correlations

rm(list = ls())

## loading libraries
library(huge)
library(invgamma)
library(BDgraph)

## reading functions for the simulation study
## "PxG-GW.R" contains all functions for Gaussian-G-Wishart model (section 3.2.1)
source("PxG-GW.R")

## generating simulation data from the 3-dimensional GGM in section 4.1
## first, generating covariate x
set.seed(2020)
x = rep(0, 300)
x[1:100] = runif(100, -1, -1/3)
x[101:200] = runif(100, -1/3, 1/3)
x[201:300] = runif(100, 1/3, 1)
x = x[order(x)]

## true precision matrix as a piecewise function of x
Omega = function(i, x){
  Omega = diag(3)*1.2
  if(i == 1){
    Omega[1,2] = Omega[2,1] = 0
    Omega[1,3] = Omega[3,1] = -0.75*x+0.25
    Omega[2,3] = Omega[3,2] = 0.75*x+1.25
  }
  if(i == 2){
    Omega[1,2] = Omega[2,1] = 0.75*x+0.75
    Omega[1,3] = Omega[3,1] = 0
    Omega[2,3] = Omega[3,2] = -0.75*x+0.75
  }
  if(i == 3){
    Omega[1,2] = Omega[2,1] = -0.75*x+1.25
    Omega[1,3] = Omega[3,1] = 0.75*x+0.25
    Omega[2,3] = Omega[3,2] = 0   
  }
  return(Omega)
}

## true graph as a function of x
G0 = function(i, x){
  Omega = diag(3)
  if(i == 1){
    Omega[1,2] = Omega[2,1] = 0
    Omega[1,3] = Omega[3,1] = 1
    Omega[2,3] = Omega[3,2] = 1
  }
  if(i == 2){
    Omega[1,2] = Omega[2,1] = 1
    Omega[1,3] = Omega[3,1] = 0
    Omega[2,3] = Omega[3,2] = 1
  }
  if(i == 3){
    Omega[1,2] = Omega[2,1] = 1
    Omega[1,3] = Omega[3,1] = 1
    Omega[2,3] = Omega[3,2] = 0   
  }
  return(Omega)
}

## generating response y according to the true precision matrix
y = matrix(0, 300, 3)
for(i in 1:100){
  temp = solve(Omega(1, x[i]))
  y[i,] = t(t(chol(temp))%*%rnorm(3,0,1)) 
}
for(i in 101:200){
  temp = solve(Omega(2, x[i]))
  y[i,] = t(t(chol(temp))%*%rnorm(3,0,1)) 
}
for(i in 201:300){
  temp = solve(Omega(3, x[i]))
  y[i,] = t(t(chol(temp))%*%rnorm(3,0,1)) 
}

x = matrix(x, 300, 1)
Y = apply(y, 2, function(x) (x-mean(x))/sd(x))
X = apply(x, 2, function(x) (x-mean(x))/sd(x))



## fitting simulation data with the proposed PxG model
set.seed(2020)
A = PxG_GW(Y = Y, X = X, b = 3, D = diag(3), alphaG = 0.1,
           mu0 = as.vector(rep(0,1)), sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5, 
           burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 1)
A$DIC # 2193.304
save(A, file="Example1A.RData")

## fitting simulation data with the covariate-only model
set.seed(2020)
B = CovOnly_GW(Y = Y, X = X, b = 3, D = diag(3), alphaG = 0.1,
               mu0 = as.vector(rep(0,1)), sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5, 
               burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)))
B$DIC # 2584.573
save(B, file="Example1B.RData")

## fitting simulation data with the graph-only model
set.seed(2020)
C = GraphOnly_GW(Y = Y, X = X, b = 3, D = diag(3), alphaG = 0.1,
                 mu0 = as.vector(rep(0,1)), sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5, 
                 burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 1)
C$DIC # 3189.946
save(C, file="Example1C.RData")



## k-lasso: kernel lasso function
ekernel = function(y,x,h){
  return(dnorm(sqrt(sum((y-x)^2))/h))
} 
kernel_mc = function(y,x,e=x,h=1){
  n = nrow(y)
  m = nrow(e)
  p = ncol(y)
  q = ncol(x)
  mu = matrix(0,m,p)
  Sig = array(0,c(m,p,p))
  for (i in 1:m){
    tmp = 0 
    for (j in 1:n){
      mu[i,] =  mu[i,] + ekernel(e[i,],x[j,],h)*y[j,]
      tmp = tmp + ekernel(e[i,],x[j,],h)
    }
    mu[i,] = mu[i,]/tmp
    for (j in 1:n){
      Sig[i,,] = Sig[i,,] + ekernel(e[i,],x[j,],h)*(y[j,]-mu[i,])%*%t(y[j,]-mu[i,])
    }
    Sig[i,,] = Sig[i,,]/tmp
  }
  return(list(mu = mu,Sig = Sig))
}

## fitting simulation data with k-lasso
set.seed(2020)
nh=20
nl=20
h_grid = seq(0.1,1,length=nh)
l_grid = seq(0.1,1,length=nl)
Omega_est=array(0,c(3,3,300))
aic= matrix(0,nh,nl)
time=proc.time()
for (i in 1:300){ print(i)
  for (j in 1:nh){
    r=kernel_mc(Y,X,h=h_grid[j])
    gl = huge(r$Sig[i,,],lambda=l_grid,method="glasso",verbose=FALSE)
    aic[j,] = -2*gl$loglik+2*gl$df
  }
  lmid = arrayInd(which.min(aic),c(nh,nl))
  r=kernel_mc(Y,X,h=h_grid[lmid[1]])
  gl = huge(r$Sig[i,,],lambda=l_grid[lmid[2]],method="glasso",verbose=FALSE)
  Omega_est[,,i] = gl$icov[[1]]
}
proc.time() - time
G_est = (Omega_est!=0)
save(Omega_est, file = "klasso.RData")


## fitting k-lasso takes about 100 minutes to run
## to save time, uncomment the next two lines to load results from PxG and k-lasso
## load("Example1A.RData")
## load("klasso.RData")


## calculate partial correlations from k-lasso simulation results for Figure 4
GK = Omega_est
for(i in 1:300){
  GK[1,2,i] = GK[2,1,i] = Omega_est[1,2,i]/sqrt(Omega_est[1,1,i]*Omega_est[2,2,i])
  GK[1,3,i] = GK[3,1,i] = Omega_est[1,3,i]/sqrt(Omega_est[1,1,i]*Omega_est[3,3,i])
  GK[2,3,i] = GK[3,2,i] = Omega_est[2,3,i]/sqrt(Omega_est[2,2,i]*Omega_est[3,3,i])
}

## generating true edge and partial correlation sequences for plot purposes
xseq = seq(-1, 1, by = 0.01)
Gseq = Oseq = array(0, c(3,3,length(xseq)))
for(i in 1:length(xseq)){
  if(xseq[i]<(-1/3)) { Oseq[,,i] = Omega(1, xseq[i]); Gseq[,,i] = G0(1, xseq[i]) }
  if(xseq[i]>1/3)    { Oseq[,,i] = Omega(3, xseq[i]); Gseq[,,i] = G0(3, xseq[i]) }
  if((xseq[i]>-1/3)&&(xseq[i]<1/3)) 
  { Oseq[,,i] = Omega(2, xseq[i]); Gseq[,,i] = G0(2, xseq[i]) }
}

## true partial correlation sequences
w12 = Oseq[1,2,]/1.2
w13 = Oseq[1,3,]/1.2
w23 = Oseq[2,3,]/1.2

## true edge sequences
e12 = Gseq[1,2,]
e13 = Gseq[1,3,]
e23 = Gseq[2,3,]

## calculate true partial correlations for each covariate sample
xw12 = xw13 = xw23 = rep(0, 300)
for(i in 1:100){
  temp = Omega(1, x[i])
  xw12[i] = temp[1,2]/1.2
  xw13[i] = temp[1,3]/1.2
  xw23[i] = temp[2,3]/1.2
}
for(i in 101:200){
  temp = Omega(2, x[i])
  xw12[i] = temp[1,2]/1.2
  xw13[i] = temp[1,3]/1.2
  xw23[i] = temp[2,3]/1.2
}
for(i in 201:300){
  temp = Omega(3, x[i])
  xw12[i] = temp[1,2]/1.2
  xw13[i] = temp[1,3]/1.2
  xw23[i] = temp[2,3]/1.2
}



## the following codes are for generating Figure 3 and 4 in section 4.1
## loading libraries for plot
library(ggplot2)
library(ggpubr)

## plot Figure 3
a = A$Gsum[1,2,]
temp1 = data.frame(x1 = xseq, y1 = e12)
post = lowess(x, a, f = 0.3)
temp2 = data.frame(x2 = post$x, y2=post$y)
p1 = ggplot() + 
  geom_line(temp1, mapping = aes(x=x1,y=y1, color = "true edge inclusion probability"), size = 0.6) + 
  geom_line(temp2, mapping = aes(x=x2,y=y2, color = "PxG posterior edge inclusion probability"), size = 0.6) +
  scale_colour_manual(values = c("true edge inclusion probability" = "black", "PxG posterior edge inclusion probability" = "red")) +
  xlab("x") + ylab(expression(p[12])) + labs(color = "") +
  theme(axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14))

a = A$Gsum[1,3,]
temp1 = data.frame(x1=xseq, y1=e13)
post = lowess(x, a, f=0.3)
temp2 = data.frame(x2=post$x, y2=post$y)
p2 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true edge inclusion probability"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior edge inclusion probability"), size = 0.6) +
  scale_colour_manual(values = c("true edge inclusion probability" = "black", "PxG posterior edge inclusion probability" = "red")) +
  xlab("x") + ylab(expression(p[13])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14))

a = A$Gsum[2,3,]
temp1 = data.frame(x1=xseq, y1=e23)
post = lowess(x, a, f=0.3)
temp2 = data.frame(x2=post$x, y2=post$y)
p3 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true edge inclusion probability"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = " PxG posterior edge inclusion probability"), size = 0.6) +
  scale_colour_manual(values = c("true edge inclusion probability" = "black", " PxG posterior edge inclusion probability" = "red")) +
  xlab("x") + ylab(expression(p[23])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14))

p0 = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE); p0



## plot Figure 4 and calculate MSE
a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[1,2,i]/sqrt(A$Osum[1,1,i]*A$Osum[2,2,i])
# MSE 
mean((a - xw12)^2) # 0.0251
mean((GK[1,2,]-xw12)^2) # 0.1433
temp1 = data.frame(x1=xseq, y1=w12)
post = lowess(x, a, f=0.3)
temp2 = data.frame(x2=post$x, y2=post$y)
klasso = lowess(x, GK[1,2,], f=0.3)
temp3 = data.frame(x3=klasso$x, y3=klasso$y)
p1 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  geom_line(temp3, mapping=aes(x=x3,y=y3, color = "k-lasso estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red", "k-lasso estimate" = "blue")) +
  xlab("x") + ylab(expression(rho[12])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = -0.65, y = 0.8, label = "MSE") +
  annotate("text", x = -0.65, y = 0.7, label = "PxG  0.0251") +
  annotate("text", x = -0.65, y = 0.6, label = "k-lasso  0.1433")

a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[1,3,i]/sqrt(A$Osum[1,1,i]*A$Osum[3,3,i])
# MSE 
mean((a - xw13)^2) # 0.0221
mean((GK[1,3,]-xw13)^2) # 0.1094
temp1 = data.frame(x1=xseq, y1=w13)
post = lowess(x, a, f=0.3)
temp2 = data.frame(x2=post$x, y2=post$y)
klasso = lowess(x, GK[1,3,], f=0.3)
temp3 = data.frame(x3=klasso$x, y3=klasso$y)
p2 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  geom_line(temp3, mapping=aes(x=x3,y=y3, color = "k-lasso estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red", "k-lasso estimate" = "blue")) +
  xlab("x") + ylab(expression(rho[13])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = 0, y = 0.855, label = "MSE") +
  annotate("text", x = 0, y = 0.755, label = "PxG  0.0221") +
  annotate("text", x = 0, y = 0.655, label = "k-lasso  0.1094")

a = rep(0, length(x))
for(i in 1:length(x)) a[i] = A$Osum[2,3,i]/sqrt(A$Osum[2,2,i]*A$Osum[3,3,i])
# MSE 
mean((a - xw23)^2) # 0.0104
mean((GK[1,2,]-xw23)^2) # 0.2266
temp1 = data.frame(x1=xseq, y1=w23)
post = lowess(x, a, f=0.3)
temp2 = data.frame(x2=post$x, y2=post$y)
klasso = lowess(x, GK[2,3,], f=0.3)
temp3 = data.frame(x3=klasso$x, y3=klasso$y)
p3 = ggplot() + 
  geom_line(temp1, mapping=aes(x=x1,y=y1, color = "true partial correlation"), size = 0.6) + 
  geom_line(temp2, mapping=aes(x=x2,y=y2, color = "PxG posterior estimate"), size = 0.6) +
  geom_line(temp3, mapping=aes(x=x3,y=y3, color = "k-lasso estimate"), size = 0.6) +
  scale_colour_manual(values = c("true partial correlation" = "black", "PxG posterior estimate" = "red", "k-lasso estimate" = "blue")) +
  xlab("x") + ylab(expression(rho[23])) + labs(color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(size=14)) +
  annotate("text", x = 0.6, y = 0.8, label = "MSE") +
  annotate("text", x = 0.6, y = 0.7, label = "PxG  0.0104") +
  annotate("text", x = 0.6, y = 0.6, label = "k-lasso  0.2266")

p0 = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE); p0


