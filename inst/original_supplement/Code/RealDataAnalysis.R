### This R script reproduces section 5 "Real Data Analysis: Breast Cancer Data"
### including Figure 10: posterior clustering results for the PxG model
### Figure 11: posterior clustering results for the graph-only model
### Figure 12: posterior edge inclusion probability for each cluster under the PxG model

rm(list = ls())

## loading libraries and compiling all C++ functions
library(Rcpp)
library(RcppArmadillo)
sourceCpp("PxG-PL.cpp")

## loading TCGA breast cancer data
## mRNA: "BRCA_RNA.RData"
## protein: "BRCA_RPPA.RData"
load("BRCA_RNA.RData")
load("BRCA_RPPA.RData")

## fitting the breast cancer data with the proposed PxG model
p = 3; q = 15; n = 883
set.seed(2022)
A = PxG_PL(YData = RPPA, XData = RNA, 
           spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5, 
           mu0 = rep(0, p), sigsq0 = 1, b1 = 30, b2 = 0.1, 
           M = 10, N = 4, burnin = 1000, nmcmc = 10000, z0 = sample(1:4,n,TRUE), B = 1)
save(A, file = "RealDataAnalysisA.RData")

## fitting the breast cancer data with the graph-only model
set.seed(2022)
C = GraphOnly_PL(YData = RPPA, XData = RNA, 
                 spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5, 
                 mu0 = rep(0, p), sigsq0 = 1, b1 = 30, b2 = 0.1, 
                 M = 10, N = 4, burnin = 1000, nmcmc = 10000, z0 = sample(1:4,n,TRUE), B = 1)
save(C, file = "RealDataAnalysisC.RData")



## DIC for two models
A$DIC # 37531.48
C$DIC # 40798.02

## calculating the total sum of squares of covariates for two models
A1 = which(A$Label==1)
A2 = which(A$Label==2)
A3 = which(A$Label==3)
A4 = which(A$Label==4)

C1 = which(C$Label==1)
C2 = which(C$Label==2)
C3 = which(C$Label==3)
C4 = which(C$Label==4)

## calculating the total sum of squares for the PxG model
mA1 = apply(RNA[A1,1:3],2,mean)
mA2 = apply(RNA[A2,1:3],2,mean)
mA3 = apply(RNA[A3,1:3],2,mean)
mA4 = apply(RNA[A4,1:3],2,mean)
TSSA = 0
for(i in A1){ TSSA = TSSA + sum((RNA[A1,1:3]-mA1)^2) }
for(i in A2){ TSSA = TSSA + sum((RNA[A2,1:3]-mA2)^2) }
for(i in A3){ TSSA = TSSA + sum((RNA[A3,1:3]-mA3)^2) }
for(i in A4){ TSSA = TSSA + sum((RNA[A4,1:3]-mA4)^2) }
TSSA # 509029.4
TSSA/883 # 576.4773

## calculating the total sum of squares for the graph-only model
mC1 = apply(RNA[C1,1:3],2,mean)
mC2 = apply(RNA[C2,1:3],2,mean)
mC3 = apply(RNA[C3,1:3],2,mean)
mC4 = apply(RNA[C4,1:3],2,mean)
TSSC = 0
for(i in C1){ TSSC = TSSC + sum((RNA[C1,1:3]-mC1)^2) }
for(i in C2){ TSSC = TSSC + sum((RNA[C2,1:3]-mC2)^2) }
for(i in C3){ TSSC = TSSC + sum((RNA[C3,1:3]-mC3)^2) }
for(i in C4){ TSSC = TSSC + sum((RNA[C4,1:3]-mC4)^2) }
TSSC # 645584.9
TSSC/883 # 731.1267



## loading libraries for plot
library(ggplot2)
library(ggpubr)

## plot Figure 10
subX = as.data.frame(RNA[,1:2])
Aind = A$Label
Aind[which(A$Label==1)] = "Cluster 1"
Aind[which(A$Label==2)] = "Cluster 2"
Aind[which(A$Label==3)] = "Cluster 3"
Aind[which(A$Label==4)] = "Cluster 4"
subX$Cluster = Aind
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
fg1 = ggplot(subX, aes(x = ERBB2, y = ESR1, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=10, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = "ERBB2 (HER2)", y = "ESR1 (ER)", caption="")

subX = as.data.frame(RNA[,c(1,3)])
Aind = A$Label
Aind[which(A$Label==1)] = "Cluster 1"
Aind[which(A$Label==2)] = "Cluster 2"
Aind[which(A$Label==3)] = "Cluster 3"
Aind[which(A$Label==4)] = "Cluster 4"
subX$Cluster = Aind
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
fg2 = ggplot(subX, aes(x = ERBB2, y = PGR, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=10, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = "ERBB2 (HER2)", y = "PGR (PR)", caption="")

subX = as.data.frame(RNA[,c(2,3)])
Aind = A$Label
Aind[which(A$Label==1)] = "Cluster 1"
Aind[which(A$Label==2)] = "Cluster 2"
Aind[which(A$Label==3)] = "Cluster 3"
Aind[which(A$Label==4)] = "Cluster 4"
subX$Cluster = Aind
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
fg3 = ggplot(subX, aes(x = ESR1, y = PGR, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=10, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = "ESR1 (ER)", y = "PGR (PR)", caption="")

fg0 = ggarrange(fg1, fg2, fg3, ncol = 3, nrow = 1, common.legend = TRUE); fg0



## plot Figure 11
subX = as.data.frame(RNA[,1:2])
Cind = C$Label
Cind[which(C$Label==1)] = "Cluster 1"
Cind[which(C$Label==2)] = "Cluster 2"
Cind[which(C$Label==3)] = "Cluster 3"
Cind[which(C$Label==4)] = "Cluster 4"
subX$Cluster = Cind
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
fg1 = ggplot(subX, aes(x = ERBB2, y = ESR1, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=10, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = "ERBB2 (HER2)", y = "ESR1 (ER)", caption="")

subX = as.data.frame(RNA[,c(1,3)])
Cind = C$Label
Cind[which(C$Label==1)] = "Cluster 1"
Cind[which(C$Label==2)] = "Cluster 2"
Cind[which(C$Label==3)] = "Cluster 3"
Cind[which(C$Label==4)] = "Cluster 4"
subX$Cluster = Cind
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
fg2 = ggplot(subX, aes(x = ERBB2, y = PGR, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=10, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = "ERBB2 (HER2)", y = "PGR (PR)", caption="")

subX = as.data.frame(RNA[,c(2,3)])
Cind = C$Label
Cind[which(C$Label==1)] = "Cluster 1"
Cind[which(C$Label==2)] = "Cluster 2"
Cind[which(C$Label==3)] = "Cluster 3"
Cind[which(C$Label==4)] = "Cluster 4"
subX$Cluster = Cind
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
fg3 = ggplot(subX, aes(x = ESR1, y = PGR, color = Cluster, shape = Cluster) ) +
  geom_point(size = 2.5) + scale_color_manual(values = cols) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        plot.caption = element_text(color="black", size=10, face="bold", hjust = 0.5, vjust = 3)) + 
  labs(x = "ESR1 (ER)", y = "PGR (PR)", caption="")

fg0 = ggarrange(fg1, fg2, fg3, ncol = 3, nrow = 1, common.legend = TRUE); fg0



## getting graph estimate for each cluster under the PxG model
set.seed(2022)
probG1 = Graph_NS(YData = RPPA[which(A$Label==1),], spike = 0.1, slab = 5, 
                  a1 = 1, a2 = 1, alphaG = 0.5, burnin = 100, nmcmc = 1000)

set.seed(2022)
probG2 = Graph_NS(YData = RPPA[which(A$Label==2),], spike = 0.1, slab = 5, 
                  a1 = 1, a2 = 1, alphaG = 0.5, burnin = 100, nmcmc = 1000)

set.seed(2022)
probG3 = Graph_NS(YData = RPPA[which(A$Label==3),], spike = 0.1, slab = 5, 
                  a1 = 1, a2 = 1, alphaG = 0.5, burnin = 100, nmcmc = 1000)

set.seed(2022)
probG4 = Graph_NS(YData = RPPA[which(A$Label==4),], spike = 0.1, slab = 5, 
                  a1 = 1, a2 = 1, alphaG = 0.5, burnin = 100, nmcmc = 1000)

save(probG1, probG2, probG3, probG4, file = "RealDataAnalysisAG.RData")

## making pseudo-likelihood graph estimates symmetric
probG11 = probG22 = probG33 = probG44 = matrix(0,15,15)
for(i in 1:15){
  for(j in 1:15){
    probG11[i,j] = max(probG1[i,j], probG1[j,i])
    probG22[i,j] = max(probG2[i,j], probG2[j,i])
    probG33[i,j] = max(probG3[i,j], probG3[j,i])
    probG44[i,j] = max(probG4[i,j], probG4[j,i])
  }
}
diag(probG11) = 1
diag(probG22) = 1
diag(probG33) = 1
diag(probG44) = 1

G1 = probG11
G1[which(G1>=0.5)] = 1
G1[which(G1<0.5)] = 0
diag(G1) = 0

G2 = probG22
G2[which(G2>=0.5)] = 1
G2[which(G2<0.5)] = 0
diag(G2) = 0

G3 = probG33
G3[which(G3>=0.5)] = 1
G3[which(G3<0.5)] = 0
diag(G3) = 0

G4 = probG44
G4[which(G4>=0.5)] = 1
G4[which(G4<0.5)] = 0
diag(G4) = 0



## plot Figure 12
xaxis = colnames(RPPA)
yaxis = xaxis[15:1]
data = expand.grid(X=xaxis, Y=yaxis)

gg1 = probG11[15:1,]
data$freq = as.vector(t(gg1))
data$edge = as.vector(t(G1[15:1,]))
fg1 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "Cluster 1") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 45, hjust = 0.75, vjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white") + 
  theme(axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=12, face="bold"), 
        legend.position="right",
        axis.title.x = element_text(margin = margin(c(1,0,0,0)), vjust = 4),
        aspect.ratio=1, plot.margin=unit(c(0,-3,0,-3.65), "cm")) +
  geom_point(col = "#EF6F6A", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=2.5, no_dot=NA), guide="none")

gg2 = probG22[15:1,]
data$freq = as.vector(t(gg2))
data$edge = as.vector(t(G2[15:1,]))
fg2 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "Cluster 2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 45, hjust = 0.75, vjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white") + 
  theme(axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=12, face="bold"), 
        legend.position="right",
        axis.title.x = element_text(margin = margin(c(1,0,0,0)), vjust = 4),
        aspect.ratio=1, plot.margin=unit(c(0,-3,0,-3.65), "cm")) +
  geom_point(col = "#EF6F6A", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=2.5, no_dot=NA), guide="none")

gg3 = probG33[15:1,]
data$freq = as.vector(t(gg3))
data$edge = as.vector(t(G3[15:1,]))
fg3 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "Cluster 3") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 45, hjust = 0.75, vjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white") + 
  theme(axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=12, face="bold"), 
        legend.position="right",
        axis.title.x = element_text(margin = margin(c(1,0,0,0)), vjust = 4),
        aspect.ratio=1, plot.margin=unit(c(0,-3,0,-3.65), "cm")) +
  geom_point(col = "#EF6F6A", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=2.5, no_dot=NA), guide="none")

gg4 = probG44[15:1,]
data$freq = as.vector(t(gg4))
data$edge = as.vector(t(G4[15:1,]))
fg4 = ggplot(data, aes(X, Y, fill = freq)) + geom_tile() + labs(y="", x = "Cluster 4") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 45, hjust = 0.75, vjust = 0.5)) + 
  scale_fill_continuous(high = "black", low = "white") + 
  theme(axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=12, face="bold"), 
        legend.position="right",
        axis.title.x = element_text(margin = margin(c(1,0,0,0)), vjust = 4),
        aspect.ratio=1, plot.margin=unit(c(0,-3,0,-3.65), "cm")) +
  geom_point(col = "#EF6F6A", aes(size=ifelse(edge, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=2.5, no_dot=NA), guide="none")

fg0 = ggarrange(fg1, fg2, fg3, fg4, ncol = 2, nrow = 2, common.legend = TRUE, legend="right"); fg0


