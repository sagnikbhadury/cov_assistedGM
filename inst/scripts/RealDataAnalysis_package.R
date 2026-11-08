### Section 5 "Real Data Analysis: Breast Cancer Data" of Niu et al. (2024), run with the
### covassistedGM package. This is the authors' RealDataAnalysis.R with:
###   library(Rcpp); sourceCpp("PxG-PL.cpp"); load("BRCA_*.RData")  ->  library(covassistedGM); data(...)
### and the post-processing loops replaced by the package helpers that reproduce them.
### Every model call keeps the authors' arguments and seeds. Full run: about 7 minutes.

library(covassistedGM)
data(brca_rna); data(brca_rppa)
RNA <- brca_rna; RPPA <- brca_rppa          # the authors' object names
p <- 3; q <- 15; n <- 883

## PxG model
set.seed(2022)
A <- PxG_PL(YData = RPPA, XData = RNA,
            spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
            mu0 = rep(0, p), sigsq0 = 1, b1 = 30, b2 = 0.1,
            M = 10, N = 4, burnin = 1000, nmcmc = 10000, z0 = sample(1:4, n, TRUE), B = 1)
A$DIC   # 37531.48 in the paper

## graph-only model
set.seed(2022)
C <- GraphOnly_PL(YData = RPPA, XData = RNA,
                  spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
                  mu0 = rep(0, p), sigsq0 = 1, b1 = 30, b2 = 0.1,
                  M = 10, N = 4, burnin = 1000, nmcmc = 10000, z0 = sample(1:4, n, TRUE), B = 1)
C$DIC   # 40798.02 in the paper

## "total sum of squares" of the covariates: the script's quantity is the size-weighted one
within_cluster_ss(RNA, A$Label)$weighted_ss   # 509029.4
within_cluster_ss(RNA, C$Label)$weighted_ss   # 645584.9

## Figures 10 and 11 (ggplot2 and ggpubr as in the script)
if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggpubr", quietly = TRUE)) {
  library(ggplot2)
  cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000")
  cluster_plot <- function(lab, x, y, xl, yl) {
    d <- data.frame(x = RNA[, x], y = RNA[, y], Cluster = paste("Cluster", lab))
    ggplot(d, aes(x = x, y = y, color = Cluster, shape = Cluster)) + geom_point(size = 2.5) +
      scale_color_manual(values = cols) + labs(x = xl, y = yl) +
      theme(axis.title = element_text(size = 13, face = "bold"))
  }
  fig10 <- ggpubr::ggarrange(cluster_plot(A$Label, "ERBB2", "ESR1", "ERBB2 (HER2)", "ESR1 (ER)"),
                             cluster_plot(A$Label, "ERBB2", "PGR", "ERBB2 (HER2)", "PGR (PR)"),
                             cluster_plot(A$Label, "ESR1", "PGR", "ESR1 (ER)", "PGR (PR)"),
                             ncol = 3, nrow = 1, common.legend = TRUE); print(fig10)
  fig11 <- ggpubr::ggarrange(cluster_plot(C$Label, "ERBB2", "ESR1", "ERBB2 (HER2)", "ESR1 (ER)"),
                             cluster_plot(C$Label, "ERBB2", "PGR", "ERBB2 (HER2)", "PGR (PR)"),
                             cluster_plot(C$Label, "ESR1", "PGR", "ESR1 (ER)", "PGR (PR)"),
                             ncol = 3, nrow = 1, common.legend = TRUE); print(fig11)
}

## graph estimate for each PxG cluster (Graph_NS, seed 2022 before each, max-rule symmetrisation,
## cutoff 0.5) = the script's probG1..probG4 / G1..G4
graphs <- cluster_graphs(RPPA, A$Label, spike = 0.1, slab = 5, a1 = 1, a2 = 1, alphaG = 0.5,
                         burnin = 100, nmcmc = 1000, cutoff = 0.5, seed = 2022)

## Figure 12
op <- par(mfrow = c(2, 2))
for (k in seq_along(graphs)) plot_edge_heatmap(graphs[[k]]$prob, cutoff = 0.5, main = names(graphs)[k])
par(op)

## networks drawn with Network_VisualizationSB() (package addition): one ggplot per cluster
plots <- plot_pxg_networks(graphs, title = "PxG cluster networks, TCGA breast cancer", thres = 0.5)
for (k in names(plots)) print(plots[[k]])
if (requireNamespace("ggpubr", quietly = TRUE)) print(ggpubr::ggarrange(plotlist = plots, ncol = 2, nrow = 2))
