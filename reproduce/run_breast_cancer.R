# Section 5 of Niu et al. (2024) with the packaged (unchanged) samplers on the TCGA breast cancer
# data, and the compact results used by the README and vignette("pxg-breast-cancer").
#
# Memory: the authors' *_PL samplers keep an n x n x nmcmc integer cube, i.e. n^2 * nmcmc * 8
# bytes. The paper's chain (883 tumours, 10,000 kept draws) needs 62 GB, so by default this
# script runs a reduced chain (200 burn-in, 1,000 kept: 6.2 GB) on all 883 tumours, and keeps
# the authors' saved full-chain results next to it. Set PXG_FULL=TRUE for the paper's chain on a
# machine with enough memory.
#
# Run from the package root:   Rscript reproduce/run_breast_cancer.R

suppressPackageStartupMessages(library(covassistedGM))
FULL <- as.logical(Sys.getenv("PXG_FULL", "FALSE"))
out_dir <- "inst/extdata"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- "man/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("reproduce/output", showWarnings = FALSE, recursive = TRUE)

data(brca_rna); data(brca_rppa)
p <- 3; q <- 15; n <- 883
burnin <- if (FULL) 1000 else 200
nmcmc  <- if (FULL) 10000 else 1000
cat(sprintf("chain: burnin %d, nmcmc %d; membership cube %.1f GB\n", burnin, nmcmc, n^2 * nmcmc * 8 / 1e9))

## authors' saved full-chain results (same code, 1000 + 10000 iterations, seed 2022)
sup <- system.file("original_supplement/RData", package = "covassistedGM")
e <- new.env()
load(file.path(sup, "RealDataAnalysisA.RData"), envir = e)
load(file.path(sup, "RealDataAnalysisC.RData"), envir = e)
load(file.path(sup, "RealDataAnalysisAG.RData"), envir = e)

## 1. PxG model: RealDataAnalysis.R call, only the chain length differs when FULL is FALSE
set.seed(2022)
t0 <- proc.time()
A <- PxG_PL(YData = brca_rppa, XData = brca_rna,
            spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
            mu0 = rep(0, p), sigsq0 = 1, b1 = 30, b2 = 0.1,
            M = 10, N = 4, burnin = burnin, nmcmc = nmcmc, z0 = sample(1:4, n, TRUE), B = 1,
            verbose = FALSE)
tA <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("PxG_PL: %.0f s, DIC %.2f, cluster sizes %s\n", tA, A$DIC, paste(table(A$Label), collapse = " ")))
gc()

## 2. Graph-only model
set.seed(2022)
t0 <- proc.time()
C <- GraphOnly_PL(YData = brca_rppa, XData = brca_rna,
                  spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
                  mu0 = rep(0, p), sigsq0 = 1, b1 = 30, b2 = 0.1,
                  M = 10, N = 4, burnin = burnin, nmcmc = nmcmc, z0 = sample(1:4, n, TRUE), B = 1,
                  verbose = FALSE)
tC <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("GraphOnly_PL: %.0f s, DIC %.2f, cluster sizes %s\n", tC, C$DIC, paste(table(C$Label), collapse = " ")))
gc()

## 3. Within-cluster sums of squares of the covariates (the script reports the size-weighted form)
ssA <- within_cluster_ss(brca_rna, A$Label); ssC <- within_cluster_ss(brca_rna, C$Label)
ssA_auth <- within_cluster_ss(brca_rna, e$A$Label); ssC_auth <- within_cluster_ss(brca_rna, e$C$Label)
cat(sprintf("within-cluster SS, package chain: PxG %.1f, graph-only %.1f\n", ssA$ss, ssC$ss))
cat(sprintf("within-cluster SS, authors' labels: PxG %.1f, graph-only %.1f\n", ssA_auth$ss, ssC_auth$ss))
cat(sprintf("script's expression, authors' labels: PxG %.1f, graph-only %.1f (script: 509029.4, 645584.9)\n",
            ssA_auth$script_value, ssC_auth$script_value))

## 4. Agreement between the reduced-chain partition and the authors' full-chain partition
agree <- function(a, b) { tb <- table(a, b); sum(apply(tb, 1, max)) / length(a) }
cat(sprintf("PxG label agreement with authors (best matching): %.3f\n", agree(A$Label, e$A$Label)))
print(table(package = A$Label, authors = e$A$Label))

## 5. One graph per cluster with Graph_NS, the script's settings and seed, for both partitions
t0 <- proc.time()
graphs_pkg  <- cluster_graphs(brca_rppa, A$Label,   spike = 0.1, slab = 5, a1 = 1, a2 = 1, alphaG = 0.5,
                              burnin = 100, nmcmc = 1000, cutoff = 0.5, seed = 2022)
graphs_auth <- cluster_graphs(brca_rppa, e$A$Label, spike = 0.1, slab = 5, a1 = 1, a2 = 1, alphaG = 0.5,
                              burnin = 100, nmcmc = 1000, cutoff = 0.5, seed = 2022)
tG <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Graph_NS x 8: %.0f s; edges per cluster (package labels): %s; (authors' labels): %s\n", tG,
            paste(sapply(graphs_pkg, function(g) sum(g$adj) / 2), collapse = " "),
            paste(sapply(graphs_auth, function(g) sum(g$adj) / 2), collapse = " ")))

## the authors' own Graph_NS output (probG1..4) should match the refit on their labels exactly
auth_prob <- lapply(list(e$probG1, e$probG2, e$probG3, e$probG4), function(G) {
  dimnames(G) <- list(colnames(brca_rppa), colnames(brca_rppa)); symmetrize_graph(G) })
match_auth <- sapply(1:4, function(k) max(abs(graphs_auth[[k]]$prob - auth_prob[[k]])))
cat("max |refit - authors' saved Graph_NS| per cluster:", format(match_auth, digits = 3), "\n")

## 6. Compact results (no MCMC draws)
res <- list(
  settings = list(burnin = burnin, nmcmc = nmcmc, full = FULL, seed = 2022, n = n, q = q,
                  R = R.version.string, package = as.character(packageVersion("covassistedGM")), date = Sys.Date()),
  pxg = list(Label = A$Label, DIC = A$DIC, logTotLike = A$logTotLike, Var = A$Var, seconds = tA),
  graph_only = list(Label = C$Label, DIC = C$DIC, logTotLike = C$logTotLike, Var = C$Var, seconds = tC),
  ss = list(pxg = ssA, graph_only = ssC, pxg_authors = ssA_auth, graph_only_authors = ssC_auth),
  cluster_graphs = lapply(graphs_pkg, function(g) g[c("prob", "adj", "n")]),
  cluster_graphs_authors_labels = lapply(graphs_auth, function(g) g[c("prob", "adj", "n")]),
  gsum_edge_probabilities = cluster_edge_probabilities(A, nmcmc = nmcmc, node_names = colnames(brca_rppa)),
  authors = list(pxg_Label = e$A$Label, pxg_DIC = e$A$DIC, graph_only_Label = e$C$Label, graph_only_DIC = e$C$DIC,
                 cluster_prob = auth_prob, burnin = 1000, nmcmc = 10000),
  agreement = agree(A$Label, e$A$Label), refit_matches_authors = match_auth
)
saveRDS(res, file.path(out_dir, "brca_pxg_results.rds"), compress = "xz")
cat("saved", file.path(out_dir, "brca_pxg_results.rds"), "\n")

## 7. README figures: clusters (Figure 10) for both partitions, edge heat maps (Figure 12)
cols <- c("#1170AA", "#EF6F6A", "#009E73", "#000000"); pchs <- c(16, 17, 15, 18)
cluster_panel <- function(lab, main) {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  plot(brca_rna[, "ERBB2"], brca_rna[, "ESR1"], col = cols[lab], pch = pchs[lab], cex = 0.8, xlab = "ERBB2 (HER2)", ylab = "ESR1 (ER)", main = main)
  plot(brca_rna[, "ERBB2"], brca_rna[, "PGR"],  col = cols[lab], pch = pchs[lab], cex = 0.8, xlab = "ERBB2 (HER2)", ylab = "PGR (PR)")
  plot(brca_rna[, "ESR1"],  brca_rna[, "PGR"],  col = cols[lab], pch = pchs[lab], cex = 0.8, xlab = "ESR1 (ER)", ylab = "PGR (PR)")
  legend("bottomright", paste("Cluster", 1:4), col = cols, pch = pchs, bty = "n", cex = 0.9)
}
png(file.path(fig_dir, "brca_clusters.png"), width = 1500, height = 520, res = 150)
cluster_panel(A$Label, sprintf("PxG, package run (%d + %d iter.)", burnin, nmcmc)); dev.off()
png(file.path(fig_dir, "brca_clusters_authors.png"), width = 1500, height = 520, res = 150)
cluster_panel(e$A$Label, "PxG, authors' chain (1000 + 10000 iter.)"); dev.off()
png(file.path(fig_dir, "brca_edge_heatmaps.png"), width = 1800, height = 1800, res = 150)
par(mfrow = c(2, 2))
for (k in seq_along(graphs_auth)) plot_edge_heatmap(graphs_auth[[k]]$prob, cutoff = 0.5,
                                                    main = sprintf("%s (n = %d), authors' partition", names(graphs_auth)[k], graphs_auth[[k]]$n))
dev.off()
cat("figures written to", fig_dir, "; now run reproduce/make_network_figures.R\n")
