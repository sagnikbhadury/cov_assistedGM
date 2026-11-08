# Draw the breast cancer cluster networks for the README with Network_VisualizationSB(), from the
# results saved by reproduce/run_breast_cancer.R (no MCMC here). Two sets: the networks of the
# authors' full-chain partition (the paper's clusters, graphs refitted with the authors' Graph_NS
# settings) and those of the package's reduced-chain partition.
# Run from the package root:   Rscript reproduce/make_network_figures.R

suppressPackageStartupMessages({library(covassistedGM); library(ggplot2)})
res <- readRDS("inst/extdata/brca_pxg_results.rds")
fig_dir <- "man/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("reproduce/output", showWarnings = FALSE, recursive = TRUE)

draw_set <- function(graphs, title, stem) {
  plots <- plot_pxg_networks(graphs, title = title, thres = 0.5, label_size = 3.2)
  for (k in names(plots))
    ggsave(file.path(fig_dir, sprintf("%s_%s.png", stem, gsub(" ", "", tolower(k)))), plots[[k]],
           width = 6, height = 6.5, dpi = 130, bg = "white")
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    panel <- ggpubr::ggarrange(plotlist = lapply(plots, function(p) p + labs(title = NULL)), ncol = 2, nrow = 2)
    ggsave(file.path(fig_dir, paste0(stem, "_panel.png")), panel, width = 11, height = 12, dpi = 120, bg = "white")
  }
  # check: the edges drawn are exactly the selected edges, with their posterior probabilities
  do.call(rbind, lapply(names(graphs), function(k) {
    g <- graphs[[k]]
    M <- pxg_network_matrix(g$prob, g$adj)
    r <- Network_VisualizationSB(M, NULL, 0.5, k)
    d <- igraph::as_data_frame(r$graph, "edges")
    stopifnot(max(abs(d$Correlation - g$prob[cbind(d$from, d$to)])) < 1e-12)
    data.frame(set = stem, cluster = k, n = g$n, edges_selected = sum(g$adj) / 2, edges_drawn = nrow(d),
               proteins_drawn = igraph::vcount(r$graph), proteins_total = nrow(g$prob))
  }))
}

chk1 <- draw_set(res$cluster_graphs_authors_labels, "PxG network, authors' partition", "brca_network_authors")
chk2 <- draw_set(res$cluster_graphs, sprintf("PxG network, package run (%d + %d iter.)", res$settings$burnin, res$settings$nmcmc),
                 "brca_network_package")
chk <- rbind(chk1, chk2)
print(chk)
write.csv(chk, file.path("reproduce", "output", "network_plot_check.csv"), row.names = FALSE)
cat("network figures written to", fig_dir, "\n")
