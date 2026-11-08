# covassistedGM 0.1.0 (2026-11-08)

First release. The package is built from the JASA reproducibility supplement of Niu, Ni, Pati
and Mallick (2024), "Covariate-Assisted Bayesian Graph Learning for Heterogeneous Data"
(doi:10.1080/01621459.2023.2233744).

* PxG model and its covariate-only and graph-only comparison models in both specifications of
  the paper: Gaussian / G-Wishart in R (`PxG_GW`, `CovOnly_GW`, `GraphOnly_GW`, `GraphEst`) and
  pseudo-likelihood spike-and-slab in C++ (`PxG_PL`, `CovOnly_PL`, `GraphOnly_PL`, `Graph_NS`).
  The authors' code is unchanged apart from printing (`verbose`), includes and compile flags.
* Comparison methods compiled or written from the same supplement: BM13 (`LSGraph`, `MCMCGraph`,
  `HIWsim`, `meanOmegaHIW`), JRNS (`JRNS`, `StepB`, `StepOmega`, `jrns_initial_values`,
  `jrns_edge_frequencies`), k-lasso (`klasso`, `kernel_mc`, `ekernel`); mSSL wrapped
  (`fit_mSSL`, needs the GitHub package).
* TCGA breast cancer data of the paper (`brca_rna`, `brca_rppa`) and the authors' saved
  intermediate results (`inst/original_supplement/RData`).
* Data generators for the three simulation studies (`simulate_example1/2/3`), post-processing
  helpers (`symmetrize_graph`, `cluster_edge_probabilities`, `cluster_graphs`,
  `partial_correlation_path`, `within_cluster_ss`), heat maps (`plot_edge_heatmap`) and network
  drawing with `Network_VisualizationSB` (plus `pxg_network_matrix`, `plot_pxg_networks`,
  `plot_cluster_network`).
* Four vignettes, package-based versions of the four original scripts (`inst/scripts`), a
  reproduction folder (`reproduce/`) and the complete original supplement
  (`inst/original_supplement`).
