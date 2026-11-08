#' covassistedGM: Covariate-Assisted Bayesian Graph Learning for Heterogeneous Data
#'
#' An R package built from the JASA reproducibility supplement of Niu, Ni, Pati and
#' Mallick (2024). The PxG ("partition times graph") model partitions samples with a
#' truncated Dirichlet process whose allocation is informed jointly by a response matrix
#' \code{Y} (the variables whose conditional-independence graph is of interest) and a
#' covariate matrix \code{X}, and fits a sparse Gaussian graphical model within each
#' cluster. Two specifications from the paper are provided exactly as the authors wrote
#' them:
#' \describe{
#'   \item{Gaussian / G-Wishart (R)}{\code{\link{PxG_GW}}, \code{\link{CovOnly_GW}},
#'     \code{\link{GraphOnly_GW}}, \code{\link{GraphEst}}; paper Section 3.2.1, used in
#'     simulation Examples 1 and 2 (3 and 5 nodes).}
#'   \item{Pseudo-likelihood spike-and-slab (C++)}{\code{\link{PxG_PL}},
#'     \code{\link{CovOnly_PL}}, \code{\link{GraphOnly_PL}}, \code{\link{Graph_NS}}; paper
#'     Section 3.2.2, used in simulation Example 3 (50 nodes) and the breast cancer analysis
#'     (15 proteins, 883 tumours).}
#' }
#' The package adds the TCGA breast cancer data (\code{\link{brca_rna}},
#' \code{\link{brca_rppa}}), data generators for the three simulation examples
#' (\code{\link{simulate_example1}}, \code{\link{simulate_example2}},
#' \code{\link{simulate_example3}}), post-processing helpers
#' (\code{\link{symmetrize_graph}}, \code{\link{cluster_edge_probabilities}},
#' \code{\link{cluster_graphs}}, \code{\link{plot_cluster_network}},
#' \code{\link{plot_edge_heatmap}}, \code{\link{within_cluster_ss}}), and ships the
#' complete original supplement under \code{system.file("original_supplement", package =
#' "covassistedGM")}.
#'
#' @references Niu, Y., Ni, Y., Pati, D. and Mallick, B. K. (2024). Covariate-Assisted Bayesian
#'   Graph Learning for Heterogeneous Data. \emph{Journal of the American Statistical Association},
#'   119(547), 1985-1999. \doi{10.1080/01621459.2023.2233744}. arXiv:2308.07806.
#' @useDynLib covassistedGM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm runif rbinom rbeta var lowess
#' @importFrom graphics image axis par plot points text legend title
#' @importFrom grDevices gray colorRampPalette
#' @importFrom utils head
#' @keywords internal
"_PACKAGE"

#' PxG model, pseudo-likelihood specification (C++)
#'
#' Blocked Gibbs sampler for the PxG model of Niu et al. (2024), Section 3.2.2. Within each
#' cluster the graph is learned by Bayesian neighbourhood selection: every node is regressed on
#' the others with a spike-and-slab prior on the regression coefficients, so the sampler scales
#' to the 50-node simulation and the 15-protein breast cancer analysis of the paper. The
#' covariates enter through the same Gaussian auxiliary model as in \code{\link{PxG_GW}}.
#'
#' This is the authors' \code{PxG_PL} from \code{PxG-PL.cpp}, compiled into the package. The
#' only changes are the replacement of \code{std::cout} by \code{Rcpp::Rcout} behind the added
#' \code{verbose} argument and C++14 compilation; see the README for details.
#'
#' @param YData Response (graphical) data, an \code{n x q} numeric matrix.
#' @param XData Covariate data, an \code{n x p} numeric matrix.
#' @param spike,slab Variances of the spike and slab components of the prior on the
#'   neighbourhood-regression coefficients (scaled by the node variance).
#' @param a1,a2 Shape and rate of the inverse-gamma prior on the node (residual) variances.
#' @param alphaG Prior edge inclusion probability.
#' @param mu0 Prior mean (length \code{p}) of the cluster means in the covariate model.
#' @param sigsq0 Prior variance scale of the cluster means in the covariate model.
#' @param b1,b2 Shape and rate of the inverse-gamma prior on the covariate variance.
#' @param M Concentration parameter of the truncated stick-breaking prior.
#' @param N Maximum number of clusters (truncation level).
#' @param burnin Number of burn-in iterations.
#' @param nmcmc Number of MCMC iterations kept after burn-in.
#' @param z0 Initial cluster allocation, an integer vector of length \code{n} with values in \code{1:N}.
#' @param B Number of neighbourhood-regression sweeps per iteration.
#' @param verbose Logical; print per-iteration progress as the original code did. Default \code{TRUE}.
#' @return A list with components
#' \describe{
#'   \item{Label}{posterior cluster labels chosen by Dahl's (2006) least-squares criterion}
#'   \item{DIC}{deviance information criterion}
#'   \item{logTotLike}{log total likelihood at the selected MCMC sample}
#'   \item{Var}{variance of the log total likelihood over the kept samples}
#'   \item{Gsum}{integer array of dimension \code{(q-1) x q x n}: for every sample, the number of
#'     kept iterations in which each neighbourhood-regression indicator of its cluster's graph
#'     was on. Divide by \code{nmcmc} for probabilities, and see
#'     \code{\link{cluster_edge_probabilities}} to turn it into \code{q x q} matrices.}
#' }
#' @references Niu, Y., Ni, Y., Pati, D. and Mallick, B. K. (2024). Covariate-Assisted Bayesian
#'   Graph Learning for Heterogeneous Data. \emph{Journal of the American Statistical Association},
#'   119(547), 1985-1999. \doi{10.1080/01621459.2023.2233744}
#' @seealso \code{\link{CovOnly_PL}}, \code{\link{GraphOnly_PL}}, \code{\link{Graph_NS}},
#'   \code{\link{cluster_graphs}}
#' @examples
#' \donttest{
#' data(brca_rna); data(brca_rppa)
#' ## the paper uses burnin = 1000, nmcmc = 10000 (about 3 minutes); this is a short run
#' set.seed(2022)
#' fit <- PxG_PL(YData = brca_rppa, XData = brca_rna,
#'               spike = 0.01, slab = 1, a1 = 1, a2 = 1, alphaG = 0.5,
#'               mu0 = rep(0, 3), sigsq0 = 1, b1 = 30, b2 = 0.1,
#'               M = 10, N = 4, burnin = 50, nmcmc = 200,
#'               z0 = sample(1:4, 883, TRUE), B = 1, verbose = FALSE)
#' table(fit$Label)
#' }
#' @name PxG_PL
#' @export
NULL

#' Covariate-only model, pseudo-likelihood specification (C++)
#'
#' Comparison model of Niu et al. (2024): the partition is driven by the covariates alone and
#' one neighbourhood-selection graph is fitted to all of \code{YData} for the DIC.
#' Authors' code; arguments as in \code{\link{PxG_PL}} (there is no \code{B}).
#'
#' @inheritParams PxG_PL
#' @return A list with \code{Label}, \code{DIC}, \code{logTotLike} and \code{Var}.
#' @seealso \code{\link{PxG_PL}}, \code{\link{GraphOnly_PL}}
#' @name CovOnly_PL
#' @export
NULL

#' Graph-only model, pseudo-likelihood specification (C++)
#'
#' Comparison model of Niu et al. (2024): the partition is driven by the graphical data
#' alone; the covariates enter the DIC through a single-cluster Gaussian fit.
#' Authors' code; arguments as in \code{\link{PxG_PL}}.
#'
#' @inheritParams PxG_PL
#' @return A list with \code{Label}, \code{DIC}, \code{logTotLike} and \code{Var}.
#' @seealso \code{\link{PxG_PL}}, \code{\link{CovOnly_PL}}
#' @name GraphOnly_PL
#' @export
NULL

#' Bayesian neighbourhood selection for one graph (C++)
#'
#' Spike-and-slab neighbourhood-selection sampler for a single Gaussian graphical model, used
#' by the authors to estimate the graph of each PxG cluster after the partition has been
#' chosen (paper Figures 8 and 12). The returned matrix is not symmetric: entry \code{(i, j)}
#' is the posterior frequency with which node \code{j} is selected in the regression of node
#' \code{i}. The authors symmetrise it with the maximum rule, see \code{\link{symmetrize_graph}}.
#'
#' Authors' code. Note that the sampler runs \code{burnin + nmcmc - 1} iterations and divides the
#' counts by \code{nmcmc}, so the frequencies are bounded by \code{(nmcmc - 1) / nmcmc}; this is
#' kept as in the original.
#'
#' @inheritParams PxG_PL
#' @return A \code{q x q} matrix of posterior neighbourhood-selection frequencies.
#' @seealso \code{\link{cluster_graphs}}, \code{\link{symmetrize_graph}}
#' @examples
#' \donttest{
#' data(brca_rppa)
#' set.seed(2022)
#' G <- Graph_NS(YData = brca_rppa, spike = 0.1, slab = 5, a1 = 1, a2 = 1, alphaG = 0.5,
#'               burnin = 100, nmcmc = 1000, verbose = FALSE)
#' round(symmetrize_graph(G), 2)
#' }
#' @name Graph_NS
#' @export
NULL
