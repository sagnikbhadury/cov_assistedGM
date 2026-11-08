# Provenance of the C++ sources

All files in this directory come from the JASA reproducibility supplement of Niu, Ni, Pati and
Mallick (2024), "Covariate-Assisted Bayesian Graph Learning for Heterogeneous Data"
(doi:10.1080/01621459.2023.2233744). The originals are shipped verbatim under
`inst/original_supplement/Code/` (`PxG-PL.cpp`, `LSGraph_0.1.tar.gz`, `JRNS.zip`).

| Package file | Original | Method |
|---|---|---|
| `PxG_PL.cpp` | `Code/PxG-PL.cpp` | PxG, pseudo-likelihood specification (paper Section 3.2.2) |
| `lsgraph_*.cpp`, `lsgraph_*.h` | `LSGraph_0.1.tar.gz` (`src/`) | BM13 competitor (Bhadra and Mallick 2013; Niu, Pati and Mallick 2020) |
| `jrns_*.cpp` | `JRNS.zip` | JRNS competitor (Samanta, Khare and Michailidis 2022) |

Changes made for packaging, and nothing else:

1. `// [[Rcpp::plugins("cpp11")]]` removed and C++14 selected in `Makevars`; current
   RcppArmadillo refuses to compile in C++11 mode, which is why the original `sourceCpp()` calls
   fail on an up-to-date R installation.
2. `#include <R.h>` and `#include <Rmath.h>` placed before `RcppArmadillo.h` were removed; Rcpp
   includes them itself and warns when they come first.
3. `std::cout` replaced by `Rcpp::Rcout` and placed behind a trailing `verbose` argument (default
   `TRUE`, so the default output matches the original) in `PxG_PL`, `CovOnly_PL`, `GraphOnly_PL`,
   `Graph_NS` and `LSGraph`. R packages may not write to `std::cout`.
4. In `lsgraph_utility.cpp`, `makesplinematrix()` returned the integer `0` from a function whose
   return type is `arma::mat` after printing "error in Kmat"; both branches now call
   `Rcpp::stop("error in Kmat")`. This path is only reached with `which = 3` and a misshaped
   `Kmat`, which no script in the supplement uses.
5. `Rcpp::export` tags inside LSGraph header files were dropped (they were ignored by Rcpp anyway).

Every arithmetic statement, prior, update and random-number call is as in the originals, so
results obtained with the same seeds are the same as with the supplement's `sourceCpp()`
workflow.

Memory note (no change made): the three `*_PL` samplers store an `n x n x nmcmc` co-clustering
cube of 64-bit integers for Dahl's (2006) least-squares partition. For the paper's breast cancer
chain (n = 883, 10,000 kept iterations) this is 62 GB; the authors ran it on a 32 GB Mac with
memory compression. On machines with less memory, use a shorter chain (about
`n^2 * nmcmc * 8` bytes are needed) and rely on the authors' saved results, shipped in
`inst/original_supplement/RData`, for the full-length numbers. The samplers are deliberately
left exactly as the authors wrote them.
