#' TCGA breast cancer mRNA expression of three receptor genes
#'
#' Standardised mRNA expression (column means 0, standard deviations 1) of the oestrogen
#' receptor gene \code{ESR1}, the progesterone receptor gene \code{PGR} and \code{ERBB2}
#' (HER2) for 883 TCGA breast tumours. These are the covariates \code{X} of the real data
#' analysis in Section 5 of Niu et al. (2024): the three receptors define the clinical
#' subtypes of breast cancer, and the paper asks whether the protein network differs
#' between the tumour groups they induce. The matrix is the authors' \code{RNA} object from
#' \code{BRCA_RNA.RData}, unchanged apart from the name.
#'
#' @format A numeric matrix with 883 rows (tumours, in the same order as
#'   \code{\link{brca_rppa}}) and 3 columns: \code{ERBB2}, \code{ESR1}, \code{PGR}.
#' @source The Cancer Genome Atlas, downloaded by the authors with TCGA-Assembler 2
#'   (\url{https://github.com/compgenome365/TCGA-Assembler-2}) and distributed with the
#'   JASA supplement of Niu et al. (2024), \doi{10.1080/01621459.2023.2233744}.
#' @seealso \code{\link{brca_rppa}}
"brca_rna"

#' TCGA breast cancer reverse-phase protein array data for 15 proteins
#'
#' Standardised protein expression (column means 0, standard deviations 1) of 15 proteins
#' of the PI3K/AKT/mTOR signalling axis and the three receptors for the same 883 TCGA breast
#' tumours as \code{\link{brca_rna}}. These are the graphical data \code{Y} of the real data
#' analysis in Section 5 of Niu et al. (2024). The matrix is the authors' \code{RPPA} object
#' from \code{BRCA_RPPA.RData}, unchanged apart from the name.
#'
#' @format A numeric matrix with 883 rows (tumours) and 15 columns: \code{PTEN},
#'   \code{PIK3CA}, \code{BAD}, \code{NFKB1}, \code{TSC1}, \code{IRS1}, \code{EIF4EBP1},
#'   \code{EGFR}, \code{AKT1}, \code{TSC2}, \code{MTOR}, \code{DVL3}, \code{ERBB2},
#'   \code{ESR1}, \code{PGR}.
#' @source The Cancer Genome Atlas, downloaded by the authors with TCGA-Assembler 2
#'   (\url{https://github.com/compgenome365/TCGA-Assembler-2}) and distributed with the
#'   JASA supplement of Niu et al. (2024), \doi{10.1080/01621459.2023.2233744}.
#' @seealso \code{\link{brca_rna}}
"brca_rppa"
