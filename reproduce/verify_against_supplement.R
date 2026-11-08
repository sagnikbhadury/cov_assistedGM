# Verify that the packaged samplers reproduce the authors' saved results bit for bit when run
# with the same seeds. Two checks that finish in reasonable time:
#   (1) Example 3, PxG_PL with the paper's chain (500 + 1000 iterations, n = 500, q = 50);
#   (2) Example 1, PxG_GW with the paper's chain (500 + 1000 iterations, n = 300, q = 3).
# Run from the package root:  Rscript reproduce/verify_against_supplement.R
# Output: reproduce/verify_against_supplement.log (written by the caller) and a small rds.

suppressPackageStartupMessages(library(covassistedGM))
sup <- system.file("original_supplement/RData", package = "covassistedGM")
load_sup <- function(f) { e <- new.env(); load(file.path(sup, f), envir = e); as.list(e) }
out <- list(R = R.version.string, package = as.character(packageVersion("covassistedGM")), date = Sys.Date())

## (1) Example 3, PxG_PL
sim3 <- simulate_example3(seed = 2022)
set.seed(2022)
t0 <- proc.time()
A3 <- PxG_PL(YData = sim3$Y, XData = sim3$X, spike = 0.01, slab = 5, a1 = 1, a2 = 1, alphaG = 0.1,
             mu0 = rep(0, 10), sigsq0 = 1, b1 = 1, b2 = 0.1, M = 1, N = 3,
             burnin = 500, nmcmc = 1000, z0 = sample(1:3, 500, TRUE, c(1/3, 1/3, 1/3)), B = 1, verbose = FALSE)
t3 <- (proc.time() - t0)[["elapsed"]]
S3 <- load_sup("Example3A.RData")$A
out$example3 <- list(seconds = t3, DIC = A3$DIC, DIC_authors = S3$DIC,
                     labels_identical = identical(as.integer(A3$Label), as.integer(S3$Label)),
                     DIC_identical = isTRUE(all.equal(A3$DIC, S3$DIC)),
                     Gsum_identical = identical(as.integer(A3$Gsum), as.integer(S3$Gsum)),
                     label_table = table(A3$Label, S3$Label))
cat(sprintf("Example 3 PxG_PL: %.0f s; DIC %.4f vs authors %.4f; labels identical %s; Gsum identical %s\n",
            t3, A3$DIC, S3$DIC, out$example3$labels_identical, out$example3$Gsum_identical))

## (2) Example 1, PxG_GW
sim1 <- simulate_example1(seed = 2020)
set.seed(2020)
t0 <- proc.time()
A1 <- PxG_GW(Y = sim1$Y, X = sim1$X, b = 3, D = diag(3), alphaG = 0.1,
             mu0 = 0, sigsq0 = 1, b1 = 1, b2 = 1, M = 1, N = 5,
             burnin = 500, nmcmc = 1000, z0 = sample(1:15, 300, TRUE, rep(1/15, 15)), B = 1, verbose = FALSE)
t1 <- (proc.time() - t0)[["elapsed"]]
S1 <- load_sup("Example1A.RData")$A
out$example1 <- list(seconds = t1, DIC = A1$DIC, DIC_authors = S1$DIC,
                     labels_identical = identical(as.integer(A1$Label), as.integer(S1$Label)),
                     DIC_identical = isTRUE(all.equal(A1$DIC, S1$DIC)),
                     Gsum_max_abs_diff = max(abs(A1$Gsum - S1$Gsum)),
                     Osum_max_abs_diff = max(abs(A1$Osum - S1$Osum)),
                     label_table = table(A1$Label, S1$Label))
cat(sprintf("Example 1 PxG_GW: %.0f s; DIC %.4f vs authors %.4f; labels identical %s; max |Gsum diff| %.3g\n",
            t1, A1$DIC, S1$DIC, out$example1$labels_identical, out$example1$Gsum_max_abs_diff))

dir.create("reproduce/output", showWarnings = FALSE, recursive = TRUE)
saveRDS(out, "reproduce/output/verify_against_supplement.rds")
cat("saved reproduce/output/verify_against_supplement.rds\n")
