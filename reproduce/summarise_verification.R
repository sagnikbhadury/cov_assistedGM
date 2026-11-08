# Print a Markdown table from reproduce/output/verify_against_supplement.rds for the README.
# Run from the package root after reproduce/verify_against_supplement.R has finished.
v <- readRDS("reproduce/output/verify_against_supplement.rds")
yn <- function(x) if (isTRUE(x)) "yes" else "no"
cat(sprintf("Run on %s, %s, covassistedGM %s.\n\n", v$date, v$R, v$package))
cat("| Check | Package | Authors' saved | Labels identical | DIC identical | Edge summaries | Time |\n")
cat("|---|---|---|---|---|---|---|\n")
cat(sprintf("| Example 3, `PxG_PL` (500 + 1000, n = 500, q = 50) | DIC %.4f | DIC %.4f | %s | %s | `Gsum` identical: %s | %.0f min |\n",
            v$example3$DIC, v$example3$DIC_authors, yn(v$example3$labels_identical), yn(v$example3$DIC_identical),
            yn(v$example3$Gsum_identical), v$example3$seconds / 60))
cat(sprintf("| Example 1, `PxG_GW` (500 + 1000, n = 300, q = 3) | DIC %.4f | DIC %.4f | %s | %s | max abs. diff. `Gsum` %.3g, `Osum` %.3g | %.0f min |\n",
            v$example1$DIC, v$example1$DIC_authors, yn(v$example1$labels_identical), yn(v$example1$DIC_identical),
            v$example1$Gsum_max_abs_diff, v$example1$Osum_max_abs_diff, v$example1$seconds / 60))
cat("\nLabel cross-tables (package rows, authors' columns):\n")
print(v$example3$label_table); print(v$example1$label_table)
