# Coverage and convergence of the REML FA results (reml/results/), and a
# check that the testing cells match our partitions (reml/data).
#
#   Rscript reml_check.R

out_dir <- "."
pdir <- "reml/data"
for (nm in c("Wheat_599", "S1_bilinear_J4", "S3_rnorm_J4", "Japonica")) {
  r <- read.csv(file.path(out_dir, "reml/results", paste0("reml_fa_", nm, ".csv")), colClasses = c(line = "character", env = "character"))
  cat("\n=====", nm, ":", nrow(r), "rows\n")
  s <- unique(r[, c("scheme", "partition", "fold", "model", "converged")])
  print(table(s$scheme, s$model))
  cat("partitions:", sort(unique(r$partition)), "\n")
  cat("testing sets not converged:", sum(!s$converged), "of", nrow(s), "\n")
  cat("rows with NA blup:", sum(is.na(r$blup)), "\n")
  # testing cells must match our CV2 / CV1 folds
  pa <- read.csv(file.path(pdir, paste0("partitions_", nm, ".csv")), colClasses = c(line = "character", env = "character"))
  for (sc in c("CV1", "CV2")) {
    x <- r[r$scheme == sc & r$model == "REML FA(1)", ]
    m <- merge(x, pa, by = c("partition", "line", "env"))
    f <- if (sc == "CV1") m$cv1_fold else m$cv2_fold
    cat(sc, ": rows", nrow(x), " matched", nrow(m), " fold agrees", mean(f == m$fold), "\n")
  }
}
