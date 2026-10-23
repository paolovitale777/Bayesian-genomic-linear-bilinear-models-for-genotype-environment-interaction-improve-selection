# Exports the data and the exact partitions of the decision (CV1/CV2) and sparse-testing
# analyses, so that the REML FA models can be fitted in ASReml-R on identical folds.
# Reproduces the random draws of gsreg_decision.R (seed 1000 + partition) and
# gsreg_sparse.R (seed 2000 + 10 * partition + k) step by step.
#
#   Rscript gsreg_export_for_reml.R
# Writes to reml/data: pheno_<name>.csv, G_<name>.csv, partitions_<name>.csv,
# sparse_<name>.csv; the true values of the simulations stay in the working folder.

out_dir <- "."
dest <- "reml/data"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
source(file.path(out_dir, "gsreg_functions.R"))

sets <- list(Wheat_599 = 10, EYT_16_17 = 3, EYT_22_23 = 3, Japonica = 10, S1_bilinear_J4 = 10, S3_rnorm_J4 = 10)
sparse_sets <- c("Wheat_599", "S1_bilinear_J4", "S3_rnorm_J4")

for (nm in names(sets)) {
  cs <- load_case(nm); Y <- cs$Y; G <- cs$G
  n <- nrow(Y); J <- ncol(Y); lines <- rownames(Y); envs <- colnames(Y)
  dimnames(G) <- list(lines, lines)
  obs <- which(!is.na(Y), arr.ind = TRUE)
  write.csv(data.frame(line = lines[obs[, 1]], env = envs[obs[, 2]], y = Y[obs]),
            file.path(dest, paste0("pheno_", nm, ".csv")), row.names = FALSE)
  write.csv(G, file.path(dest, paste0("G_", nm, ".csv")))
  if (!is.null(cs$truth)) {
    tv <- data.frame(line = lines[obs[, 1]], env = envs[obs[, 2]], true_value = cs$truth[obs])
    write.csv(tv, file.path(out_dir, paste0("truth_", nm, ".csv")), row.names = FALSE)   # not sent
  }

  # decision analysis: CV1 (lines) and CV2 (cells) folds, partitions 1..R
  parts <- list()
  for (rep in 1:sets[[nm]]) {
    set.seed(1000 + rep)
    fold_line <- sample(rep(1:5, length.out = n))
    fold_cell <- matrix(sample(rep(1:5, length.out = n * J)), n, J)
    parts[[rep]] <- data.frame(partition = rep, line = lines[obs[, 1]], env = envs[obs[, 2]],
                               cv1_fold = fold_line[obs[, 1]], cv2_fold = fold_cell[obs])
  }
  write.csv(do.call(rbind, parts), file.path(dest, paste0("partitions_", nm, ".csv")), row.names = FALSE)

  # sparse testing: environments kept for training, partitions 1..5, k = 1..J-1
  if (nm %in% sparse_sets) {
    sp <- list()
    for (rep in 1:5) for (k in 1:(J - 1)) {
      set.seed(2000 + 10 * rep + k)
      Ytr <- Y
      for (i in 1:n) {
        o <- which(!is.na(Y[i, ]))
        if (length(o) > k) Ytr[i, o[-sample.int(length(o), k)]] <- NA
      }
      sp[[length(sp) + 1]] <- data.frame(partition = rep, k = k, line = lines[obs[, 1]], env = envs[obs[, 2]],
                                         train = !is.na(Ytr[obs]))
    }
    write.csv(do.call(rbind, sp), file.path(dest, paste0("sparse_", nm, ".csv")), row.names = FALSE)
  }
  cat(nm, ": ", n, " lines, ", J, " environments, ", nrow(obs), " records\n", sep = "")
}
