# REML factor-analytic models with a genomic relationship matrix, FA(1) and FA(2), fitted with
# ASReml-R 4 on exactly the same partitions as the Bayesian analyses.
#
#   y_ij = mu_j + g_ij + e_ij,   vec(g) ~ N(0, (Lambda Lambda' + Psi) (x) G),   e_ij ~ N(0, s2_j)
#
# For every partition and testing set, the phenotypes of the testing cells are set to missing,
# the model is fitted, and the BLUPs and prediction error variances (PEV) of the testing cells
# are written to results/reml_fa_<data set>.csv. The script can be stopped and restarted:
# testing sets already in the output file are skipped.
#
#   Rscript run_reml_fa.R                         # all data sets
#   Rscript run_reml_fa.R Wheat_599 Japonica      # selected data sets
#   Rscript run_reml_fa.R --partitions=3 Wheat_599   # only partitions 1-3 (CV1, CV2 and sparse)
# Run from the reml folder (it reads data/ and writes results/).

library(asreml)
asreml.options(workspace = "2gb", pworkspace = "2gb", trace = FALSE)

args <- commandArgs(trailingOnly = TRUE)
MAXPART <- Inf
if (any(grepl("^--partitions=", args))) {
  MAXPART <- as.integer(sub("^--partitions=", "", grep("^--partitions=", args, value = TRUE)[1]))
  args <- args[!grepl("^--partitions=", args)]
}
SETS <- if (length(args)) args else c("Wheat_599", "S1_bilinear_J4", "S3_rnorm_J4", "Japonica", "EYT_16_17", "EYT_22_23")
dir.create("results", showWarnings = FALSE)

fit_fa <- function(dat, k, G) {
  fm <- asreml(fixed = y ~ env, random = ~ fa(env, k):vm(line, G),
               residual = ~ dsum(~ units | env), data = dat,
               na.action = na.method(y = "include"), maxit = 50)
  tries <- 0
  while (!isTRUE(fm$converge) && tries < 5) { fm <- update(fm); tries <- tries + 1 }
  pv <- predict(fm, classify = "env:line")$pvals
  list(pv = pv, converged = isTRUE(fm$converge))
}

for (nm in SETS) {
  pheno <- read.csv(file.path("data", paste0("pheno_", nm, ".csv")), colClasses = c("character", "character", "numeric"))
  G <- as.matrix(read.csv(file.path("data", paste0("G_", nm, ".csv")), row.names = 1, check.names = FALSE))
  G <- G + diag(1e-4 * mean(diag(G)), nrow(G))            # keeps G positive definite
  attr(G, "INVERSE") <- FALSE
  pheno$line <- factor(pheno$line, levels = rownames(G))
  pheno$env <- factor(pheno$env)
  cell <- paste(pheno$line, pheno$env)

  # testing sets: decision analysis (CV1, CV2) and, where available, sparse testing (k = 1, 2, 3)
  sets <- list()
  pa <- read.csv(file.path("data", paste0("partitions_", nm, ".csv")), colClasses = c(line = "character", env = "character"))
  for (p in sort(unique(pa$partition[pa$partition <= MAXPART]))) for (f in 1:5) {
    z <- pa[pa$partition == p, ]
    sets[[length(sets) + 1]] <- list(scheme = "CV1", partition = p, fold = f, test = paste(z$line, z$env)[z$cv1_fold == f])
    sets[[length(sets) + 1]] <- list(scheme = "CV2", partition = p, fold = f, test = paste(z$line, z$env)[z$cv2_fold == f])
  }
  sp_file <- file.path("data", paste0("sparse_", nm, ".csv"))
  if (file.exists(sp_file)) {
    sp <- read.csv(sp_file, colClasses = c(line = "character", env = "character"))
    for (p in sort(unique(sp$partition[sp$partition <= MAXPART]))) for (k in sort(unique(sp$k))) {
      z <- sp[sp$partition == p & sp$k == k, ]
      sets[[length(sets) + 1]] <- list(scheme = paste0("sparse_k", k), partition = p, fold = 0, test = paste(z$line, z$env)[!z$train])
    }
  }

  out_file <- file.path("results", paste0("reml_fa_", nm, ".csv"))
  done <- if (file.exists(out_file)) unique(with(read.csv(out_file), paste(scheme, partition, fold, model))) else character(0)
  t0 <- Sys.time()
  for (s in sets) {
    dat <- pheno; dat$y[cell %in% s$test] <- NA
    for (k in 1:2) {
      model <- paste0("REML FA(", k, ")")
      if (paste(s$scheme, s$partition, s$fold, model) %in% done) next
      res <- tryCatch(fit_fa(dat, k, G), error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
      if (is.null(res)) {
        out <- data.frame(data = nm, scheme = s$scheme, partition = s$partition, fold = s$fold, model = model,
                          line = NA, env = NA, blup = NA, pev = NA, converged = FALSE)
      } else {
        pv <- res$pv; key <- paste(pv$line, pv$env); keep <- key %in% s$test
        out <- data.frame(data = nm, scheme = s$scheme, partition = s$partition, fold = s$fold, model = model,
                          line = pv$line[keep], env = pv$env[keep], blup = pv$predicted.value[keep],
                          pev = pv$std.error[keep]^2, converged = res$converged)
      }
      write.table(out, out_file, sep = ",", row.names = FALSE, col.names = !file.exists(out_file), append = file.exists(out_file))
    }
    cat(nm, s$scheme, "partition", s$partition, "fold", s$fold, "-",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
}
