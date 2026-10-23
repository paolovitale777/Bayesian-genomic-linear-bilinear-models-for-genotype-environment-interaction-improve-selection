# Data scarcity: does the Bayesian genomic FA gain over the REML FA when data are few and
# unbalanced (as in Japonica)? Subsets of n lines are drawn from a data set, and in each subset
#   sparse_k1: every line is phenotyped in only one of the J environments, or
#   unbal_50:  every cell is observed with probability 0.5 (at least one cell per line);
# the models predict the remaining cells. For each environment: selection differential of the
# top 10% of the predicted cells, precision, correlation and Brier score of P(top 10%), against the
# true genetic values in simulations and the phenotypes otherwise.
#
# The design is written first to reml/data/scarcity_<name>.csv, so that the REML FA models
# (run_reml_scarcity.R, ASReml-R) are fitted on exactly the same training sets.
#
#   Rscript gsreg_scarcity.R Wheat_599 10 2000 500          # name, partitions, iterations, burn-in
#   Rscript gsreg_scarcity.R S1_bilinear_J4 10 2000 500 design   # only write the design
# Writes scarcity_<name>.csv (one row per n x scheme x partition x environment x model); restarts skip
# the testing sets already written.

args <- commandArgs(trailingOnly = TRUE)
NAME  <- args[1]
REPS  <- if (length(args) >= 2) as.integer(args[2]) else 10
NITER <- if (length(args) >= 3) as.integer(args[3]) else 2000
BURN  <- if (length(args) >= 4) as.integer(args[4]) else 500
ONLY_DESIGN <- length(args) >= 5 && args[5] == "design"
NSUB  <- c(75, 150, 300)
SCHEMES <- c("sparse_k1", "unbal_50")
TOP   <- 0.10
out_dir <- "."
paolo_dir <- "reml/data"
source(file.path(out_dir, "gsreg_functions.R"))

cs <- load_case(NAME); Y <- cs$Y; G <- cs$G; truth <- cs$truth
V <- if (is.null(truth)) Y else truth
J <- ncol(Y); envs <- colnames(Y)

# ---------------------------------------------------------------- design (same random draws for every model)
design <- list()
for (a in seq_along(NSUB)) for (rep in 1:REPS) {
  set.seed(3000 + 100 * a + rep)
  sub <- sort(sample(nrow(Y), NSUB[a]))
  for (sc in SCHEMES) {
    O <- !is.na(Y[sub, , drop = FALSE])
    train <- matrix(FALSE, length(sub), J)
    for (i in seq_along(sub)) {
      o <- which(O[i, ])
      if (sc == "sparse_k1") train[i, o[sample.int(length(o), 1)]] <- TRUE
      else { t <- o[runif(length(o)) < 0.5]; if (!length(t)) t <- o[sample.int(length(o), 1)]; train[i, t] <- TRUE }
    }
    idx <- which(O, arr.ind = TRUE)
    design[[length(design) + 1]] <- data.frame(n_sub = NSUB[a], scheme = sc, partition = rep,
      line = rownames(Y)[sub][idx[, 1]], env = envs[idx[, 2]], train = train[idx])
  }
}
design <- do.call(rbind, design)
write.csv(design, file.path(paolo_dir, paste0("scarcity_", NAME, ".csv")), row.names = FALSE)
cat(NAME, "- design written:", nrow(design), "rows\n")
if (ONLY_DESIGN) quit(save = "no")

# ---------------------------------------------------------------- fits
out_file <- file.path(out_dir, paste0("scarcity_", NAME, ".csv"))
done <- if (file.exists(out_file)) unique(with(read.csv(out_file), paste(n_sub, scheme, partition))) else character(0)
t0 <- Sys.time()
for (key in unique(paste(design$n_sub, design$scheme, design$partition))) {
  if (key %in% done) next
  z <- design[paste(design$n_sub, design$scheme, design$partition) == key, ]
  sub <- match(unique(z$line), rownames(Y))
  Ys <- Y[sub, , drop = FALSE]; Gs <- G[sub, sub]; Vs <- V[sub, , drop = FALSE]; n <- length(sub)
  Ytr <- matrix(NA, n, J, dimnames = dimnames(Ys))
  tr <- z[z$train, ]; Ytr[cbind(match(tr$line, rownames(Ys)), match(tr$env, envs))] <- Ys[cbind(match(tr$line, rownames(Ys)), match(tr$env, envs))]
  test <- is.na(Ytr) & !is.na(Ys)
  models <- make_fitters(Ys, Gs, niter = NITER, burn = BURN, label = paste0("scar_", NAME))
  out <- list()
  for (m in names(models)) {
    fit <- models[[m]](Ytr)
    A <- if (is.null(truth)) fit$pred else fit$draws        # predictive vs. phenotype; genetic vs. truth
    for (j in 1:J) {
      T <- which(test[, j]); nt <- length(T)
      if (nt < 20) next
      kk <- ceiling(TOP * nt)
      mu <- fit$mean[T, j]; v <- Vs[T, j]
      sel <- rank(-mu, ties.method = "first") <= kk
      top <- rank(-v, ties.method = "first") <= kk
      p_top <- rowMeans(apply(matrix(A[T, j, ], nt), 2, function(x) rank(-x, ties.method = "random") <= kk))
      out[[length(out) + 1]] <- data.frame(Data = NAME, n_sub = z$n_sub[1], scheme = z$scheme[1], partition = z$partition[1],
        Env = envs[j], Model = m, n_test = nt, r = cor(mu, v), Precision = mean(top[sel]),
        Sel_diff = mean(v[sel]) - mean(v), Brier = mean((p_top - top)^2))
    }
  }
  o <- do.call(rbind, out)
  write.table(o, out_file, sep = ",", row.names = FALSE, col.names = !file.exists(out_file), append = file.exists(out_file))
  cat(NAME, key, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
