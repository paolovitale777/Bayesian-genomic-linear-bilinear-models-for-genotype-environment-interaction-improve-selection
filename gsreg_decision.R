# Selection decisions and calibration of selection probabilities for untested lines (CV1)
# and incomplete trials (CV2), with repeated random partitions.
#
# For every environment and testing set, each model
#   - selects the top 10% of the testing lines by its predicted mean, and
#   - gives each testing line a probability P(top 10%) of ranking among the best 10%.
# Two probabilities: P_top_gen (genetic value in the top 10%, compared with the true values in
# simulations) and P_top (posterior predictive: phenotype in the top 10%, draws plus residual
# error, compared with the observed phenotypes).
# Bayesian genomic models (GAMMI K = 2, genomic FA K = 2): probabilities from posterior draws.
# GBLUP and reaction norm (BGLR): BGLR keeps no draws of the predictions, so probabilities
# use independent normal draws N(yHat, SD.yHat) (+ varE for P_top) - a plug-in approximation.
#
# Outputs per data set: decision_<name>_lines.csv (one row per line x environment x fold x
# partition x model: probability and outcomes) and decision_<name>_folds.csv (precision,
# Brier score and accuracy per environment x fold x partition x model).
#
#   Rscript gsreg_decision.R Wheat_599 10 2000 500       # name, partitions, iterations, burn-in
#   Rscript gsreg_decision.R EYT_16_17 3 2000 500 2      # resume: run partitions 2..3, keep partition 1
#   names: Wheat_599, EYT_16_17, EYT_22_23, Japonica, S1_bilinear_J4, S3_rnorm_J4

args <- commandArgs(trailingOnly = TRUE)
NAME  <- args[1]
REPS  <- if (length(args) >= 2) as.integer(args[2]) else 10
NITER <- if (length(args) >= 3) as.integer(args[3]) else 2000
BURN  <- if (length(args) >= 4) as.integer(args[4]) else 500
START <- if (length(args) >= 5) as.integer(args[5]) else 1
THIN  <- 5
TOP   <- 0.10
out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))

cs <- load_case(NAME); Y <- cs$Y; G <- cs$G; truth <- cs$truth
n <- nrow(Y); J <- ncol(Y); envs <- colnames(Y)
cat(NAME, "- lines:", n, " environments:", J, " partitions:", REPS, "\n")

eigG <- eigen(G, symmetric = TRUE)
Zg <- model.matrix(~ factor(rep(1:n, J)) - 1); Ze <- model.matrix(~ factor(rep(1:J, each = n)) - 1)
KG <- Zg %*% G %*% t(Zg)
eig <- function(K) { e <- eigen(K, symmetric = TRUE); keep <- e$values > 1e-8; list(V = e$vectors[, keep], d = e$values[keep]) }
EG <- eig(KG); EGE <- eig(KG * tcrossprod(Ze)); rm(KG); invisible(gc())
fit_bglr <- function(yNA, gxe) {
  ETA <- list(E = list(X = Ze[, -1, drop = FALSE], model = "FIXED"), G = list(V = EG$V, d = EG$d, model = "RKHS"))
  if (gxe) ETA$GE <- list(V = EGE$V, d = EGE$d, model = "RKHS")
  fm <- BGLR(y = yNA, ETA = ETA, nIter = NITER * 2, burnIn = BURN, verbose = FALSE,
             saveAt = file.path(tempdir(), paste0("dec_", NAME, "_", gxe, "_")))
  S <- (NITER - BURN) / THIN
  mu <- matrix(fm$yHat, n, J); sdv <- matrix(fm$SD.yHat, n, J)
  list(mean = mu,
       draws = array(rnorm(n * J * S, mu, sdv), c(n, J, S)),                        # genetic values
       pred  = array(rnorm(n * J * S, mu, sqrt(sdv^2 + fm$varE)), c(n, J, S)))     # + residual: phenotypes
}
fit_bayes <- function(Ytr, specific) {
  P <- gsreg(Ytr, G, K = 2, specific = specific, main = !specific, niter = NITER, burn = BURN,
             eigK = eigG, keep_draws = TRUE, thin = THIN)
  D <- attr(P, "draws"); s2 <- attr(P, "s2"); S <- dim(D)[3]
  E <- array(rnorm(n * J * S), c(n, J, S)) * aperm(array(sqrt(s2), c(S, J, n)), c(3, 2, 1))
  list(mean = as.matrix(P[, ]), draws = D, pred = D + E)                            # posterior predictive
}
models <- list("GBLUP"                 = function(Ytr) fit_bglr(as.vector(Ytr), FALSE),
               "Reaction norm"         = function(Ytr) fit_bglr(as.vector(Ytr), TRUE),
               "GAMMI K=2"             = function(Ytr) fit_bayes(Ytr, FALSE),
               "Bayes genomic FA K=2"  = function(Ytr) fit_bayes(Ytr, TRUE))

lines_out <- list(); folds_out <- list()
if (START > 1) {   # keep the partitions already saved (seeds are per partition, so splits are unchanged)
  prev <- function(kind) {
    cc <- if (kind == "lines") c(Line = "character") else NA
    d <- read.csv(file.path(out_dir, paste0("decision_", NAME, "_", kind, ".csv")), colClasses = cc)
    d[d$Rep < START, ]
  }
  lines_out[[1]] <- prev("lines"); folds_out[[1]] <- prev("folds")
  cat("Resuming at partition", START, "- kept", nrow(folds_out[[1]]), "fold rows\n")
}
t0 <- Sys.time()
for (rep in START:REPS) {
  set.seed(1000 + rep)
  fold_line <- sample(rep(1:5, length.out = n))
  fold_cell <- matrix(sample(rep(1:5, length.out = n * J)), n, J)
  for (cv in c("CV1", "CV2")) {
    for (f in 1:5) {
      test <- if (cv == "CV1") matrix(fold_line == f, n, J) else (fold_cell == f)
      Ytr <- Y; Ytr[test] <- NA
      for (m in names(models)) {
        fit <- models[[m]](Ytr)
        for (j in 1:J) {
          T <- which(test[, j] & !is.na(Y[, j])); nt <- length(T)
          if (nt < 20) next
          k <- ceiling(TOP * nt)
          yo <- Y[T, j]; mu <- fit$mean[T, j]
          top_obs  <- rank(-yo, ties.method = "first") <= k
          top_true <- if (is.null(truth)) rep(NA, nt) else rank(-truth[T, j], ties.method = "first") <= k
          sel <- rank(-mu, ties.method = "first") <= k
          ptop <- function(A) rowMeans(apply(matrix(A[T, j, ], nt), 2, function(v) rank(-v, ties.method = "random") <= k))
          p_pred <- ptop(fit$pred)     # posterior predictive: probability that the PHENOTYPE is in the top 10%
          p_gen  <- ptop(fit$draws)    # probability that the GENETIC VALUE is in the top 10%
          lines_out[[length(lines_out) + 1]] <- data.frame(Data = NAME, CV = cv, Rep = rep, Fold = f, Env = envs[j],
            Model = m, Line = rownames(Y)[T], P_top = round(p_pred, 4), P_top_gen = round(p_gen, 4),
            Selected = sel, Top_obs = top_obs, Top_true = top_true)
          folds_out[[length(folds_out) + 1]] <- data.frame(Data = NAME, CV = cv, Rep = rep, Fold = f, Env = envs[j],
            Model = m, n_test = nt, k = k,
            Precision_obs = mean(top_obs[sel]),
            Precision_true = if (is.null(truth)) NA else mean(top_true[sel]),
            Brier_obs = mean((p_pred - top_obs)^2),
            Brier_true = if (is.null(truth)) NA else mean((p_gen - top_true)^2),
            r_obs = cor(mu, yo), r_true = if (is.null(truth)) NA else cor(mu, truth[T, j]))
        }
      }
      cat(NAME, "partition", rep, cv, "fold", f, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
    }
  }
  # write after every partition, so partial results survive an interruption
  write.csv(do.call(rbind, lines_out), file.path(out_dir, paste0("decision_", NAME, "_lines.csv")), row.names = FALSE)
  write.csv(do.call(rbind, folds_out), file.path(out_dir, paste0("decision_", NAME, "_folds.csv")), row.names = FALSE)
}
fo <- do.call(rbind, folds_out)
print(aggregate(cbind(Precision_obs, Brier_obs, r_obs) ~ CV + Model, fo, mean), digits = 3)
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
