# Convergence diagnostics of the Gibbs sampler for GAMMI (K = 2) and the Bayesian genomic FA
# (K = 2): 4 independent chains on Wheat 599 with 20% of the cells masked (CV2, partition 1,
# fold 1), so that predictions of untested cells are monitored too.
#
# Monitored quantities are identifiable: variance components, genetic correlations between
# environments (from s2u 11' + C C' + Psi, invariant to rotation of the components) and the
# predictions of 10 masked cells. Individual scores a_ik, c_jk are not identifiable and are
# not monitored.
# Diagnostics: split R-hat (Gelman-Rubin) and effective sample size (initial positive sequence).
#
#   Rscript gsreg_convergence.R 6000 1000          # iterations, burn-in (thin = 5)
# Writes convergence_<model>.csv and convergence_traces.png.

args <- commandArgs(trailingOnly = TRUE)
NITER <- if (length(args) >= 1) as.integer(args[1]) else 6000
BURN  <- if (length(args) >= 2) as.integer(args[2]) else 1000
THIN  <- 5; CHAINS <- 4
out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))

cs <- load_case("Wheat_599"); Y <- cs$Y; G <- cs$G; n <- nrow(Y); J <- ncol(Y)
set.seed(1001); invisible(sample(rep(1:5, length.out = n)))      # same draws as partition 1
fold_cell <- matrix(sample(rep(1:5, length.out = n * J)), n, J)
Ytr <- Y; Ytr[fold_cell == 1] <- NA
set.seed(7); cells <- which(fold_cell == 1)[sample(sum(fold_cell == 1), 10)]   # 10 masked cells
eigG <- eigen(G, symmetric = TRUE)

ess <- function(x) {                                 # effective sample size, initial positive sequence
  n <- length(x); if (sd(x) == 0) return(NA)
  r <- acf(x, lag.max = min(n - 1, 500), plot = FALSE)$acf[-1]
  s <- 0; for (l in seq(1, length(r) - 1, 2)) { p <- r[l] + r[l + 1]; if (p < 0) break; s <- s + p }
  n / (1 + 2 * s)
}
rhat <- function(M) {                                # split R-hat; M: draws x chains
  h <- floor(nrow(M) / 2); M <- cbind(M[1:h, , drop = FALSE], M[(h + 1):(2 * h), , drop = FALSE])
  m <- ncol(M); n <- nrow(M); W <- mean(apply(M, 2, var)); B <- n * var(colMeans(M))
  sqrt(((n - 1) / n * W + B / n) / W)
}

models <- list("GAMMI K=2" = list(specific = FALSE, main = TRUE),
               "Bayes genomic FA K=2" = list(specific = TRUE, main = FALSE))
traces <- list()
for (m in names(models)) {
  ch <- list()
  for (c in 1:CHAINS) {
    set.seed(100 * c)
    t0 <- Sys.time()
    P <- gsreg(Ytr, G, K = 2, specific = models[[m]]$specific, main = models[[m]]$main, niter = NITER, burn = BURN,
               eigK = eigG, keep_draws = TRUE, keep_trace = TRUE, thin = THIN)
    D <- attr(P, "draws"); pr <- t(apply(D, 3, function(F) F[cells]))
    colnames(pr) <- paste0("pred_", rownames(Y)[(cells - 1) %% n + 1], "_E", (cells - 1) %/% n + 1)
    ch[[c]] <- cbind(attr(P, "trace"), pr)
    cat(m, "chain", c, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
  traces[[m]] <- ch
  par_names <- colnames(ch[[1]])
  res <- t(sapply(par_names, function(p) {
    M <- sapply(ch, function(x) x[, p]); v <- as.vector(M)
    c(mean = mean(v), sd = sd(v), q2.5 = unname(quantile(v, 0.025)), q97.5 = unname(quantile(v, 0.975)),
      Rhat = rhat(M), ESS = sum(apply(M, 2, ess)))
  }))
  res <- data.frame(parameter = rownames(res), round(res, 4), row.names = NULL)
  write.csv(res, file.path(out_dir, paste0("convergence_", gsub("[ =]", "_", m), ".csv")), row.names = FALSE)
  cat("\n", m, ": max R-hat =", round(max(res$Rhat, na.rm = TRUE), 3), " min ESS =", round(min(res$ESS, na.rm = TRUE)),
      " (draws kept =", CHAINS * nrow(ch[[1]]), ")\n"); print(res, row.names = FALSE)
}

# trace plots of six quantities per model, chains in colour
cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
lab <- c("GAMMI K=2" = "GAMMI", "Bayes genomic FA K=2" = "GSREG-ESE")   # model names used in the manuscript
png(file.path(out_dir, "convergence_traces.png"), width = 2400, height = 1500, res = 170)
par(mfrow = c(2, 6), mar = c(3.5, 3.5, 2.5, 0.5), mgp = c(2.2, 0.7, 0))
for (m in names(traces)) {
  ch <- traces[[m]]; nm <- colnames(ch[[1]])
  show <- c("sigma2_1", "tau2_1", grep("^rg_12$|^rg_13$|^rg_23$", nm, value = TRUE)[1:3], grep("^pred_", nm, value = TRUE)[1])
  for (p in show) {
    rg <- range(sapply(ch, function(x) range(x[, p])))
    plot(ch[[1]][, p], type = "l", col = cols[1], ylim = rg, xlab = "Saved draw", ylab = p,
         main = paste0(lab[m], ": ", p), cex.main = 0.85)
    for (c in 2:length(ch)) lines(ch[[c]][, p], col = cols[c])
  }
}
dev.off()
cat("\nWrote convergence_*.csv and convergence_traces.png\n")
