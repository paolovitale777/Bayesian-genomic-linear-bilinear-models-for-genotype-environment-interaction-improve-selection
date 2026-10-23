# Sparse testing: each line is phenotyped in only k of the J environments (chosen at random,
# among the environments where it was observed); the models predict the environments where
# the line was not tested. Hypothesis: the advantage of the genomic bilinear / FA models over
# the reaction norm grows as k decreases.
#
# For every k, partition and environment: correlation, precision of selecting the top 10% of
# the untested lines, selection differential of the selected lines (SD units; true genetic
# values in simulations, phenotypes otherwise) and Brier score of P(top 10%).
#
#   Rscript gsreg_sparse.R Wheat_599 5 2000 500            # name, partitions, iterations, burn-in
#   Rscript gsreg_sparse.R EYT_16_17 3 2000 500 1,2,3      # optional: values of k
# Writes sparse_<name>.csv (one row per k x partition x environment x model).

args <- commandArgs(trailingOnly = TRUE)
NAME  <- args[1]
REPS  <- if (length(args) >= 2) as.integer(args[2]) else 5
NITER <- if (length(args) >= 3) as.integer(args[3]) else 2000
BURN  <- if (length(args) >= 4) as.integer(args[4]) else 500
TOP   <- 0.10
out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))

cs <- load_case(NAME); Y <- cs$Y; G <- cs$G; truth <- cs$truth
n <- nrow(Y); J <- ncol(Y); envs <- colnames(Y)
KS <- if (length(args) >= 5) as.integer(strsplit(args[5], ",")[[1]]) else 1:(J - 1)
V <- if (is.null(truth)) Y else truth                       # value the selection is judged on
cat(NAME, "- lines:", n, " environments:", J, " k:", KS, " partitions:", REPS, "\n")
models <- make_fitters(Y, G, niter = NITER, burn = BURN, label = paste0("sparse_", NAME))

out <- list(); t0 <- Sys.time()
for (rep in 1:REPS) {
  for (k in KS) {
    set.seed(2000 + 10 * rep + k)
    Ytr <- Y
    for (i in 1:n) {                                         # keep k observed environments per line
      o <- which(!is.na(Y[i, ]))
      if (length(o) > k) Ytr[i, o[-sample.int(length(o), k)]] <- NA
    }
    test <- is.na(Ytr) & !is.na(Y)
    for (m in names(models)) {
      fit <- models[[m]](Ytr)
      for (j in 1:J) {
        T <- which(test[, j]); nt <- length(T)
        if (nt < 20) next
        kk <- ceiling(TOP * nt)
        mu <- fit$mean[T, j]; yo <- Y[T, j]; v <- V[T, j]
        sel <- rank(-mu, ties.method = "first") <= kk
        top_v <- rank(-v, ties.method = "first") <= kk
        A <- if (is.null(truth)) fit$pred else fit$draws     # predictive vs. phenotype; genetic vs. truth
        p_top <- rowMeans(apply(matrix(A[T, j, ], nt), 2, function(x) rank(-x, ties.method = "random") <= kk))
        out[[length(out) + 1]] <- data.frame(Data = NAME, k = k, Tested = round(k / J, 2), Rep = rep, Env = envs[j],
          Model = m, n_test = nt, r = cor(mu, v), Precision = mean(top_v[sel]),
          Sel_diff = mean(v[sel]) - mean(v), Brier = mean((p_top - top_v)^2))
      }
    }
    cat(NAME, "partition", rep, "k =", k, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
    write.csv(do.call(rbind, out), file.path(out_dir, paste0("sparse_", NAME, ".csv")), row.names = FALSE)
  }
}

# summary: means by k and model, and paired differences against the reaction norm
o <- do.call(rbind, out)
s <- aggregate(cbind(r, Precision, Sel_diff, Brier) ~ k + Tested + Model, o, mean)
print(s[order(s$k, -s$Sel_diff), ], digits = 3, row.names = FALSE)
key <- c("k", "Rep", "Env")
rn <- o[o$Model == "Reaction norm", c(key, "Sel_diff", "r")]
pr <- list()
for (m in setdiff(names(models), "Reaction norm")) {
  x <- merge(o[o$Model == m, c(key, "Sel_diff", "r")], rn, by = key, suffixes = c("", "_RN"))
  for (z in split(x, x$k)) {
    d <- z$Sel_diff - z$Sel_diff_RN
    pr[[length(pr) + 1]] <- data.frame(k = z$k[1], Model = m, n_pairs = nrow(z), dSel = round(mean(d), 3),
      rel_gain_pct = round(100 * mean(d) / mean(z$Sel_diff_RN), 1), p = signif(t.test(d)$p.value, 2),
      dr = round(mean(z$r - z$r_RN), 3))
  }
}
cat("\nPaired difference against the reaction norm (selection differential and correlation):\n")
print(do.call(rbind, pr), row.names = FALSE)
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
