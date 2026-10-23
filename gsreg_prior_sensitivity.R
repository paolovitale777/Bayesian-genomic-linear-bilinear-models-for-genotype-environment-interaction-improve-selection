# Prior sensitivity of GAMMI (K = 2) and the Bayesian genomic FA (K = 2) on Wheat 599:
# the same partitions are refitted with the default, a diffuse and an informative set of
# inverse-gamma priors, and accuracy and selection differential are compared.
#
#   default:     tau2 ~ IG(3, 1),    sigma2_j ~ IG(1, 0.5),    psi_j ~ IG(3, 0.2),  sigma2_u ~ IG(3, 0.5)
#   diffuse:     tau2 ~ IG(1, 0.1),  sigma2_j ~ IG(0.01, 0.01), psi_j ~ IG(1, 0.01), sigma2_u ~ IG(1, 0.1)
#   informative: tau2 ~ IG(10, 5),   sigma2_j ~ IG(10, 5),      psi_j ~ IG(10, 1),   sigma2_u ~ IG(10, 3)
# Schemes: CV1 and CV2 of partition 1 (5 folds each) and sparse testing with k = 1 (partitions
# 1-3), where the priors matter most because each line has a single environment.
#
#   Rscript gsreg_prior_sensitivity.R 2000 500
# Writes prior_sensitivity.csv and prints the paired differences against the default priors.

args <- commandArgs(trailingOnly = TRUE)
NITER <- if (length(args) >= 1) as.integer(args[1]) else 2000
BURN  <- if (length(args) >= 2) as.integer(args[2]) else 500
out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))

cs <- load_case("Wheat_599"); Y <- cs$Y; G <- cs$G; n <- nrow(Y); J <- ncol(Y); envs <- colnames(Y)
eigG <- eigen(G, symmetric = TRUE)
priors <- list(default     = list(at = 3,  bt = 1,   as0 = 1,    bs0 = 0.5,  ap = 3,  bp = 0.2,  ag = 3,  bg = 0.5),
               diffuse     = list(at = 1,  bt = 0.1, as0 = 0.01, bs0 = 0.01, ap = 1,  bp = 0.01, ag = 1,  bg = 0.1),
               informative = list(at = 10, bt = 5,   as0 = 10,   bs0 = 5,    ap = 10, bp = 1,    ag = 10, bg = 3))
models <- list("GAMMI K=2" = c(specific = FALSE, main = TRUE), "Bayes genomic FA K=2" = c(specific = TRUE, main = FALSE))

# testing sets: identical to gsreg_decision.R (partition 1) and gsreg_sparse.R (k = 1)
sets <- list()
set.seed(1001); fold_line <- sample(rep(1:5, length.out = n)); fold_cell <- matrix(sample(rep(1:5, length.out = n * J)), n, J)
for (f in 1:5) {
  sets[[length(sets) + 1]] <- list(scheme = "CV1", id = f, test = matrix(fold_line == f, n, J))
  sets[[length(sets) + 1]] <- list(scheme = "CV2", id = f, test = fold_cell == f)
}
for (rep in 1:3) {
  set.seed(2000 + 10 * rep + 1); Ytr <- Y
  for (i in 1:n) { o <- which(!is.na(Y[i, ])); if (length(o) > 1) Ytr[i, o[-sample.int(length(o), 1)]] <- NA }
  sets[[length(sets) + 1]] <- list(scheme = "sparse_k1", id = rep, test = is.na(Ytr) & !is.na(Y))
}

out <- list(); t0 <- Sys.time()
for (s in sets) {
  Ytr <- Y; Ytr[s$test] <- NA
  for (m in names(models)) for (pn in names(priors)) {
    pr <- priors[[pn]]
    P <- do.call(gsreg, c(list(Y = Ytr, Kmat = G, K = 2, specific = models[[m]][["specific"]], main = models[[m]][["main"]],
                               niter = NITER, burn = BURN, eigK = eigG), pr))
    for (j in 1:J) {
      T <- which(s$test[, j] & !is.na(Y[, j])); if (length(T) < 20) next
      kk <- ceiling(0.1 * length(T)); sel <- rank(-P[T, j], ties.method = "first") <= kk
      out[[length(out) + 1]] <- data.frame(Scheme = s$scheme, Set = s$id, Env = envs[j], Model = m, Prior = pn,
        r = cor(P[T, j], Y[T, j]), Sel_diff = mean(Y[T, j][sel]) - mean(Y[T, j]))
    }
  }
  cat(s$scheme, s$id, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  write.csv(do.call(rbind, out), file.path(out_dir, "prior_sensitivity.csv"), row.names = FALSE)
}

o <- do.call(rbind, out)
s <- aggregate(cbind(r, Sel_diff) ~ Scheme + Model + Prior, o, mean)
print(s[order(s$Scheme, s$Model, s$Prior), ], digits = 3, row.names = FALSE)
key <- c("Scheme", "Set", "Env", "Model")
d0 <- o[o$Prior == "default", c(key, "r", "Sel_diff")]
res <- list()
for (pn in c("diffuse", "informative")) {
  x <- merge(o[o$Prior == pn, c(key, "r", "Sel_diff")], d0, by = key, suffixes = c("", "_def"))
  for (z in split(x, list(x$Scheme, x$Model), drop = TRUE))
    res[[length(res) + 1]] <- data.frame(Scheme = z$Scheme[1], Model = z$Model[1], Prior = pn, n = nrow(z),
      dr = round(mean(z$r - z$r_def), 4), p_r = signif(t.test(z$r - z$r_def)$p.value, 2),
      dSel = round(mean(z$Sel_diff - z$Sel_diff_def), 4), p_Sel = signif(t.test(z$Sel_diff - z$Sel_diff_def)$p.value, 2))
}
cat("\nDifferences against the default priors:\n"); print(do.call(rbind, res), row.names = FALSE)
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
