# Data-scarcity experiment: compares the REML FA(1) and FA(2)  (run_reml_scarcity.R,
# reml/results/reml_scarcity_<data set>.csv) with the Bayesian and BGLR models of gsreg_scarcity.R
# (scarcity_<data set>.csv) on identical testing sets, with the same metrics: selection differential
# of the top 10% of the predicted cells in each environment (SD units), precision and correlation,
# against the true genetic values in simulations and the phenotypes otherwise.
#
#   Rscript reml_scarcity_compare.R                 # REML files in reml/results/
#   Rscript reml_scarcity_compare.R <reml dir>      # REML files elsewhere (e.g. a test folder)
# REML FA(2) often fails to converge in the smallest sets, so "REML FA(2)/FA(1)" takes FA(2) where it
# converged and FA(1) otherwise, and every model is also compared with it.
# Without REML files for a data set, only the Bayesian and BGLR models are summarised.
# Writes reml_scarcity_comparison.csv.

out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))
args <- commandArgs(trailingOnly = TRUE)
reml_dir <- if (length(args)) args[1] else file.path(out_dir, "reml/results")
TOP <- 0.10
sets <- c("Wheat_599", "S1_bilinear_J4")
key <- c("Data", "n_sub", "scheme", "partition", "Env")

paired <- function(x, ref) {                 # paired differences of every model against `ref`
  b <- x[x$Model == ref, c(key, "S", "r")]
  out <- list()
  for (m in setdiff(unique(x$Model), c(ref, "GBLUP"))) {
    z <- merge(x[x$Model == m, c(key, "S", "r")], b, by = key, suffixes = c("", "_ref"))
    for (g in split(z, z[, c("Data", "n_sub", "scheme")], drop = TRUE)) {
      d <- g$S - g$S_ref
      out[[length(out) + 1]] <- data.frame(g[1, c("Data", "n_sub", "scheme")], Model = m, vs = ref,
        n_pairs = nrow(g), dS = round(mean(d), 3), rel_pct = round(100 * mean(d) / mean(g$S_ref), 1),
        p_S = signif(t.test(d)$p.value, 2), wins = round(mean(d > 0), 2), dr = round(mean(g$r - g$r_ref, na.rm = TRUE), 4),
        p_r = signif(t.test(g$r - g$r_ref)$p.value, 2))
    }
  }
  do.call(rbind, out)
}

res <- list()
for (nm in sets) {
  ba <- read.csv(file.path(out_dir, paste0("scarcity_", nm, ".csv")))
  ba$S <- ba$Sel_diff
  x <- ba[, c(key, "Model", "n_test", "S", "Precision", "r")]
  f <- file.path(reml_dir, paste0("reml_scarcity_", nm, ".csv"))
  if (file.exists(f)) {
    cs <- load_case(nm); Y <- cs$Y
    V <- if (is.null(cs$truth)) Y else cs$truth
    re <- read.csv(f, colClasses = c(line = "character", env = "character"))
    u <- unique(re[, c("n_sub", "scheme", "partition", "model", "converged")])
    cat(nm, "- REML fits:", nrow(u), "| not converged or failed:", sum(!u$converged), "\n")
    print(aggregate(cbind(failed = !converged) ~ n_sub + scheme + model, u, sum), row.names = FALSE)
    re <- re[re$converged & !is.na(re$line), ]
    re$v <- V[cbind(match(re$line, rownames(Y)), match(re$env, colnames(Y)))]
    stopifnot(!anyNA(re$v))
    rm_ <- do.call(rbind, lapply(split(re, list(re$n_sub, re$scheme, re$partition, re$model, re$env), drop = TRUE),
      function(z) {
        nt <- nrow(z); if (nt < 20) return(NULL); k <- ceiling(TOP * nt)
        sel <- rank(-z$blup, ties.method = "first") <= k
        top <- rank(-z$v, ties.method = "first") <= k
        data.frame(Data = nm, n_sub = z$n_sub[1], scheme = z$scheme[1], partition = z$partition[1], Env = z$env[1],
                   Model = z$model[1], n_test = nt, S = mean(z$v[sel]) - mean(z$v), Precision = mean(top[sel]),
                   r = cor(z$blup, z$v))
      }))
    # the testing cells must be the same as in the Bayesian runs
    chk <- merge(rm_[, c(key, "Model", "n_test")], unique(x[, c(key, "n_test")]), by = key, suffixes = c("", "_bayes"))
    if (any(chk$n_test != chk$n_test_bayes)) warning(nm, ": REML and Bayesian testing sets differ in ",
                                                     sum(chk$n_test != chk$n_test_bayes), " cases")
    # what a breeder would use: FA(2) where it converged, otherwise FA(1) of the same testing set
    ts <- paste(rm_$n_sub, rm_$scheme, rm_$partition)
    fb <- rm_[rm_$Model == "REML FA(2)" | (rm_$Model == "REML FA(1)" & !(ts %in% ts[rm_$Model == "REML FA(2)"])), ]
    fb$Model <- "REML FA(2)/FA(1)"
    x <- rbind(x, rm_, fb)
  } else cat(nm, "- no REML file in", reml_dir, "\n")
  res[[nm]] <- x
}
x <- do.call(rbind, res)

cat("\nMean selection differential, precision and correlation:\n")
cat("Testing sets with constant predictions (r undefined):\n")
print(aggregate(cbind(const = is.na(r)) ~ Data + n_sub + scheme + Model, x[is.na(x$r), ], sum), row.names = FALSE)
m <- aggregate(cbind(S, Precision, r) ~ Data + n_sub + scheme + Model, x, mean, na.rm = TRUE, na.action = na.pass)
print(m[order(m$Data, m$scheme, m$n_sub, -m$S), ], digits = 3, row.names = FALSE)
pd <- rbind(paired(x, "Reaction norm"), paired(x, "Bayes genomic FA K=2"))
if (any(x$Model == "REML FA(2)/FA(1)")) pd <- rbind(pd, paired(x, "REML FA(2)/FA(1)"))
cat("\nPaired differences (rel_pct = Model minus vs, % of vs):\n")
print(pd[order(pd$vs, pd$Data, pd$scheme, pd$n_sub, pd$Model), ], row.names = FALSE)
write.csv(merge(m, pd, all = TRUE), file.path(out_dir, "reml_scarcity_comparison.csv"), row.names = FALSE)
