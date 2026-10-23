# Selection differential of the lines each model selects (top 10% of the testing lines):
# mean value of the selected lines minus the mean of all testing lines, in units of the
# phenotypic SD of the environment (phenotypes are standardised within environment).
# Uses the true genetic values in simulations and the observed phenotypes otherwise.
# Less noisy than counting hits in the top 10%, because every selected line contributes its value.
#
#   Rscript gsreg_decision_gain.R Japonica Wheat_599 S1_bilinear_J4 S3_rnorm_J4

out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))
args <- commandArgs(trailingOnly = TRUE)
res <- list()
for (nm in args) {
  li <- read.csv(file.path(out_dir, paste0("decision_", nm, "_lines.csv")))
  cs <- load_case(nm)
  V <- if (is.null(cs$truth)) cs$Y else cs$truth                         # value the selection is judged on
  li$value <- V[cbind(match(li$Line, rownames(cs$Y)), match(li$Env, colnames(cs$Y)))]
  g <- aggregate(value ~ Data + CV + Rep + Fold + Env + Model, li[li$Selected, ], mean)
  b <- aggregate(value ~ Data + CV + Rep + Fold + Env + Model, li, mean)
  g$S <- g$value - b$value[match(do.call(paste, g[, 1:6]), do.call(paste, b[, 1:6]))]
  res[[nm]] <- g
}
g <- do.call(rbind, res)
s <- aggregate(S ~ Data + CV + Model, g, mean); s$S <- round(s$S, 3)
s <- s[order(s$Data, s$CV, -s$S), ]
cat("Mean selection differential of the selected lines (SD units):\n"); print(s, row.names = FALSE)
key <- c("Data", "CV", "Rep", "Fold", "Env")
rn <- g[g$Model == "Reaction norm", c(key, "S")]
out <- list()
for (m in setdiff(unique(g$Model), "Reaction norm")) {
  x <- merge(g[g$Model == m, c(key, "S")], rn, by = key, suffixes = c("", "_RN"))
  for (z in split(x, list(x$Data, x$CV), drop = TRUE)) {
    d <- z$S - z$S_RN
    out[[length(out) + 1]] <- data.frame(Data = z$Data[1], CV = z$CV[1], Model = m, n_pairs = nrow(z),
      dS = round(mean(d), 3), rel_gain_pct = round(100 * mean(d) / mean(z$S_RN), 1),
      p = signif(t.test(d)$p.value, 2), wins = round(mean(d > 0), 2))
  }
}
out <- do.call(rbind, out); out <- out[order(out$Data, out$CV, out$Model), ]
cat("\nPaired difference against the reaction norm:\n"); print(out, row.names = FALSE)
write.csv(out, file.path(out_dir, "decision_gain_paired.csv"), row.names = FALSE)
