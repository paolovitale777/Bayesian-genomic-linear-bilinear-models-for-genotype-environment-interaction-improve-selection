# Summary of the decision and calibration analysis (gsreg_decision.R) over all data sets found.
#
#   Rscript gsreg_decision_summary.R
#
# Writes decision_summary.csv (per data set x CV x model), decision_paired.csv (paired
# differences against the reaction norm over partitions x folds x environments) and
# decision_calibration.png (reliability curves of P(top 10%)).

out_dir <- "."
args <- commandArgs(trailingOnly = TRUE)                  # optional: data sets to include, e.g. Japonica Wheat_599
ff <- list.files(out_dir, "^decision_.*_folds\\.csv$", full.names = TRUE)
if (length(args)) ff <- ff[sub("^decision_(.*)_folds\\.csv$", "\\1", basename(ff)) %in% args]
fo <- do.call(rbind, lapply(ff, read.csv))
li <- do.call(rbind, lapply(sub("_folds", "_lines", ff), read.csv))
sim <- grepl("^S[0-9]_", fo$Data)
# outcome used everywhere: the true genetic value in simulations, the phenotype otherwise
fo$Precision <- ifelse(sim, fo$Precision_true, fo$Precision_obs)
fo$Brier     <- ifelse(sim, fo$Brier_true, fo$Brier_obs)
fo$r         <- ifelse(sim, fo$r_true, fo$r_obs)
simL <- grepl("^S[0-9]_", li$Data)
li$Top <- ifelse(simL, li$Top_true, li$Top_obs)
li$P_top <- ifelse(simL, li$P_top_gen, li$P_top)        # genetic probability vs. truth; predictive vs. phenotype
base <- 0.1 * (1 - 0.1)                                   # Brier score of always predicting 0.10

s <- aggregate(cbind(Precision, Brier, r) ~ Data + CV + Model, fo, mean)
s$BSS <- 1 - s$Brier / base                               # Brier skill score vs. the base rate
s$n_cases <- aggregate(Precision ~ Data + CV + Model, fo, length)$Precision
s <- s[order(s$Data, s$CV, -s$Precision), ]
s[, c("Precision", "Brier", "r", "BSS")] <- round(s[, c("Precision", "Brier", "r", "BSS")], 3)
write.csv(s, file.path(out_dir, "decision_summary.csv"), row.names = FALSE)
print(s, row.names = FALSE)

# paired differences against the reaction norm on identical partitions, folds and environments
key <- c("Data", "CV", "Rep", "Fold", "Env")
rn <- fo[fo$Model == "Reaction norm", c(key, "Precision", "Brier", "r")]
pr <- list()
for (m in setdiff(unique(fo$Model), "Reaction norm")) {
  x <- merge(fo[fo$Model == m, c(key, "Precision", "Brier", "r")], rn, by = key, suffixes = c("", "_RN"))
  for (g in split(x, list(x$Data, x$CV), drop = TRUE)) {
    test <- function(a, b) { d <- a - b; if (sd(d) == 0) return(c(mean(d), NA, mean(d > 0))); c(mean(d), t.test(d)$p.value, mean(d > 0)) }
    tp <- test(g$Precision, g$Precision_RN); tb <- test(-g$Brier, -g$Brier_RN); tr <- test(g$r, g$r_RN)
    pr[[length(pr) + 1]] <- data.frame(Data = g$Data[1], CV = g$CV[1], Model = m, n_pairs = nrow(g),
      dPrecision = tp[1], p_Precision = tp[2], win_Precision = tp[3],
      dBrier = -tb[1], p_Brier = tb[2], dr = tr[1], p_r = tr[2], win_r = tr[3])
  }
}
pr <- do.call(rbind, pr); num <- sapply(pr, is.numeric); pr[num] <- lapply(pr[num], function(v) signif(v, 3))
write.csv(pr, file.path(out_dir, "decision_paired.csv"), row.names = FALSE)
cat("\nPaired differences against the reaction norm (positive dPrecision and dr, negative dBrier = better):\n")
print(pr, row.names = FALSE)

# reliability curves: predicted P(top 10%) in bins vs. observed frequency, by model (all data pooled)
bins <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1)
li$bin <- cut(li$P_top, bins, include.lowest = TRUE)
cal <- aggregate(cbind(P_top, Top) ~ Model + CV + bin, li, mean)
cal$n <- aggregate(P_top ~ Model + CV + bin, li, length)$P_top
write.csv(cal, file.path(out_dir, "decision_calibration.csv"), row.names = FALSE)
cols <- c("GBLUP" = "#999999", "Reaction norm" = "#0072B2", "GAMMI K=2" = "#D55E00", "Bayes genomic FA K=2" = "#009E73")
lab <- c("GBLUP", "RN", "GAMMI (K = 2)", "GSREG-ESE (K = 2)")      # model names used in the manuscript
png(file.path(out_dir, "decision_calibration.png"), width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (cv in c("CV1", "CV2")) {
  plot(c(0, 1), c(0, 1), type = "n", xlab = "Predicted P(top 10%)", ylab = "Observed frequency in top 10%",
       main = paste(cv, "- calibration (all data sets)"), asp = 1)
  abline(0, 1, col = "grey70", lty = 2)
  for (m in names(cols)) { z <- cal[cal$Model == m & cal$CV == cv, ]; z <- z[order(z$P_top), ]
    lines(z$P_top, z$Top, col = cols[m], lwd = 2, type = "b", pch = 16) }
  legend("topleft", lab, col = cols, lwd = 2, bty = "n", cex = 0.8)
}
dev.off()
cat("\nWrote decision_summary.csv, decision_paired.csv, decision_calibration.csv and decision_calibration.png\n")
