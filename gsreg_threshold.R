# Selection of the top 10% as a classification problem (Montesinos-López et al., 2023a, 2025), using only the
# saved results of gsreg_decision.R (decision_<data set>_lines.csv): no model is refitted. Three rules:
#   topk:    the 10% of lines with the highest predicted mean (as in the paper so far)
#   balance: lines with P(top 10%) > t, with t chosen so that sensitivity = specificity (as the RO method)
#   maxF1:   lines with P(top 10%) > t, with t chosen to maximise F1
# The thresholds are chosen on the other partitions of the same data set, CV scheme and model (leave one
# partition out), so they never see the testing set. Outcome: the true genetic value in simulations and the
# phenotype otherwise, with the matching probability (genetic or predictive), as in gsreg_decision_summary.R.
#
#   Rscript gsreg_threshold.R
# Writes threshold_by_set.csv (per testing set), threshold_summary.csv and threshold_paired.csv.

out_dir <- "."
ff <- list.files(out_dir, "^decision_.*_lines\\.csv$", full.names = TRUE)
li <- do.call(rbind, lapply(ff, read.csv, colClasses = c(Line = "character")))
sim <- grepl("^S[0-9]_", li$Data)
li$Top <- ifelse(sim, li$Top_true, li$Top_obs)
li$P <- ifelse(sim, li$P_top_gen, li$P_top)
li <- li[!is.na(li$Top) & !is.na(li$P), ]
grid <- seq(0.01, 0.99, 0.01)

conf <- function(sel, top) {                      # sensitivity, specificity, F1, Cohen's kappa
  tp <- sum(sel & top); fp <- sum(sel & !top); fn <- sum(!sel & top); tn <- sum(!sel & !top); n <- tp + fp + fn + tn
  se <- tp / (tp + fn); sp <- tn / (tn + fp); pr <- if (tp + fp > 0) tp / (tp + fp) else 0
  f1 <- if (pr + se > 0) 2 * pr * se / (pr + se) else 0
  po <- (tp + tn) / n; pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
  c(Se = se, Sp = sp, F1 = f1, Kappa = if (pe < 1) (po - pe) / (1 - pe) else 0, Selected = (tp + fp) / n)
}

res <- list()
for (g in split(li, list(li$Data, li$CV, li$Model), drop = TRUE)) {
  reps <- sort(unique(g$Rep))
  for (r in reps) {
    tr <- g[g$Rep != r, ]; te <- g[g$Rep == r, ]
    # thresholds from the other partitions, pooled over their folds and environments
    m <- sapply(grid, function(t) conf(tr$P > t, tr$Top))
    t_bal <- grid[which.min(abs(m["Se", ] - m["Sp", ]))]
    t_f1 <- grid[which.max(m["F1", ])]
    for (s in split(te, list(te$Fold, te$Env), drop = TRUE)) {
      for (rule in c("topk", "balance", "maxF1")) {
        sel <- switch(rule, topk = s$Selected, balance = s$P > t_bal, maxF1 = s$P > t_f1)
        res[[length(res) + 1]] <- data.frame(Data = s$Data[1], CV = s$CV[1], Model = s$Model[1], Rep = r,
          Fold = s$Fold[1], Env = s$Env[1], Rule = rule,
          t = switch(rule, topk = NA, balance = t_bal, maxF1 = t_f1), t(conf(sel, s$Top)))
      }
    }
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_dir, "threshold_by_set.csv"), row.names = FALSE)

met <- c("Se", "Sp", "F1", "Kappa", "Selected")
sm <- aggregate(res[, met], res[, c("CV", "Model", "Rule")], mean)
sm <- sm[order(sm$CV, sm$Rule, -sm$F1), ]
sm[met] <- round(sm[met], 3)
write.csv(sm, file.path(out_dir, "threshold_summary.csv"), row.names = FALSE)
cat("Mean over the six data sets and all testing sets:\n"); print(sm, row.names = FALSE)
sd_ <- aggregate(res[, met], res[, c("Data", "CV", "Model", "Rule")], mean)
sd_[met] <- round(sd_[met], 3)
write.csv(sd_, file.path(out_dir, "threshold_by_data.csv"), row.names = FALSE)

# paired comparisons on identical testing sets: (i) each rule against topk within a model;
# (ii) each model against the reaction norm within a rule
key <- c("Data", "CV", "Rep", "Fold", "Env")
pt <- function(a, b, what) { d <- a - b; c(d = mean(d), p = if (sd(d) > 0) t.test(d)$p.value else NA) }
out <- list()
for (cv in unique(res$CV)) for (m in unique(res$Model)) for (rule in c("balance", "maxF1")) {
  x <- merge(res[res$CV == cv & res$Model == m & res$Rule == rule, c(key, "F1", "Kappa")],
             res[res$CV == cv & res$Model == m & res$Rule == "topk", c(key, "F1", "Kappa")], by = key, suffixes = c("", "_ref"))
  out[[length(out) + 1]] <- data.frame(CV = cv, Comparison = paste(rule, "vs topk"), Model = m, n = nrow(x),
    dF1 = pt(x$F1, x$F1_ref)["d"], p_F1 = pt(x$F1, x$F1_ref)["p"], dKappa = pt(x$Kappa, x$Kappa_ref)["d"],
    p_Kappa = pt(x$Kappa, x$Kappa_ref)["p"])
}
for (cv in unique(res$CV)) for (rule in c("topk", "balance", "maxF1")) for (m in setdiff(unique(res$Model), "Reaction norm")) {
  x <- merge(res[res$CV == cv & res$Model == m & res$Rule == rule, c(key, "F1", "Kappa")],
             res[res$CV == cv & res$Model == "Reaction norm" & res$Rule == rule, c(key, "F1", "Kappa")], by = key,
             suffixes = c("", "_ref"))
  out[[length(out) + 1]] <- data.frame(CV = cv, Comparison = paste(m, "vs RN,", rule), Model = m, n = nrow(x),
    dF1 = pt(x$F1, x$F1_ref)["d"], p_F1 = pt(x$F1, x$F1_ref)["p"], dKappa = pt(x$Kappa, x$Kappa_ref)["d"],
    p_Kappa = pt(x$Kappa, x$Kappa_ref)["p"])
}
out <- do.call(rbind, out); num <- sapply(out, is.numeric); out[num] <- lapply(out[num], function(v) signif(v, 3))
write.csv(out, file.path(out_dir, "threshold_paired.csv"), row.names = FALSE)
cat("\nPaired comparisons (positive = better than the reference):\n"); print(out, row.names = FALSE)
