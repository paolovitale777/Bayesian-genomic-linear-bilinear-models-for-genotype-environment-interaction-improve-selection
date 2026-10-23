# Summary of the sparse-testing analysis (gsreg_sparse.R) for all data sets found:
# mean selection differential and correlation by k (environments tested per line) and model,
# and paired differences against the reaction norm over partitions x environments.
#
#   Rscript gsreg_sparse_summary.R

out_dir <- "."
ff <- list.files(out_dir, "^sparse_.*\\.csv$", full.names = TRUE)
o <- do.call(rbind, lapply(ff[!grepl("sparse_summary", ff)], read.csv))
key <- c("Data", "k", "Rep", "Env")
rn <- o[o$Model == "Reaction norm", c(key, "Sel_diff", "r", "Brier")]
pr <- list()
for (m in setdiff(unique(o$Model), "Reaction norm")) {
  x <- merge(o[o$Model == m, c(key, "Tested", "Sel_diff", "r", "Brier")], rn, by = key, suffixes = c("", "_RN"))
  for (z in split(x, list(x$Data, x$k), drop = TRUE)) {
    d <- z$Sel_diff - z$Sel_diff_RN
    pr[[length(pr) + 1]] <- data.frame(Data = z$Data[1], k = z$k[1], Tested = z$Tested[1], Model = m,
      n_pairs = nrow(z), Sel_RN = round(mean(z$Sel_diff_RN), 3), Sel_model = round(mean(z$Sel_diff), 3),
      gain_pct = round(100 * mean(d) / mean(z$Sel_diff_RN), 1),
      p = if (nrow(z) > 2 && sd(d) > 0) signif(t.test(d)$p.value, 2) else NA,
      wins = round(mean(d > 0), 2), dr = round(mean(z$r - z$r_RN), 3), dBrier = signif(mean(z$Brier - z$Brier_RN), 2))
  }
}
pr <- do.call(rbind, pr); pr <- pr[order(pr$Data, pr$Model, pr$k), ]
write.csv(pr, file.path(out_dir, "sparse_summary.csv"), row.names = FALSE)
print(pr[pr$Model != "GBLUP", ], row.names = FALSE)
