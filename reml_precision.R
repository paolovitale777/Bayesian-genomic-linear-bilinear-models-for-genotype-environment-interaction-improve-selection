# Precision of selecting the top 10% (fraction of the selected lines that belong to the true top 10%)
# for the REML FA models , compared with the models of gsreg_decision.R on the same
# testing sets. True genetic values in simulations, phenotypes otherwise.
#
#   Rscript reml_precision.R

out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))
res <- list()
for (nm in c("Wheat_599", "S1_bilinear_J4", "S3_rnorm_J4", "Japonica", "EYT_16_17", "EYT_22_23")) {
  cs <- load_case(nm); Y <- cs$Y; V <- if (is.null(cs$truth)) Y else cs$truth
  re <- read.csv(file.path(out_dir, "reml/results", paste0("reml_fa_", nm, ".csv")), colClasses = c(line = "character", env = "character"))
  re <- re[re$scheme %in% c("CV1", "CV2"), ]
  re$v <- V[cbind(match(re$line, rownames(Y)), match(re$env, colnames(Y)))]
  re <- re[!is.na(re$v), ]
  pr <- do.call(rbind, lapply(split(re, list(re$scheme, re$partition, re$fold, re$model, re$env), drop = TRUE), function(z) {
    nt <- nrow(z); if (nt < 20) return(NULL); k <- ceiling(0.1 * nt)
    data.frame(CV = z$scheme[1], Rep = z$partition[1], Fold = z$fold[1], Env = z$env[1], Model = z$model[1],
               Precision = mean((rank(-z$v, ties.method = "first") <= k)[rank(-z$blup, ties.method = "first") <= k]))
  }))
  fo <- read.csv(file.path(out_dir, paste0("decision_", nm, "_folds.csv")))
  fo$Precision <- if (is.null(cs$truth)) fo$Precision_obs else fo$Precision_true
  x <- rbind(fo[, c("CV", "Rep", "Fold", "Env", "Model", "Precision")], pr)
  x$Data <- nm
  res[[nm]] <- x
}
x <- do.call(rbind, res)
m <- aggregate(Precision ~ Data + CV + Model, x, mean)
w <- reshape(m, idvar = c("Data", "CV"), timevar = "Model", direction = "wide")
names(w) <- sub("Precision.", "", names(w))
print(w[, c("Data", "CV", "GBLUP", "Reaction norm", "REML FA(1)", "REML FA(2)", "GAMMI K=2", "Bayes genomic FA K=2")],
      digits = 3, row.names = FALSE)
key <- c("Data", "CV", "Rep", "Fold", "Env")
a <- merge(x[x$Model == "REML FA(2)", c(key, "Precision")], x[x$Model == "GAMMI K=2", c(key, "Precision")], by = key, suffixes = c("", "_G"))
b <- merge(x[x$Model == "REML FA(2)", c(key, "Precision")], x[x$Model == "Bayes genomic FA K=2", c(key, "Precision")], by = key, suffixes = c("", "_F"))
cat("\nREML FA(2) minus GAMMI and minus Bayes FA (paired):\n")
for (g in split(a, a[, c("Data", "CV")], drop = TRUE)) {
  h <- b[b$Data == g$Data[1] & b$CV == g$CV[1], ]
  cat(sprintf("%-15s %s  vs GAMMI %+.3f (p = %.2g)   vs Bayes FA %+.3f (p = %.2g)\n", g$Data[1], g$CV[1],
      mean(g$Precision - g$Precision_G), t.test(g$Precision - g$Precision_G)$p.value,
      mean(h$Precision - h$Precision_F), t.test(h$Precision - h$Precision_F)$p.value))
}
