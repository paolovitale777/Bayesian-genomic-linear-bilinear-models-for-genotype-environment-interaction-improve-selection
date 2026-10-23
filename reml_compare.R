# Compares the REML FA(1) and FA(2) models fitted with ASReml-R (reml/results/) with
# the Bayesian models, on identical testing sets, using the same metrics as gsreg_decision_gain.R
# and gsreg_sparse.R: selection differential of the top 10% (SD units), precision and correlation,
# against the true genetic values in simulations and the phenotypes otherwise.
#
#   Rscript reml_compare.R
# Writes reml_comparison.csv (CV1/CV2) and reml_sparse_comparison.csv (sparse testing).

out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))
TOP <- 0.10
sets <- c("Wheat_599", "S1_bilinear_J4", "S3_rnorm_J4", "Japonica", "EYT_16_17", "EYT_22_23")
sparse_sets <- c("Wheat_599", "S1_bilinear_J4", "S3_rnorm_J4")

# metrics of one testing set in one environment, given predictions mu and values v
metrics <- function(mu, v) {
  nt <- length(v); k <- ceiling(TOP * nt)
  sel <- rank(-mu, ties.method = "first") <= k
  top <- rank(-v, ties.method = "first") <= k
  c(n_test = nt, S = mean(v[sel]) - mean(v), Precision = mean(top[sel]), r = cor(mu, v))
}

paired <- function(x, ref, key, what) {      # paired differences of every model against `ref`
  b <- x[x$Model == ref, c(key, what)]
  out <- list()
  for (m in setdiff(unique(x$Model), ref)) {
    z <- merge(x[x$Model == m, c(key, what)], b, by = key, suffixes = c("", "_ref"))
    for (g in split(z, z[, intersect(c("Data", "CV", "k"), key)], drop = TRUE)) {
      d <- g$S - g$S_ref
      out[[length(out) + 1]] <- data.frame(g[1, intersect(c("Data", "CV", "k"), key), drop = FALSE], Model = m,
        vs = ref, n_pairs = nrow(g), dS = round(mean(d), 3), rel_pct = round(100 * mean(d) / mean(g$S_ref), 1),
        p_S = signif(t.test(d)$p.value, 2), wins = round(mean(d > 0), 2), dr = round(mean(g$r - g$r_ref), 4),
        p_r = signif(t.test(g$r - g$r_ref)$p.value, 2))
    }
  }
  do.call(rbind, out)
}

dec <- list(); spa <- list()
for (nm in sets) {
  cs <- load_case(nm); Y <- cs$Y
  V <- if (is.null(cs$truth)) Y else cs$truth
  re <- read.csv(file.path(out_dir, "reml/results", paste0("reml_fa_", nm, ".csv")), colClasses = c(line = "character", env = "character"))
  re$i <- match(re$line, rownames(Y)); re$j <- match(re$env, colnames(Y))
  stopifnot(!anyNA(re$i), !anyNA(re$j))
  re$v <- V[cbind(re$i, re$j)]
  re <- re[!is.na(Y[cbind(re$i, re$j)]), ]

  # REML metrics per testing set x environment
  g <- split(re, list(re$scheme, re$partition, re$fold, re$model, re$env), drop = TRUE)
  rm_ <- do.call(rbind, lapply(g, function(z) {
    if (nrow(z) < 20) return(NULL)
    data.frame(Data = nm, scheme = z$scheme[1], Rep = z$partition[1], Fold = z$fold[1], Env = z$env[1],
               Model = z$model[1], t(metrics(z$blup, z$v)))
  }))

  # Bayesian and BGLR models, CV1/CV2: selection differential from the selected lines, r from the folds file
  li <- read.csv(file.path(out_dir, paste0("decision_", nm, "_lines.csv")), colClasses = c(Line = "character"))
  li$v <- V[cbind(match(li$Line, rownames(Y)), match(li$Env, colnames(Y)))]
  a <- aggregate(v ~ Data + CV + Rep + Fold + Env + Model, li[li$Selected, ], mean)
  b <- aggregate(v ~ Data + CV + Rep + Fold + Env + Model, li, mean)
  a$S <- a$v - b$v[match(do.call(paste, a[, 1:6]), do.call(paste, b[, 1:6]))]
  fo <- read.csv(file.path(out_dir, paste0("decision_", nm, "_folds.csv")))
  fo$r <- if (is.null(cs$truth)) fo$r_obs else fo$r_true
  a$r <- fo$r[match(do.call(paste, a[, 1:6]), do.call(paste, fo[, c("Data", "CV", "Rep", "Fold", "Env", "Model")]))]
  rd <- rm_[rm_$scheme %in% c("CV1", "CV2"), ]; rd$CV <- rd$scheme
  dec[[nm]] <- rbind(a[, c("Data", "CV", "Rep", "Fold", "Env", "Model", "S", "r")],
                     rd[, c("Data", "CV", "Rep", "Fold", "Env", "Model", "S", "r")])

  if (nm %in% sparse_sets) {
    sp <- read.csv(file.path(out_dir, paste0("sparse_", nm, ".csv")))
    sp$S <- sp$Sel_diff
    rs <- rm_[grepl("^sparse_k", rm_$scheme), ]; rs$k <- as.integer(sub("sparse_k", "", rs$scheme))
    spa[[nm]] <- rbind(sp[, c("Data", "k", "Rep", "Env", "Model", "S", "r")], rs[, c("Data", "k", "Rep", "Env", "Model", "S", "r")])
  }
}
dec <- do.call(rbind, dec); spa <- do.call(rbind, spa)

cat("Mean selection differential and correlation, CV1/CV2:\n")
m <- aggregate(cbind(S, r) ~ Data + CV + Model, dec, mean)
print(m[order(m$Data, m$CV, -m$S), ], digits = 3, row.names = FALSE)
k1 <- c("Data", "CV", "Rep", "Fold", "Env")
pd <- rbind(paired(dec, "Reaction norm", k1, c("S", "r")), paired(dec, "Bayes genomic FA K=2", k1, c("S", "r")))
pd <- pd[pd$Model != "GBLUP", ]
cat("\nPaired differences, CV1/CV2:\n"); print(pd[order(pd$vs, pd$Data, pd$CV, pd$Model), ], row.names = FALSE)

cat("\nMean selection differential and correlation, sparse testing:\n")
ms <- aggregate(cbind(S, r) ~ Data + k + Model, spa, mean)
print(ms[order(ms$Data, ms$k, -ms$S), ], digits = 3, row.names = FALSE)
k2 <- c("Data", "k", "Rep", "Env")
ps <- rbind(paired(spa, "Reaction norm", k2, c("S", "r")), paired(spa, "Bayes genomic FA K=2", k2, c("S", "r")))
ps <- ps[ps$Model != "GBLUP", ]
cat("\nPaired differences, sparse testing:\n"); print(ps[order(ps$vs, ps$Data, ps$k, ps$Model), ], row.names = FALSE)

write.csv(merge(m, pd, all = TRUE), file.path(out_dir, "reml_comparison.csv"), row.names = FALSE)
write.csv(merge(ms, ps, all = TRUE), file.path(out_dir, "reml_sparse_comparison.csv"), row.names = FALSE)
