# Data-scarcity experiment: REML FA(1) and FA(2) with a genomic relationship matrix, fitted with
# ASReml-R 4 on exactly the same training sets as the Bayesian models. For each data set, subsets
# of 75, 150 and 300 lines are analysed under two schemes:
#   sparse_k1: every line is phenotyped in only one environment;
#   unbal_50:  every cell is observed with probability 0.5 (at least one cell per line).
# The design (which cells are used for training) is in data/scarcity_<data set>.csv; the cells with
# train = FALSE are predicted. Same model as run_reml_fa.R:
#
#   y_ij = mu_j + g_ij + e_ij,   vec(g) ~ N(0, (Lambda Lambda' + Psi) (x) G_sub),   e_ij ~ N(0, s2_j)
#
# where G_sub is G restricted to the lines of the subset. Writes results/reml_scarcity_<data set>.csv;
# the script can be stopped and restarted (testing sets already written are skipped).
#
#   Rscript run_reml_scarcity.R                        # Wheat_599 and S1_bilinear_J4
#   Rscript run_reml_scarcity.R Wheat_599              # one data set
# Run from the reml folder (it reads data/ and writes results/).

library(asreml)
asreml.options(workspace = "1gb", pworkspace = "1gb", trace = FALSE)

args <- commandArgs(trailingOnly = TRUE)
SETS <- if (length(args)) args else c("Wheat_599", "S1_bilinear_J4")
dir.create("results", showWarnings = FALSE)

fit_fa <- function(dat, k, G) {
  fm <- asreml(fixed = y ~ env, random = ~ fa(env, k):vm(line, G),
               residual = ~ dsum(~ units | env), data = dat,
               na.action = na.method(y = "include"), maxit = 50)
  tries <- 0
  while (!isTRUE(fm$converge) && tries < 5) { fm <- update(fm); tries <- tries + 1 }
  pv <- predict(fm, classify = "env:line")$pvals
  list(pv = pv, converged = isTRUE(fm$converge))
}

for (nm in SETS) {
  pheno <- read.csv(file.path("data", paste0("pheno_", nm, ".csv")), colClasses = c("character", "character", "numeric"))
  Gall <- as.matrix(read.csv(file.path("data", paste0("G_", nm, ".csv")), row.names = 1, check.names = FALSE))
  des <- read.csv(file.path("data", paste0("scarcity_", nm, ".csv")), colClasses = c(line = "character", env = "character"))
  out_file <- file.path("results", paste0("reml_scarcity_", nm, ".csv"))
  done <- if (file.exists(out_file)) unique(with(read.csv(out_file), paste(n_sub, scheme, partition, model))) else character(0)
  t0 <- Sys.time()
  for (key in unique(paste(des$n_sub, des$scheme, des$partition))) {
    z <- des[paste(des$n_sub, des$scheme, des$partition) == key, ]
    lines <- sort(unique(z$line))
    G <- Gall[lines, lines]
    G <- G + diag(1e-4 * mean(diag(G)), nrow(G))         # keeps G positive definite
    attr(G, "INVERSE") <- FALSE
    dat <- pheno[pheno$line %in% lines, ]
    trn <- paste(z$line, z$env)[z$train]
    dat$y[!(paste(dat$line, dat$env) %in% trn)] <- NA
    dat$line <- factor(dat$line, levels = lines); dat$env <- factor(dat$env)
    test <- paste(z$line, z$env)[!z$train]
    for (k in 1:2) {
      model <- paste0("REML FA(", k, ")")
      if (paste(key, model) %in% done) next
      res <- tryCatch(fit_fa(dat, k, G), error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
      if (is.null(res)) {
        out <- data.frame(data = nm, n_sub = z$n_sub[1], scheme = z$scheme[1], partition = z$partition[1], model = model,
                          line = NA, env = NA, blup = NA, pev = NA, converged = FALSE)
      } else {
        pv <- res$pv; keep <- paste(pv$line, pv$env) %in% test
        out <- data.frame(data = nm, n_sub = z$n_sub[1], scheme = z$scheme[1], partition = z$partition[1], model = model,
                          line = pv$line[keep], env = pv$env[keep], blup = pv$predicted.value[keep],
                          pev = pv$std.error[keep]^2, converged = res$converged)
      }
      write.table(out, out_file, sep = ",", row.names = FALSE, col.names = !file.exists(out_file), append = file.exists(out_file))
    }
    cat(nm, key, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
}
