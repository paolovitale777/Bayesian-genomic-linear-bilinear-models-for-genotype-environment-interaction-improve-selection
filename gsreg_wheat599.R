# Bayesian genomic SREG and genomic FA vs. GBLUP and reaction norm - Wheat 599
# (599 lines, 4 environments, genomic relationship matrix G).
#
#   Rscript gsreg_wheat599.R 3000 1000      # iterations, burn-in
#
# Writes accuracy_by_fold.csv and accuracy_summary.csv (mean r per environment).

args <- commandArgs(trailingOnly = TRUE)
NITER <- if (length(args) >= 1) as.integer(args[1]) else 3000
BURN  <- if (length(args) >= 2) as.integer(args[2]) else 1000
out_dir <- "."
data_dir <- "data"
source(file.path(out_dir, "gsreg_functions.R"))
source(file.path(data_dir, "load_dataset.R"))

d <- load_dataset("Wheat_599", dir = data_dir)
envs <- sort(d$env_names)                                  # E1..E4
Y <- as.matrix(d$pheno[, envs]); rownames(Y) <- d$pheno$Line
stopifnot(all(as.character(d$pheno$Line) == rownames(d$A)))   # G, A and phenotypes share the line order
cat("Lines:", nrow(Y), " Environments:", ncol(Y), "\n")

t0 <- Sys.time()
acc <- run_cv(Y, d$G, niter = NITER, burn = BURN, seed = 2026, label = "W599")
write.csv(acc, file.path(out_dir, "accuracy_by_fold.csv"), row.names = FALSE)

s <- aggregate(r ~ CV + Model + Env, acc, mean)
tab <- reshape(s, idvar = c("CV", "Model"), timevar = "Env", direction = "wide")
names(tab) <- sub("^r\\.", "", names(tab))
tab$Mean <- rowMeans(tab[, envs])
tab <- tab[order(tab$CV, -tab$Mean), ]
tab[, c(envs, "Mean")] <- round(tab[, c(envs, "Mean")], 3)
write.csv(tab, file.path(out_dir, "accuracy_summary.csv"), row.names = FALSE)
print(tab, row.names = FALSE)
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
