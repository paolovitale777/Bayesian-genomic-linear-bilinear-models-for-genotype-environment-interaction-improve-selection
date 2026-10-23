# Runs the model comparison on one of the data sets of the record-order repository.
#
#   Rscript gsreg_run.R EYT_16_17 3000 1000
#   Rscript gsreg_run.R Japonica 3000 1000          # rice: with environmental covariables and CV0
#
# Phenotypes are standardised within environment. For the rice sets, G is computed
# from the SNPs (centred markers, scaled to mean diagonal 1) and Omega from the 54
# environmental covariables of each year. Writes <name>_by_fold.csv and <name>_summary.csv.

args <- commandArgs(trailingOnly = TRUE)
NAME  <- args[1]
NITER <- if (length(args) >= 2) as.integer(args[2]) else 3000
BURN  <- if (length(args) >= 3) as.integer(args[3]) else 1000
out_dir <- "."
data_dir <- "data"
source(file.path(out_dir, "gsreg_functions.R"))
source(file.path(data_dir, "load_dataset.R"))
d <- load_dataset(NAME, dir = data_dir)

Omega <- NULL; cvs <- c("CV1", "CV2")
if (d$crop == "rice") {
  ph <- d$pheno; ph$Env <- as.character(ph$Env)
  envs <- sort(unique(ph$Env)); lines <- rownames(d$markers)
  if (is.null(lines)) lines <- d$markers_row_names
  Y <- matrix(NA, length(lines), length(envs), dimnames = list(lines, envs))
  Y[cbind(match(ph$Line, lines), match(ph$Env, envs))] <- ph$GY
  M <- scale(d$markers, scale = FALSE)
  G <- tcrossprod(M); G <- G / mean(diag(G)); dimnames(G) <- list(lines, lines)
  W <- d$EC[match(envs, as.character(d$EC$Year)), setdiff(names(d$EC), "Year")]
  Omega <- make_omega(W); dimnames(Omega) <- list(envs, envs)
  cvs <- c("CV1", "CV2", "CV0")
  cat("Environmental similarity (Omega) between years:\n"); print(round(cov2cor(Omega), 2))
} else {
  envs <- d$env_names
  Y <- as.matrix(d$pheno[, envs]); rownames(Y) <- d$pheno$Line
  G <- d$G
  if (!is.null(d$A)) stopifnot(all(as.character(d$pheno$Line) == rownames(d$A)))
}
Y <- scale(Y)                                            # standardise within environment
keep <- rowSums(!is.na(Y)) > 0; Y <- Y[keep, , drop = FALSE]; G <- G[keep, keep]
cat(NAME, "- lines:", nrow(Y), " environments:", ncol(Y), " observed cells:", sum(!is.na(Y)), "\n")
cat("Phenotypic correlations between environments:\n"); print(round(cor(Y, use = "pairwise"), 2))

t0 <- Sys.time()
acc <- run_cv(Y, G, Omega = Omega, cvs = cvs, niter = NITER, burn = BURN, seed = 2026, label = NAME)
acc$Data <- NAME
write.csv(acc, file.path(out_dir, paste0(NAME, "_by_fold.csv")), row.names = FALSE)

s <- aggregate(r ~ Data + CV + Model, acc, mean)
s$sd_folds <- aggregate(r ~ Data + CV + Model, acc, sd)$r
s <- s[order(s$CV, -s$r), ]; s[, c("r", "sd_folds")] <- round(s[, c("r", "sd_folds")], 3)
write.csv(s, file.path(out_dir, paste0(NAME, "_summary.csv")), row.names = FALSE)
print(s, row.names = FALSE)
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
