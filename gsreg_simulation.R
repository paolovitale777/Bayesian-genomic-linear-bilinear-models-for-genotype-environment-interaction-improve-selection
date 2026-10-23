# Simulation study for the Bayesian genomic SREG / FA models, using the real
# genomic relationship matrix G of Wheat 599 (599 lines), heritability 0.5.
#
# Scenarios (argument 1):
#   S1_bilinear_J4 : true G x E bilinear, rank 2, 4 environments with crossover
#   S2_bilinear_J8 : as S1 with 8 environments
#   S3_rnorm_J4    : true G x E follows the reaction-norm model (g + gE, no low-rank structure)
# Accuracy: r_true = correlation with the TRUE genetic values; r = with simulated phenotypes.
#
#   Rscript gsreg_simulation.R S1_bilinear_J4 2000 500

args <- commandArgs(trailingOnly = TRUE)
SCEN  <- if (length(args) >= 1) args[1] else "S1_bilinear_J4"
NITER <- if (length(args) >= 2) as.integer(args[2]) else 2000
BURN  <- if (length(args) >= 3) as.integer(args[3]) else 500
out_dir <- "."
data_dir <- "data"
source(file.path(out_dir, "gsreg_functions.R"))
source(file.path(data_dir, "load_dataset.R"))
G <- load_dataset("Wheat_599", dir = data_dir)$G
n <- nrow(G)

# ---- simulate ---------------------------------------------------------------
set.seed(switch(SCEN, S1_bilinear_J4 = 11, S2_bilinear_J8 = 22, S3_rnorm_J4 = 33))
Lg <- t(chol(G + diag(1e-6, n)))
rG <- function() as.vector(Lg %*% rnorm(n))              # draw from N(0, G)
if (grepl("bilinear", SCEN)) {
  J <- if (grepl("J8", SCEN)) 8 else 4
  c1 <- rep(0.8, J)                                                        # general adaptation
  c2 <- 0.7 * rep(c(1, -1), length.out = J) * seq(1, 0.5, length.out = J)  # crossover
  Gtrue <- cbind(rG(), rG()) %*% t(cbind(c1, c2))
} else {
  J <- 4
  g <- rG(); Gtrue <- sapply(1:J, function(j) 0.6 * g + 0.6 * rG())        # g + independent gE
}
Gtrue <- scale(Gtrue, scale = FALSE)
h2 <- 0.5
s2e <- apply(Gtrue, 2, var) * (1 - h2) / h2
Y <- Gtrue + sapply(1:J, function(j) rnorm(n, 0, sqrt(s2e[j])))
colnames(Y) <- colnames(Gtrue) <- paste0("E", 1:J)
cat("Scenario", SCEN, "- lines", n, "envs", J, "\nTrue genetic correlations between environments:\n")
print(round(cor(Gtrue), 2))

t0 <- Sys.time()
acc <- run_cv(Y, G, truth = Gtrue, niter = NITER, burn = BURN, seed = 2026, label = SCEN)
acc$Scenario <- SCEN
write.csv(acc, file.path(out_dir, paste0("sim_", SCEN, "_by_fold.csv")), row.names = FALSE)
s <- aggregate(cbind(r_true, r) ~ Scenario + CV + Model, acc, mean)
s <- s[order(s$CV, -s$r_true), ]; s[, c("r_true", "r")] <- round(s[, c("r_true", "r")], 3)
write.csv(s, file.path(out_dir, paste0("sim_", SCEN, "_summary.csv")), row.names = FALSE)
print(s, row.names = FALSE)
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
