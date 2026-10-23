# Bayesian version of the crossover G x E analysis of Crossa, Yang & Cornelius (2004, JABES 9:362-380):
# instead of clustering sites with SREG/SHMM and testing perfect genetic correlation with REML, the
# posterior draws of the Bayesian genomic FA (K = 2) give directly, with their uncertainty:
#   1. genetic correlations between environments, r_jj', and P(r_jj' > RHO);
#   2. the split of the G x E variance into heterogeneity of genetic SDs (no crossover) and imperfect
#      genetic correlation (crossover): for each pair of environments
#        (s_j - s_j')^2 + 2 s_j s_j' (1 - r_jj')     (Cockerham 1963; Muir et al. 1992)
#   3. groups of environments without crossover: complete-linkage clusters of 1 - E(r), keeping the
#      coarsest partition whose groups have P(all within-group r > RHO) >= 0.9;
#   4. P(crossover) for pairs of the TOP best lines (tested and untested): the probability that their
#      ranking differs between two environments, averaged within and between groups, and for each
#      untested line against the best tested line.
# 10% of the lines are treated as never phenotyped (the same 60 lines as the Wheat 599 biplot).
# Genetic values are the posterior draws of mu_j + sum_k a_ik c_jk + d_ij, centred within environment.
#
#   Rscript gsreg_crossover.R Wheat_599 6000 1000        # name, iterations, burn-in
# Writes crossover_<name>_rg.csv, crossover_<name>_summary.csv, crossover_<name>_pairs.csv and
# crossover_<name>_untested.csv.

args <- commandArgs(trailingOnly = TRUE)
NAME  <- args[1]
NITER <- if (length(args) >= 2) as.integer(args[2]) else 6000
BURN  <- if (length(args) >= 3) as.integer(args[3]) else 1000
THIN  <- 5; RHO <- 0.8; TOP <- 20
out_dir <- "."
source(file.path(out_dir, "gsreg_functions.R"))

cs <- load_case(NAME); Y <- cs$Y; G <- cs$G
n <- nrow(Y); J <- ncol(Y); envs <- sub("_599", "", colnames(Y))
set.seed(599)
untested <- sort(sample(n, round(0.1 * n)))
Ytr <- Y; Ytr[untested, ] <- NA
t0 <- Sys.time()
fit <- gsreg(Ytr, G, K = 2, specific = TRUE, niter = NITER, burn = BURN, keep_draws = TRUE, thin = THIN)
cat(NAME, "- fit done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
D <- attr(fit, "draws"); S <- dim(D)[3]
for (s in 1:S) D[, , s] <- scale(D[, , s], scale = FALSE)          # centre within environment

# ---------------------------------------------------------------- 1-2. correlations and G x E decomposition
pairs <- which(upper.tri(diag(J)), arr.ind = TRUE)
R <- array(NA, c(J, J, S)); H <- L <- numeric(S)
for (s in 1:S) {
  Sg <- cov(D[, , s]); sd <- sqrt(diag(Sg)); R[, , s] <- cov2cor(Sg)
  h <- (sd[pairs[, 1]] - sd[pairs[, 2]])^2
  l <- 2 * sd[pairs[, 1]] * sd[pairs[, 2]] * (1 - R[, , s][pairs])
  H[s] <- sum(h); L[s] <- sum(l)
}
rg <- data.frame(Data = NAME, env1 = envs[pairs[, 1]], env2 = envs[pairs[, 2]],
                 mean = apply(R, 1:2, mean)[pairs], q2.5 = apply(R, 1:2, quantile, 0.025)[pairs],
                 q97.5 = apply(R, 1:2, quantile, 0.975)[pairs],
                 P_pos = apply(R > 0, 1:2, mean)[pairs], P_high = apply(R > RHO, 1:2, mean)[pairs])
write.csv(rg, file.path(out_dir, paste0("crossover_", NAME, "_rg.csv")), row.names = FALSE)
coi <- L / (H + L)
cat("\nGenetic correlations between environments (posterior mean, 95% CI, P(r > 0), P(r >", RHO, ")):\n")
print(rg[, -1], digits = 2, row.names = FALSE)
cat(sprintf("\nShare of G x E due to imperfect correlation (crossover): %.2f (95%% CI %.2f-%.2f)\n",
            mean(coi), quantile(coi, 0.025), quantile(coi, 0.975)))

# ---------------------------------------------------------------- 3. groups of environments without crossover
Rm <- apply(R, 1:2, mean); dimnames(Rm) <- list(envs, envs)
hc <- hclust(as.dist(1 - Rm), method = "complete")
p_group <- function(idx) if (length(idx) < 2) NA else
  mean(apply(R[idx, idx, , drop = FALSE], 3, function(r) all(r[upper.tri(r)] > RHO)))
best <- NULL
for (k in 1:J) {                                 # coarsest partition with every group P >= 0.9
  cl <- cutree(hc, k)
  pg <- sapply(1:k, function(g) p_group(which(cl == g)))
  if (all(is.na(pg) | pg >= 0.9)) { best <- list(k = k, cl = cl, pg = pg); break }
}
groups <- sapply(1:best$k, function(g) paste(envs[best$cl == g], collapse = "+"))
cat("\nGroups of environments without crossover (P(all within-group r >", RHO, ")):\n")
for (g in 1:best$k) cat(sprintf("  %-20s %s\n", groups[g], ifelse(is.na(best$pg[g]), "single", sprintf("%.2f", best$pg[g]))))

# ---------------------------------------------------------------- 4. P(crossover) among the best lines
gm <- rowMeans(apply(D, 1:2, mean))                       # posterior mean genetic value across environments
top <- order(-gm)[1:TOP]
lp <- t(combn(top, 2))
pc <- matrix(NA, nrow(lp), nrow(pairs))                    # P(order of the two lines differs between j and j')
for (p in 1:nrow(lp)) {
  dif <- D[lp[p, 1], , ] - D[lp[p, 2], , ]                 # J x S
  for (q in 1:nrow(pairs)) pc[p, q] <- mean(sign(dif[pairs[q, 1], ]) != sign(dif[pairs[q, 2], ]))
}
same <- best$cl[pairs[, 1]] == best$cl[pairs[, 2]]
pair_tab <- data.frame(Data = NAME, env1 = envs[pairs[, 1]], env2 = envs[pairs[, 2]], same_group = same,
                       mean_P_crossover = colMeans(pc), r_mean = rg$mean)
write.csv(pair_tab, file.path(out_dir, paste0("crossover_", NAME, "_pairs.csv")), row.names = FALSE)
cat("\nMean P(crossover) among the", TOP, "best lines, by pair of environments:\n")
print(pair_tab[, -1], digits = 2, row.names = FALSE)

# any crossover across all environments, for each untested line in the top against the best tested line
best_tested <- top[!(top %in% untested)][1]
ut <- data.frame(Data = NAME, line = rownames(Y)[untested], rank = match(untested, order(-gm)),
                 P_better_mean = sapply(untested, function(i) mean(colMeans(D[i, , ] - D[best_tested, , ]) > 0)),
                 P_crossover_any = sapply(untested, function(i) {
                   dif <- sign(D[i, , ] - D[best_tested, , ])        # J x S
                   mean(apply(dif, 2, function(x) length(unique(x)) > 1))
                 }),
                 P_crossover_between_groups = sapply(untested, function(i) {
                   dif <- D[i, , ] - D[best_tested, , ]
                   gmeans <- apply(dif, 2, function(x) tapply(x, best$cl, mean))   # groups x S
                   if (best$k < 2) return(NA)
                   mean(apply(sign(gmeans), 2, function(x) length(unique(x)) > 1))
                 }))
ut <- ut[order(ut$rank), ]
write.csv(ut, file.path(out_dir, paste0("crossover_", NAME, "_untested.csv")), row.names = FALSE)
cat("\nUntested lines among the best", TOP, "vs. the best tested line", rownames(Y)[best_tested], ":\n")
print(ut[ut$rank <= TOP, -1], digits = 2, row.names = FALSE)

summ <- data.frame(Data = NAME, n = n, J = J, n_untested = length(untested), draws = S,
                   coi_share = mean(coi), coi_q2.5 = quantile(coi, 0.025), coi_q97.5 = quantile(coi, 0.975),
                   groups = paste(groups, collapse = " | "),
                   P_groups = paste(ifelse(is.na(best$pg), "-", sprintf("%.2f", best$pg)), collapse = " | "),
                   Pcross_within = mean(pc[, same]), Pcross_between = if (any(!same)) mean(pc[, !same]) else NA,
                   best_tested = rownames(Y)[best_tested])
write.csv(summ, file.path(out_dir, paste0("crossover_", NAME, "_summary.csv")), row.names = FALSE)
cat(sprintf("\nMean P(crossover) among the %d best lines: within groups %.2f, between groups %.2f\n",
            TOP, summ$Pcross_within, summ$Pcross_between))
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
