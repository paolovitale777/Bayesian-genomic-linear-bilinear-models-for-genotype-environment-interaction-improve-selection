# Predicted GGE biplot with 95% credible regions for UNTESTED lines, and the
# posterior probability that each line wins (or is in the top 10%) in each
# environment - Wheat 599, genomic AMMI (GAMMI) with K = 2.
#
# 60 lines (10%) are treated as never phenotyped: they enter only through G.
# Each posterior draw of the genotype values is environment-centred and
# decomposed by SVD (symmetric scaling); draws are aligned to the posterior-mean
# biplot by orthogonal Procrustes before the credible ellipses are computed.
#
#   Rscript gsreg_biplot_wheat599.R 3000 1000

args <- commandArgs(trailingOnly = TRUE)
NITER <- if (length(args) >= 1) as.integer(args[1]) else 3000
BURN  <- if (length(args) >= 2) as.integer(args[2]) else 1000
out_dir <- "."
data_dir <- "data"
source(file.path(out_dir, "gsreg_functions.R"))
source(file.path(data_dir, "load_dataset.R"))

d <- load_dataset("Wheat_599", dir = data_dir)
envs <- sort(d$env_names); Y <- as.matrix(d$pheno[, envs]); rownames(Y) <- d$pheno$Line
colnames(Y) <- sub("_599", "", envs); envs <- colnames(Y)
n <- nrow(Y); J <- ncol(Y)

set.seed(599)
untested <- sort(sample(n, 60))
Ytr <- Y; Ytr[untested, ] <- NA
fit_file <- file.path(out_dir, "wheat599_gammi_fit.rds")
if (file.exists(fit_file)) fit <- readRDS(fit_file) else {
  fit <- gsreg(Ytr, d$G, K = 2, main = TRUE, niter = NITER, burn = BURN, keep_draws = TRUE, thin = 10)
  saveRDS(fit, fit_file)
}
draws <- attr(fit, "draws"); S <- dim(draws)[3]

# validation: predicted vs observed (hidden) values of the untested lines
r_untested <- sapply(1:J, function(j) cor(fit[untested, j], Y[untested, j]))
names(r_untested) <- envs
cat("Accuracy for the 60 untested lines (r with their hidden phenotypes):\n"); print(round(r_untested, 3))

# ---- biplot coordinates per draw, aligned to the posterior mean ----------------------
bip <- function(X) {                        # genotype-focused scaling: g = U D, e = V
  X <- scale(X, scale = FALSE); s <- svd(X)
  list(g = s$u[, 1:2] %*% diag(s$d[1:2]), e = s$v[, 1:2],
       pct = 100 * s$d[1:2]^2 / sum(s$d^2))
}
ref <- bip(fit)
Gd <- array(NA, c(n, 2, S)); Ed <- array(NA, c(J, 2, S))
for (s in 1:S) {
  b <- bip(draws[, , s])
  sv <- svd(crossprod(b$g, ref$g)); Q <- sv$u %*% t(sv$v)          # orthogonal Procrustes
  Gd[, , s] <- b$g %*% Q; Ed[, , s] <- b$e %*% Q
}

# ---- probability of winning / top 10% per environment -----------------------------------
top <- ceiling(0.10 * n)
p_win <- sapply(1:J, function(j) tabulate(apply(draws[, j, ], 2, which.max), n) / S)
p_top <- sapply(1:J, function(j) rowMeans(apply(draws[, j, ], 2, function(v) rank(-v) <= top)))
dimnames(p_win) <- dimnames(p_top) <- list(rownames(Y), envs)
status <- ifelse(seq_len(n) %in% untested, "untested", "tested")
res <- data.frame(Line = rownames(Y), status, round(fit, 3), round(p_win, 3), round(p_top, 3))
names(res) <- c("Line", "Status", paste0("pred_", envs), paste0("Pwin_", envs), paste0("Ptop10_", envs))
write.csv(res, file.path(out_dir, "wheat599_biplot_probabilities.csv"), row.names = FALSE)
for (j in 1:J) {
  o <- order(-p_top[, j])[1:8]
  cat("\n", envs[j], "- lines most likely to be in the top 10%:\n")
  print(data.frame(Line = rownames(Y)[o], Status = status[o], P_top10 = round(p_top[o, j], 2), P_win = round(p_win[o, j], 2)), row.names = FALSE)
}

# ---- figure ---------------------------------------------------------------------------------------
ell <- function(P, level = 0.95) {          # 2-D credible ellipse from aligned draws
  m <- colMeans(P); V <- cov(P); r <- sqrt(qchisq(level, 2))
  t <- seq(0, 2 * pi, length.out = 100); e <- eigen(V)
  sweep(cbind(cos(t), sin(t)) %*% diag(r * sqrt(pmax(e$values, 0))) %*% t(e$vectors), 2, m, "+")
}
gm <- apply(Gd, c(1, 2), mean); em <- apply(Ed, c(1, 2), mean)
em <- em * 0.9 * max(abs(gm)) / max(abs(em))                                 # environment vectors rescaled for display
best <- apply(p_top, 1, max)                                                  # best chance of the top 10% in any environment
show_u <- untested[order(-best[untested])][1:6]
show_t <- setdiff(seq_len(n), untested)[order(-best[-untested])][1:4]
png(file.path(out_dir, "wheat599_predicted_biplot.png"), width = 1800, height = 1500, res = 200)
par(mar = c(4.5, 4.5, 3, 1))
lim <- range(c(gm, em * 1.2))
plot(gm, type = "n", xlim = lim, ylim = lim, asp = 1,
     xlab = sprintf("Component 1 (%.0f%%)", ref$pct[1]), ylab = sprintf("Component 2 (%.0f%%)", ref$pct[2]),
     main = "Wheat 599 - GAMMI (K = 2): predicted GGE biplot", cex.main = 0.95)
mtext("Labelled: the 6 untested and 4 tested lines with the highest posterior probability of the top 10%", cex = 0.7)
abline(h = 0, v = 0, col = "grey85")
points(gm[-untested, ], pch = 16, cex = 0.5, col = "grey65")
points(gm[untested, ], pch = 17, cex = 0.8, col = "#D55E00")
for (i in show_t) polygon(ell(t(Gd[i, , ])), border = "grey35", lwd = 1.2)
for (i in show_u) { polygon(ell(t(Gd[i, , ])), border = "#D55E00", lwd = 1.4); text(gm[i, 1], gm[i, 2], rownames(Y)[i], pos = 3, cex = 0.6, col = "#D55E00") }
arrows(0, 0, em[, 1], em[, 2], length = 0.08, col = "#0072B2", lwd = 2)
text(em * 1.1, labels = envs, col = "#0072B2", font = 2)
legend("bottomleft", bty = "n", cex = 0.8,
       legend = c("Tested lines", "Untested lines (markers only)", "95% credible region, untested", "95% credible region, tested", "Environments"),
       pch = c(16, 17, NA, NA, NA), lty = c(NA, NA, 1, 1, 1), lwd = c(NA, NA, 1.4, 1.2, 2),
       col = c("grey65", "#D55E00", "#D55E00", "grey35", "#0072B2"))
dev.off()
cat("\nFigure: wheat599_predicted_biplot.png; table: wheat599_biplot_probabilities.csv\n")
