# Bayesian genomic SREG and genomic factor-analytic (FA) models, with an optional
# enviromic prior on the environment scores, and the cross-validation used to
# compare them with GBLUP and reaction-norm models.
# Sourced by gsreg_wheat599.R, gsreg_simulation.R and gsreg_run.R.
#
# Genomic SREG (rank K), genotype main effect absorbed in the bilinear term:
#   y_ij = mu_j + sum_k a_ik c_jk + e_ij,        e_ij ~ N(0, s2_j)
#   a_k ~ N(0, G),  c_k ~ N(0, t2_k Omega),  t2_k ~ IG(at, bt)
#   Omega = I (no covariables) or W W'/q + lambda I (environmental covariables W)
# Genomic FA (specific = TRUE) adds an environment-specific genomic term:
#   + d_ij,  d_j ~ N(0, psi_j G),  psi_j ~ IG(ap, bp)
# so that the genetic covariance is (C C' + Psi) (x) G, the FA_K structure of
# Smith, Cullis & Thompson / Tolhurst et al., fitted here by Gibbs sampling.
# Genomic AMMI (main = TRUE) adds a genomic main effect g ~ N(0, s2g G), so the
# bilinear term models only the interaction; in a new environment (CV0) the main
# effect carries the prediction and the covariable-driven scores add the specific part.
#
# Sampler: missing cells are imputed at every iteration (data augmentation), so the
# data are complete and every n x n full conditional is diagonal in the eigenbasis of
# G, computed once: each draw costs O(n^2) instead of a Cholesky factorisation.
# Predictions use mu_j + sum_k a_ik c_jk (+ d_ij), invariant to sign and rotation.

suppressMessages(library(BGLR))

gsreg <- function(Y, Kmat, K = 2, specific = FALSE, main = FALSE, Omega = NULL, niter = 3000, burn = 1000,
                  at = 3, bt = 1, as0 = 1, bs0 = 0.5, ap = 3, bp = 0.2, ag = 3, bg = 0.5, eigK = NULL,
                  keep_draws = FALSE, thin = 10, keep_trace = FALSE) {
  n <- nrow(Y); J <- ncol(Y)
  O <- !is.na(Y); nobs <- colSums(O)
  if (is.null(eigK)) eigK <- eigen(Kmat, symmetric = TRUE)
  U <- eigK$vectors; lam <- pmax(eigK$values, 0) + 1e-4 * mean(diag(Kmat))
  Oinv <- if (is.null(Omega)) diag(J) else solve(Omega)
  mu <- ifelse(nobs > 0, colMeans(Y, na.rm = TRUE), 0)
  s2 <- rep(0.5, J); t2 <- rep(0.5, K); psi <- rep(0.1, J)
  A <- matrix(rnorm(n * K, 0, 0.1), n, K); C <- matrix(rnorm(J * K, 0, 0.1), J, K)
  D <- matrix(0, n, J); g <- rep(0, n); s2g <- 0.3
  Yc <- Y; Yc[!O] <- matrix(mu, n, J, byrow = TRUE)[!O]
  # draw x ~ N(P^-1 b, P^-1) with P = U diag(1/(f lam) + s) U'
  draw_g <- function(b, s, f = 1) {
    q <- 1 / (f * lam) + s
    as.vector(U %*% ((crossprod(U, b) + sqrt(q) * rnorm(n)) / q))
  }
  pred <- matrix(0, n, J); nkeep <- 0
  draws <- if (keep_draws) array(NA_real_, c(n, J, floor((niter - burn) / thin))) else NULL
  s2draws <- if (keep_draws) matrix(NA_real_, floor((niter - burn) / thin), J) else NULL
  # trace of scalar parameters and of the identifiable genetic correlations between environments,
  # from the J x J genetic covariance s2g 11' + C C' + Psi (all terms share G, which cancels)
  pairs <- which(upper.tri(diag(J)), arr.ind = TRUE)
  trace <- if (keep_trace) list() else NULL
  for (it in 1:niter) {
    M <- matrix(mu, n, J, byrow = TRUE)
    if (main) {                                  # genomic main effect (AMMI form): g ~ N(0, s2g G)
      R <- Yc - M - D - A %*% t(C)
      g <- draw_g(as.vector(R %*% (1 / s2)), sum(1 / s2), s2g)
      s2g <- 1 / rgamma(1, ag + n / 2, bg + sum(crossprod(U, g)^2 / lam) / 2)
    }
    M <- M + g
    for (k in 1:K) {
      R <- Yc - M - D - A[, -k, drop = FALSE] %*% t(C[, -k, drop = FALSE])
      A[, k] <- draw_g(as.vector(R %*% (C[, k] / s2)), sum(C[, k]^2 / s2))
      Pc <- Oinv / t2[k] + diag(sum(A[, k]^2) / s2, J)          # joint draw of c_k
      Uc <- chol(Pc)
      bc <- as.vector(crossprod(R, A[, k])) / s2
      C[, k] <- backsolve(Uc, forwardsolve(t(Uc), bc)) + backsolve(Uc, rnorm(J))
      t2[k] <- 1 / rgamma(1, at + J / 2, bt + sum(C[, k] * (Oinv %*% C[, k])) / 2)
    }
    AC <- A %*% t(C)
    if (specific) {
      for (j in 1:J) {
        D[, j] <- draw_g((Yc[, j] - mu[j] - g - AC[, j]) / s2[j], 1 / s2[j], psi[j])
        z <- crossprod(U, D[, j])
        psi[j] <- 1 / rgamma(1, ap + n / 2, bp + sum(z^2 / lam) / 2)
      }
    }
    AC <- AC + g                                 # genotype part of the fit (main effect + bilinear)
    for (j in 1:J) {
      if (nobs[j] > 0) {
        oj <- O[, j]
        mu[j] <- rnorm(1, mean(Y[oj, j] - AC[oj, j] - D[oj, j]), sqrt(s2[j] / nobs[j]))
        e <- Y[oj, j] - mu[j] - AC[oj, j] - D[oj, j]
        s2[j] <- 1 / rgamma(1, as0 + nobs[j] / 2, bs0 + sum(e^2) / 2)
      } else {                                   # environment with no data (CV0): mean not identified
        mu[j] <- 0; s2[j] <- mean(s2[nobs > 0])
      }
    }
    Fit <- matrix(mu, n, J, byrow = TRUE) + AC + D
    Yc[!O] <- Fit[!O] + rnorm(sum(!O)) * rep(sqrt(s2), each = n)[!O]   # data augmentation
    if (it > burn) {
      pred <- pred + Fit; nkeep <- nkeep + 1
      if (keep_draws && (it - burn) %% thin == 0) {
        draws[, , (it - burn) / thin] <- Fit            # genotype values per draw
        s2draws[(it - burn) / thin, ] <- s2             # residual variances, for posterior predictive draws
      }
      if (keep_trace && (it - burn) %% thin == 0) {
        Sg <- tcrossprod(C) + (if (main) s2g else 0) + diag(if (specific) psi else 0, J)
        Rg <- cov2cor(Sg)
        trace[[length(trace) + 1]] <- c(setNames(s2, paste0("sigma2_", 1:J)), setNames(t2, paste0("tau2_", 1:K)),
          if (main) c(sigma2_u = s2g), if (specific) setNames(psi, paste0("psi_", 1:J)),
          setNames(Rg[pairs], paste0("rg_", pairs[, 1], pairs[, 2])))
      }
    }
  }
  out <- pred / nkeep
  if (keep_draws) { attr(out, "draws") <- draws; attr(out, "s2") <- s2draws }
  if (keep_trace) attr(out, "trace") <- do.call(rbind, trace)
  out
}

# Environmental similarity from covariables W (J x q): Omega = W W'/q + lambda I.
make_omega <- function(W, lambda = 0.1) {
  W <- scale(as.matrix(W)); W[is.na(W)] <- 0
  tcrossprod(W) / ncol(W) + diag(lambda, nrow(W))
}

# Cross-validation of all models on the same folds.
#   Y: lines x environments (NA allowed), G: genomic relationship matrix,
#   Omega: optional J x J environmental similarity (enables CV0 and the "+EC" models),
#   truth: optional matrix of true genetic values (simulation).
run_cv <- function(Y, G, Omega = NULL, truth = NULL, cvs = c("CV1", "CV2"), niter = 3000, burn = 1000,
                   seed = 2026, label = "", fa = TRUE) {
  n <- nrow(Y); J <- ncol(Y); envs <- colnames(Y)
  eigG <- eigen(G, symmetric = TRUE)
  Zg <- model.matrix(~ factor(rep(1:n, J)) - 1); Ze <- model.matrix(~ factor(rep(1:J, each = n)) - 1)
  KG <- Zg %*% G %*% t(Zg)
  eig <- function(K) { e <- eigen(K, symmetric = TRUE); keep <- e$values > 1e-8; list(V = e$vectors[, keep], d = e$values[keep]) }
  EG <- eig(KG); EGE <- eig(KG * tcrossprod(Ze))
  if (!is.null(Omega)) { KW <- Ze %*% Omega %*% t(Ze); EW <- eig(KW); EGW <- eig(KG * KW); rm(KW) }
  rm(KG); invisible(gc())
  fit_bglr <- function(yNA, type) {
    ETA <- list(E = list(X = Ze[, -1, drop = FALSE], model = "FIXED"), G = list(V = EG$V, d = EG$d, model = "RKHS"))
    if (type == "RN")   ETA$GE <- list(V = EGE$V, d = EGE$d, model = "RKHS")
    if (type == "RNEC") { ETA$W <- list(V = EW$V, d = EW$d, model = "RKHS"); ETA$GW <- list(V = EGW$V, d = EGW$d, model = "RKHS") }
    fm <- BGLR(y = yNA, ETA = ETA, nIter = niter * 2, burnIn = burn, verbose = FALSE,
               saveAt = file.path(tempdir(), paste0(label, "_", type, "_")))
    matrix(fm$yHat, n, J)
  }
  bglr_models <- c("GBLUP (G+E)" = "G", "Reaction norm (G+E+GxE)" = "RN")
  mk <- function(K, M = G, sp = FALSE, main = FALSE, Om = NULL) list(K = K, M = M, sp = sp, main = main, Om = Om)
  models <- list("GSREG K=1" = mk(1), "GSREG K=2" = mk(2),
                 "GAMMI K=1" = mk(1, main = TRUE), "GAMMI K=2" = mk(2, main = TRUE))
  if (fa) models <- c(models, list("Bayes genomic FA K=1" = mk(1, sp = TRUE), "Bayes genomic FA K=2" = mk(2, sp = TRUE)))
  models[["SREG K=2 (no markers)"]] <- mk(2, M = diag(n))
  if (!is.null(Omega)) {
    bglr_models <- c(bglr_models, "Reaction norm + EC (G+W+GxW)" = "RNEC")
    models[["GAMMI K=1 + EC"]] <- mk(1, main = TRUE, Om = Omega)
    models[["GAMMI K=2 + EC"]] <- mk(2, main = TRUE, Om = Omega)
  }
  set.seed(seed)
  fold_line <- sample(rep(1:5, length.out = n))                    # CV1
  fold_cell <- matrix(sample(rep(1:5, length.out = n * J)), n, J)  # CV2
  acc <- list()
  add <- function(cv, model, fold, P, test) {
    for (j in 1:J) {
      t <- test[, j] & !is.na(Y[, j])
      if (sum(t) > 5) acc[[length(acc) + 1]] <<- data.frame(CV = cv, Model = model, Fold = fold, Env = envs[j],
        n_test = sum(t), r = cor(P[t, j], Y[t, j]),
        r_true = if (is.null(truth)) NA else cor(P[t, j], truth[t, j]))
    }
  }
  # a model is skipped where it has no information: CV1 needs markers; CV0 needs covariables
  usable <- function(cv, m, mm = NULL) {
    if (cv == "CV1" && grepl("no markers", m)) return(FALSE)
    if (cv == "CV0") return(grepl("EC", m) || m == "GBLUP (G+E)")
    TRUE
  }
  t0 <- Sys.time()
  for (cv in cvs) {
    nf <- if (cv == "CV0") J else 5
    for (f in 1:nf) {
      test <- switch(cv, CV1 = matrix(fold_line == f, n, J), CV2 = (fold_cell == f),
                     CV0 = matrix(rep(1:J == f, each = n), n, J))
      Ytr <- Y; Ytr[test] <- NA
      if (cv == "CV0" && sum(!is.na(Y[, f])) < 10) next
      for (b in names(bglr_models)) if (usable(cv, b)) add(cv, b, f, fit_bglr(as.vector(Ytr), bglr_models[[b]]), test)
      for (m in names(models)) {
        if (!usable(cv, m)) next
        mm <- models[[m]]
        add(cv, m, f, gsreg(Ytr, mm$M, K = mm$K, specific = mm$sp, main = mm$main, Omega = mm$Om, niter = niter, burn = burn,
                            eigK = if (identical(mm$M, G)) eigG else NULL), test)
      }
      cat(label, cv, "fold", f, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
    }
  }
  do.call(rbind, acc)
}

# The four models of the decision and sparse-testing analyses, as functions of the training
# matrix Ytr. Each returns the predicted means and two arrays of draws (lines x envs x draws):
# genetic values (draws) and phenotypes (pred = posterior predictive, with residual error).
# BGLR keeps no draws of the predictions, so for GBLUP and the reaction norm the draws are
# independent normals N(yHat, SD.yHat) (+ varE for pred), a plug-in approximation.
make_fitters <- function(Y, G, niter = 2000, burn = 500, thin = 5, label = "fit") {
  n <- nrow(Y); J <- ncol(Y)
  eigG <- eigen(G, symmetric = TRUE)
  Zg <- model.matrix(~ factor(rep(1:n, J)) - 1); Ze <- model.matrix(~ factor(rep(1:J, each = n)) - 1)
  KG <- Zg %*% G %*% t(Zg)
  eig <- function(K) { e <- eigen(K, symmetric = TRUE); keep <- e$values > 1e-8; list(V = e$vectors[, keep], d = e$values[keep]) }
  EG <- eig(KG); EGE <- eig(KG * tcrossprod(Ze)); rm(KG); invisible(gc())
  S <- floor((niter - burn) / thin)
  fit_bglr <- function(Ytr, gxe) {
    ETA <- list(E = list(X = Ze[, -1, drop = FALSE], model = "FIXED"), G = list(V = EG$V, d = EG$d, model = "RKHS"))
    if (gxe) ETA$GE <- list(V = EGE$V, d = EGE$d, model = "RKHS")
    fm <- BGLR(y = as.vector(Ytr), ETA = ETA, nIter = niter * 2, burnIn = burn, verbose = FALSE,
               saveAt = file.path(tempdir(), paste0(label, "_", gxe, "_")))
    mu <- matrix(fm$yHat, n, J); sdv <- matrix(fm$SD.yHat, n, J)
    list(mean = mu, draws = array(rnorm(n * J * S, mu, sdv), c(n, J, S)),
         pred = array(rnorm(n * J * S, mu, sqrt(sdv^2 + fm$varE)), c(n, J, S)))
  }
  fit_bayes <- function(Ytr, specific) {
    P <- gsreg(Ytr, G, K = 2, specific = specific, main = !specific, niter = niter, burn = burn,
               eigK = eigG, keep_draws = TRUE, thin = thin)
    D <- attr(P, "draws"); s2 <- attr(P, "s2")
    E <- array(rnorm(n * J * S), c(n, J, S)) * aperm(array(sqrt(s2), c(S, J, n)), c(3, 2, 1))
    list(mean = as.matrix(P[, ]), draws = D, pred = D + E)
  }
  list("GBLUP"                = function(Ytr) fit_bglr(Ytr, FALSE),
       "Reaction norm"        = function(Ytr) fit_bglr(Ytr, TRUE),
       "GAMMI K=2"            = function(Ytr) fit_bayes(Ytr, FALSE),
       "Bayes genomic FA K=2" = function(Ytr) fit_bayes(Ytr, TRUE))
}

DATA_DIR <- "data"

# Simulated phenotypes with the real G (heritability h2 in every environment).
#   S1/S2: bilinear G x E of rank 2 (4 or 8 environments, crossover); S3: reaction-norm G x E.
simulate_case <- function(G, scen, h2 = 0.5) {
  n <- nrow(G)
  set.seed(switch(scen, S1_bilinear_J4 = 11, S2_bilinear_J8 = 22, S3_rnorm_J4 = 33))
  Lg <- t(chol(G + diag(1e-6, n)))
  rG <- function() as.vector(Lg %*% rnorm(n))              # draw from N(0, G)
  if (grepl("bilinear", scen)) {
    J <- if (grepl("J8", scen)) 8 else 4
    c1 <- rep(0.8, J)                                                        # general adaptation
    c2 <- 0.7 * rep(c(1, -1), length.out = J) * seq(1, 0.5, length.out = J)  # crossover
    Gtrue <- cbind(rG(), rG()) %*% t(cbind(c1, c2))
  } else {
    J <- 4
    g <- rG(); Gtrue <- sapply(1:J, function(j) 0.6 * g + 0.6 * rG())        # g + independent gE
  }
  Gtrue <- scale(Gtrue, scale = FALSE)
  s2e <- apply(Gtrue, 2, var) * (1 - h2) / h2
  Y <- Gtrue + sapply(1:J, function(j) rnorm(n, 0, sqrt(s2e[j])))
  colnames(Y) <- colnames(Gtrue) <- paste0("E", 1:J)
  list(Y = Y, truth = Gtrue)
}

# One analysis case: real data set of the record-order repository or a simulation scenario.
# Returns Y (lines x environments, standardised within environment for real data), G,
# Omega (rice only) and truth (simulation only).
load_case <- function(name) {
  source(file.path(DATA_DIR, "load_dataset.R"), local = TRUE)
  if (grepl("^S[0-9]_", name)) {
    G <- load_dataset("Wheat_599", dir = DATA_DIR)$G
    s <- simulate_case(G, name)
    rownames(s$Y) <- rownames(s$truth) <- paste0("L", seq_len(nrow(s$Y)))
    return(list(Y = s$Y, G = G, Omega = NULL, truth = s$truth))
  }
  d <- load_dataset(name, dir = DATA_DIR); Omega <- NULL
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
  } else {
    envs <- if (name == "Wheat_599") sort(d$env_names) else d$env_names
    Y <- as.matrix(d$pheno[, envs]); rownames(Y) <- d$pheno$Line
    G <- d$G
  }
  Y <- scale(Y)
  keep <- rowSums(!is.na(Y)) > 0
  list(Y = Y[keep, , drop = FALSE], G = G[keep, keep], Omega = Omega, truth = NULL)
}
