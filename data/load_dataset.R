# load_dataset.R - read one dataset of the compact data package.
#
#   source("data/load_dataset.R")          # from the repository root
#   d <- load_dataset("Wheat_1")          # or "EYT_1", "Indica", "Wheat_599", ...
#   str(d$pheno); dim(d$G)
#
# Returns a list with
#   pheno    phenotypes (long format: Env, line, traits; Cuevas data: wide, one column per environment)
#   G        genomic relationship matrix (lines x lines)
#   K        Gaussian kernel exp(-d / median(d)) from G       (Cuevas data only)
#   A        pedigree relationship matrix                     (Cuevas data only)
#   markers  marker matrix (lines x markers), when available
#   EC       environmental covariates (rice datasets)
# Record order and line order are those of the original files.

load_dataset <- function(name, dir = "data") {
  path <- file.path(dir, c(name, "GxE", "Cuevas2025", "."), paste0(name, ".rds"))
  path <- path[file.exists(path)][1]
  if (is.na(path)) stop("dataset ", name, " not found under ", dir, "/", name)
  x <- readRDS(path)
  unlower <- function(L) {
    X <- matrix(0, L$n, L$n, dimnames = L$dimnames)
    X[lower.tri(X, diag = TRUE)] <- L$values
    X[upper.tri(X)] <- t(X)[upper.tri(X)]
    X
  }
  if (!is.null(x$markers_int)) {
    M <- x$markers_int; storage.mode(M) <- "double"
    idx <- which(is.na(M), arr.ind = TRUE); M[idx] <- x$markers_fill[idx[, 2]]
    x$markers <- M; x$markers_int <- NULL; x$markers_fill <- NULL
  }
  if (!is.null(x$G_lower)) x$G <- unlower(x$G_lower)
  if (is.null(x$G) && !is.null(x$markers) && x$crop == "wheat") x$G <- tcrossprod(x$markers) / ncol(x$markers)
  if (!is.null(x$A_lower)) {
    x$A <- unlower(x$A_lower)
    d <- outer(diag(x$G), diag(x$G), "+") - 2 * x$G
    x$K <- exp(-d / median(d))
    if (!is.null(x$dimnames_K)) dimnames(x$K) <- x$dimnames_K
    if (!is.null(x$dimnames_G)) dimnames(x$G) <- x$dimnames_G
  }
  x$G_lower <- NULL; x$A_lower <- NULL
  x
}

DATASETS_GE     <- c("EYT_1", "EYT_2", "EYT_3", "Indica", "Japonica",
                     "Wheat_1", "Wheat_2", "Wheat_3", "Wheat_5", "Wheat_6")
DATASETS_CUEVAS <- c("Wheat_599", "EYT_16_17", "EYT_22_23", "EYT_23_24")
