# Internal grid helpers for the nonlinear (measure-transportation-based)
# one-sided PCA (not exported).

# Halton sequence in the first `dim` prime bases (radical inverses of
# start, ..., start + n - 1), as used for the grids of Hallin et al. (2021).
.halton <- function(n, dim = 1L, start = 0L) {
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
              67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113)
  if (dim > length(primes))
    stop("Halton sequences are only available up to dimension 30.", call. = FALSE)
  pos <- start:(start + n - 1)
  ans <- matrix(0, n, dim)
  for (k in seq_len(dim)) {
    b <- primes[k]
    nterms <- floor(log(max(pos, 1)) / log(b)) + 1
    for (j in 0:nterms) ans[, k] <- ans[, k] + ((pos %/% b^j) %% b) / b^(j + 1)
  }
  ans
}

# "Regular" grid of (about) m unit vectors on the sphere S^{p-1}:
#   p = 1: {+1, -1};  p = 2: equispaced angles 2 pi k / m;
#   p >= 3: Halton -> Phi^{-1} -> normalise, ceiling(m / 2) vectors and their
#   antipodes (so the grid is symmetric and has an even number of points).
.sphere_grid <- function(p, m) {
  if (p == 1L) return(matrix(c(1, -1), ncol = 1L))
  if (p == 2L) {
    a <- 2 * pi * (seq_len(m) - 1) / m
    return(cbind(cos(a), sin(a)))
  }
  w <- stats::qnorm(.halton(ceiling(m / 2), p, start = 1L))
  w <- w / sqrt(rowSums(w^2))
  rbind(w, -w)
}

# Target grid G^(n): nS directions x nR radii r_j = j / (nR + 1), direction by
# direction (row (k - 1) * nR + j is r_j u_k), followed by n0 copies of the origin.
.target_grid <- function(d, nR, nS, n0) {
  dirs <- .sphere_grid(d, nS)
  r <- seq_len(nR) / (nR + 1)
  G <- do.call(rbind, lapply(seq_len(nrow(dirs)), function(k) outer(r, dirs[k, ])))
  rbind(G, matrix(0, n0, d))
}

# Factorisation n = nR * nS + n0 with 0 <= n0 <= min(nR, nS), nS >= nR and
# (for d >= 3) nS even. If neither nR nor nS is given, nR is taken close to
# sqrt(n) / 2 so that nS grows faster than nR, as required for consistency.
.grid_sizes <- function(n, d, nR = NULL, nS = NULL) {
  ok <- function(nR, nS) {
    n0 <- n - nR * nS
    nR >= 2L && nS >= 2L && n0 >= 0L && n0 <= min(nR, nS) &&
      (d == 2L || nS %% 2L == 0L)
  }
  make <- function(nR, nS) list(nR = as.integer(nR), nS = as.integer(nS),
                                n0 = as.integer(n - nR * nS))
  if (!is.null(nR) && !is.null(nS)) {
    if (!ok(nR, nS))
      stop("nR and nS must satisfy n = nR * nS + n0 with 0 <= n0 <= min(nR, nS)",
           if (d > 2L) " and nS even (d > 2)" else "", ".", call. = FALSE)
    return(make(nR, nS))
  }
  if (!is.null(nS)) {
    nR <- n %/% nS
    if (!ok(nR, nS)) stop("No valid nR for the given nS.", call. = FALSE)
    return(make(nR, nS))
  }
  try_nR <- function(nR, keep_origin = FALSE) {
    nS <- (n - keep_origin) %/% nR
    if (d > 2L && nS %% 2L == 1L) nS <- nS - 1L
    if (ok(nR, nS)) make(nR, nS) else NULL
  }
  if (!is.null(nR)) {
    res <- try_nR(nR)
    if (is.null(res)) stop("No valid nS for the given nR.", call. = FALSE)
    return(res)
  }
  target <- max(2L, floor(sqrt(n) / 2))
  ## automatic choice: keep at least one copy of the origin, so that the centre
  ## Q(0) is the observation(s) assigned to it
  for (nR in c(target:2L, (target + 1L):max(target + 1L, floor(sqrt(n))))) {
    res <- try_nR(nR, keep_origin = TRUE)
    if (!is.null(res) && res$nS >= res$nR && res$n0 >= 1L) return(res)
  }
  stop("Could not factorise n = nR * nS + n0; please supply nR and nS.", call. = FALSE)
}
