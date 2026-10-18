# Internal measure-transportation machinery for nlospca() (not exported):
# optimal assignment sample <-> target grid, cyclically monotone potentials,
# and the interpolations Qbar (ball -> data) and F (data -> grid).

# Number of threads for the compiled code: option "ospca.threads", else 2 under
# R CMD check (_R_CHECK_LIMIT_CORES_), else all cores but one; in every case at
# most OMP_THREAD_LIMIT, if that environment variable is set.
.default_threads <- function() {
  nt <- getOption("ospca.threads")
  if (!is.null(nt)) {
    nt <- max(1L, as.integer(nt))
  } else {
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    nt <- if (nzchar(chk) && chk != "false") 2L
          else max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
  }
  lim <- suppressWarnings(as.integer(Sys.getenv("OMP_THREAD_LIMIT", "")))
  if (!is.na(lim) && lim >= 1L) nt <- min(nt, lim)
  nt
}

# Squared-Euclidean optimal assignment between the rows of G and X.
# Returns perm with G[i, ] <-> X[perm[i], ].
.ot_assign <- function(G, X, method = "networkflow") {
  a <- transport::transport(transport::pp(G), transport::pp(X), p = 2,
                            method = method)
  a <- a[order(a$from), ]
  as.integer(a$to)
}

# Potentials psi for the cyclically monotone map g_i -> y_i with maximal margin
#   <g_i, y_i> - psi_i >= <g_i, y_j> - psi_j + e0   for all i != j,
# by Karp's minimum-mean-cycle algorithm (src/karp.cpp, O(N^3) time; all walk
# lengths are stored, 8 N^2 bytes, for N <= 6000, else recomputed in three passes).
.ot_potentials <- function(G, Y, nthreads = .default_threads()) {
  N <- nrow(G)
  res <- karp_potentials_cpp(G, Y, nthreads = nthreads, store = N <= 6000L)
  ## attained margin, recomputed from the potentials in floating point (this is
  ## what decides whether the interpolation is exact at the grid points)
  e0 <- .margin(G, Y, res$psi)
  ## a margin that is negative only by rounding error is a zero margin (ties)
  tol <- sqrt(.Machine$double.eps) * max(rowSums(Y^2))
  if (e0 < -tol)
    warning("negative margin of the potentials: the assignment is not optimal (use an exact method).",
            call. = FALSE)
  else if (e0 <= 0)
    warning("zero margin of the potentials (tied or duplicated observations?); the interpolation is not exact at tied grid points.",
            call. = FALSE)
  list(psi = res$psi, e0 = e0)
}

# min_i [ <g_i, y_i> - psi_i - max_{j != i} (<g_i, y_j> - psi_j) ], in blocks.
.margin <- function(G, Y, psi, chunk = 1000L) {
  N <- nrow(G); out <- Inf
  for (st in seq(1L, N, by = chunk)) {
    ii <- st:min(N, st + chunk - 1L)
    S <- G[ii, , drop = FALSE] %*% t(Y) - rep(psi, each = length(ii))
    own <- S[cbind(seq_along(ii), ii)]
    S[cbind(seq_along(ii), ii)] <- -Inf
    out <- min(out, own - S[cbind(seq_along(ii), max.col(S, ties.method = "first"))])
  }
  out
}

# Largest Yosida parameter guaranteeing Qbar(g_i) = y_i at all grid points:
# kappa * R < e0 / diam.
.kappa_auto <- function(Y, e0) {
  R <- max(sqrt(rowSums(Y^2)))
  diam <- max(stats::dist(Y))
  0.99 * e0 / (R * diam)
}

# Interpolation Qbar(z) for the rows z of Z (src/yosida.cpp).
#  kappa = 0: T0(z) = y_{argmax_j <z, y_j> - psi_j}   (piecewise constant);
#  kappa > 0: gradient of the Moreau envelope (Yosida regularisation) of
#             phi(z) = max_j (<z, y_j> - psi_j), i.e. Qbar(z) = Y' lambda with
#             lambda = argmin over the simplex of
#             (kappa / 2) ||Y' lambda||^2 - <z, Y' lambda> + psi' lambda,
#             solved by a primal active-set method.
.qbar <- function(Z, Y, psi, kappa = 0, nthreads = .default_threads()) {
  Z <- matrix(as.numeric(Z), ncol = ncol(Y))
  Q <- yosida_cpp(Z, Y, as.numeric(psi), kappa, nthreads = nthreads)
  nf <- attr(Q, "not_converged")
  if (!is.null(nf) && nf > 0)
    warning(nf, " Yosida interpolation(s) did not converge within the iteration limit.",
            call. = FALSE)
  attr(Q, "not_converged") <- NULL
  Q
}

# Piecewise-constant cyclically monotone extension of the empirical
# center-outward distribution function: the interpolation-node index assigned to
# each row of X, F(x) = g_{argmax_i <x, g_i> - phi_i} with
# phi_i = <g_i, y_i> - psi_i (the conjugate potential; it reproduces F(y_i) = g_i
# with the same margin).
.fbar_index <- function(X, G, Y, psi, chunk = 2000L) {
  phi <- rowSums(G * Y) - psi
  out <- integer(nrow(X))
  for (st in seq(1L, nrow(X), by = chunk)) {
    ii <- st:min(nrow(X), st + chunk - 1L)
    S <- X[ii, , drop = FALSE] %*% t(G) - rep(phi, each = length(ii))
    out[ii] <- max.col(S, ties.method = "first")
  }
  out
}

# Fit the empirical center-outward quantile function: assignment, potentials,
# interpolation parameter and centre Qbar(0).
#
# Interpolation nodes (Gi, Yi): the nR * nS non-zero grid points with their
# observations and, if n0 > 0, one node g = 0 carrying the mean of the n0
# observations assigned to the origin. Merging keeps the coupling cyclically
# monotone (the conditions are linear in the point attached to g = 0), whereas
# separate copies of the origin would form a zero-weight cycle and force a zero
# margin everywhere.
.ot_fit <- function(X, nR, nS, n0, kappa = 0, method = "networkflow",
                    nthreads = .default_threads(), verbose = FALSE) {
  d <- ncol(X); m <- nR * nS
  dirs <- .sphere_grid(d, nS)
  G <- .target_grid(d, nR, nS, n0)
  if (verbose) message("Optimal assignment (", method, ") ...")
  perm <- .ot_assign(G, X, method)
  Y <- X[perm, , drop = FALSE]
  Gi <- G[seq_len(m), , drop = FALSE]; Yi <- Y[seq_len(m), , drop = FALSE]
  if (n0 > 0L) {
    Gi <- rbind(Gi, 0)
    Yi <- rbind(Yi, colMeans(Y[m + seq_len(n0), , drop = FALSE]))
  }
  storage.mode(Gi) <- storage.mode(Yi) <- "double"
  if (verbose) message("Cyclically monotone potentials (Karp) ...")
  pot <- .ot_potentials(Gi, Yi, nthreads)
  if (identical(kappa, "auto")) kappa <- .kappa_auto(Yi, pot$e0)
  kappa <- as.numeric(kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0)
    stop("kappa must be \"auto\" or a single non-negative number.", call. = FALSE)
  center <- .qbar(matrix(0, 1L, d), Yi, pot$psi, kappa, nthreads)[1L, ]
  obs_grid <- integer(nrow(X)); obs_grid[perm] <- seq_along(perm)
  list(G = G, dirs = dirs, Y = Y, perm = perm, obs_grid = obs_grid,
       Gi = Gi, Yi = Yi, psi = pot$psi, e0 = pot$e0, kappa = kappa,
       center = center, nR = nR, nS = nS, n0 = n0, nthreads = nthreads)
}
