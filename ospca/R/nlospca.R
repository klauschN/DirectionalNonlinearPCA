#' Nonlinear one-sided (directional) principal component analysis
#'
#' Nonlinear one-sided principal component analysis based on measure
#' transportation. In the linear one-sided analysis (see [ospca()]) a direction
#' \eqn{u} describes the half-line \eqn{\{r u: r > 0\}}; here, the half-line is
#' replaced by its image under the center-outward quantile function, the
#' so-called sign curve \eqn{r \mapsto Q_\pm(r u)}, \eqn{0 < r < 1}. The
#' empirical center-outward quantile function \eqn{Q_\pm^{(n)}} is obtained by
#' an optimal assignment of the observations to a regular grid in the unit ball
#' and is extended to a cyclically monotone interpolation
#' \eqn{\overline{Q}_\pm^{(n)}}. The principal directions maximize the expected
#' squared length of the (regularized) sign curves (see [nlos_objective()]), the
#' principal curves are the images of the corresponding radii, and the principal
#' components are curvilinear projections onto these curves.
#'
#' The procedure consists of the following steps:
#' 1. Target grid \eqn{\mathfrak{G}^{(n)}}: `nS` directions (equispaced for
#'    \eqn{d = 2}, based on a Halton sequence for \eqn{d > 2}) times `nR` radii
#'    \eqn{j / (n_R + 1)}, together with \eqn{n_0 = n - n_R n_S} copies of the
#'    origin.
#' 2. Optimal assignment of the sample to the grid (by default with the exact
#'    network simplex algorithm). The \eqn{n_0} observations assigned to the
#'    origin are merged into their mean, which keeps the coupling cyclically
#'    monotone, and cyclically monotone potentials with maximal margin are
#'    computed by Karp's minimum-mean-cycle algorithm.
#' 3. Interpolation \eqn{\overline{Q}_\pm^{(n)}}, either piecewise constant
#'    (`kappa = 0`) or its Yosida (Moreau envelope) regularization
#'    (`kappa > 0`; `"auto"` takes the largest value for which the
#'    interpolation is still exact at the grid points). The center is
#'    \eqn{\overline{Q}_\pm^{(n)}(0)}.
#' 4. Discrete sign curves at the radii \eqn{r_j = j/(m_R + 1)}, smoothed by the
#'    Tikhonov regularization \eqn{(I + \varepsilon L'L)^{-1}}, where `eps` is
#'    used for the criterion, and thus for the choice of the directions, and
#'    `eps_curves` for the returned principal curves, the scores and
#'    [reconstruct()].
#' 5. The direction \eqn{U_1} maximizes the criterion over a grid of `mS`
#'    directions. Each further direction \eqn{U_i}, \eqn{i \ge 2}, maximizes it
#'    over a grid on the unit sphere of the orthogonal complement of the
#'    previous directions, together with the antipodes of previous directions
#'    that have not been selected yet. All \eqn{2d} directions are computed, and
#'    the `n_components` with the largest criterion are kept, in decreasing
#'    order. As the empirical criterion is rough, the direction grid should be
#'    dense. By default, it has at least \eqn{n} directions, so that
#'    \eqn{m = m_R m_S \ge n}, and it contains the directions of the target
#'    grid (for \eqn{d = 2} since `mS` is a multiple of `nS`, and for
#'    \eqn{d > 2} since both grids are based on the same Halton sequence).
#' 6. Scores: column \eqn{i} contains the nonlinear principal component
#'    \eqn{P_{\pm,i}^{(n)}(X)}, that is, the sum of the positive parts of the
#'    projections of the increments of the observation's own sign curve onto
#'    the unit tangents of the \eqn{i}th principal curve. [reconstruct()]
#'    projects onto the manifold spanned by several principal curves.
#'
#' With `paths = "smooth"` (default), the observation's own sign curve is
#' regularized with the same Tikhonov filter as the principal curves before its
#' increments are projected. With the raw empirical sign curve
#' (`paths = "raw"`), the positive part of each increment makes the zigzags of
#' the piecewise constant interpolation count as progress, so that observations
#' opposite to a principal direction receive positive scores. Note that the
#' population definition, with the positive part taken at each step, is the
#' same in both cases.
#'
#' A warning is given if a direction is a near-tie, that is, if its criterion is
#' within the relative tolerance `tol_tie` of the best candidate more than 45
#' degrees away, or if two consecutive retained criteria differ by less than
#' `tol_tie` (relative).
#'
#' @param X Numeric data matrix with at least two columns, observations in rows.
#' @param n_components Number of directions to keep. Defaults to `2 * ncol(X)`,
#'   which is also the maximum.
#' @param nR,nS Number of radii and of directions of the target grid, with
#'   `nrow(X) = nR * nS + n0`, `0 <= n0 <= min(nR, nS)` and, for \eqn{d > 2},
#'   `nS` even. If both are `NULL`, `nR` is chosen close to \eqn{\sqrt{n}/2}
#'   with \eqn{n_0 \ge 1}; if only one is given, the other one is derived from
#'   it. With \eqn{n_0 = 0} the center is \eqn{\overline{Q}_\pm^{(n)}(0)}.
#' @param mR Number of radii of the grids on which the curves are evaluated
#'   (default `nR`).
#' @param mS Size of the direction grids searched for the principal directions.
#'   The default is the smallest multiple of `nS` that is at least \eqn{n}, so
#'   that the grid has at least as many directions as there are observations
#'   and contains the directions of the target grid (see Details).
#' @param eps Tikhonov regularization parameter \eqn{\varepsilon} of the curves
#'   whose lengths define the criterion (`0` for no smoothing).
#' @param eps_curves Tikhonov regularization parameter of the returned
#'   principal curves, of the observations' sign curves and of the projections
#'   (default `eps`). A smaller value gives, for example, lightly regularized
#'   curves along directions that were chosen with more regularization.
#' @param paths `"smooth"` (default) or `"raw"`: whether the observations' own
#'   sign curves are regularized before they are projected (see Details).
#' @param delta Trimming proportion \eqn{\delta} of the criterion.
#' @param kappa Yosida regularization parameter of the interpolation: `0`
#'   (piecewise constant, default), a positive number, or `"auto"`.
#' @param n_cone Size of the sphere grid whose points in the positive span of
#'   the retained directions are used as candidate directions in the
#'   projections (default `mS`). The target grid directions in that span and the
#'   retained directions themselves are always candidates.
#' @param ot_method Optimal assignment algorithm of [transport::transport()];
#'   `"networkflow"` (default) and `"shortsimplex"` are exact.
#' @param nthreads Number of threads of the compiled code. Defaults to the option
#'   `ospca.threads` if it is set, and otherwise to all cores but one (two under
#'   `R CMD check`); in both cases at most the value of the environment variable
#'   `OMP_THREAD_LIMIT`, if it is set.
#' @param tol_tie Relative tolerance for the near-tie warnings.
#' @param verbose Logical; if `TRUE`, the progress is reported.
#'
#' @return An object of class `nlospca`, a list with the following components,
#'   where \eqn{n} denotes the number of observations, \eqn{p} the number of
#'   variables and \eqn{k} the number of retained directions:
#'   \item{directions}{\eqn{p \times k}{p x k} matrix of the principal directions, in
#'     its columns and in decreasing order of the criterion.}
#'   \item{objective}{The criterion values of the \eqn{k} directions (see
#'     [nlos_objective()]).}
#'   \item{curves}{List of the \eqn{k} principal curves, each a matrix of
#'     points in the space of the data that starts at the center.}
#'   \item{scores}{\eqn{n \times k}{n x k} matrix of the nonlinear principal
#'     components, that is, the curvilinear projection lengths.}
#'   \item{center}{The center \eqn{\overline{Q}_\pm^{(n)}(0)}.}
#'   \item{X}{The data matrix.}
#'   \item{n, p}{Number of observations and of variables.}
#'   \item{transport}{List describing the fitted empirical center-outward
#'     quantile function: the target grid (`G`, `dirs`, `nR`, `nS`, `n0`), the
#'     optimal assignment (`Y`, `perm`, `obs_grid`), the interpolation nodes
#'     (`Gi`, `Yi`), the potentials (`psi`) with their margin (`e0`), `kappa`,
#'     `center` and `nthreads`.}
#'   \item{rays}{Array with the empirical sign curves along the directions of
#'     the target grid, used for the projections.}
#'   \item{settings}{List of the settings used (`nR`, `nS`, `n0`, `mR`, `mS`,
#'     `eps`, `eps_curves`, `paths`, `delta`, `kappa`, `n_cone`, `ot_method`,
#'     `nthreads`).}
#'   \item{diagnostics}{Data frame with, for each retained direction, the step
#'     of the greedy search in which it was selected, its criterion, the
#'     relative `gap` to the best candidate more than 45 degrees away, the
#'     near-tie flag `tie` and the number of candidate directions.}
#'   \item{call}{The matched call.}
#'
#' @references Hallin, M., Nordhausen, K. and Radojicic, U. Directional
#'   nonlinear principal components: a measure transportation approach.
#'
#'   Hallin, M., del Barrio, E., Cuesta-Albertos, J. and Matrán, C. (2021).
#'   Distribution and quantile functions, ranks and signs in dimension d: a
#'   measure transportation approach. *Annals of Statistics*, 49, 1139--1165.
#'   \doi{10.1214/20-AOS1996}
#'
#' @seealso [nlos_objective()], [nlos_curve()], [reconstruct()], [screeplot()],
#'   and [ospca()] for the linear version.
#' @examples
#' set.seed(1)
#' x <- rnorm(400)
#' X <- cbind(x, x^2 + rnorm(400, sd = 0.5))
#' fit <- nlospca(X)
#' fit
#' plot(fit)
#' screeplot(fit)
#' @export
nlospca <- function(X, n_components = 2L * ncol(X), nR = NULL, nS = NULL,
                   mR = NULL, mS = NULL, eps = 0.5, eps_curves = eps,
                   paths = c("smooth", "raw"), delta = 0.1, kappa = 0,
                   n_cone = NULL,
                   ot_method = c("networkflow", "shortsimplex", "auctionbf"),
                   nthreads = NULL, tol_tie = 1e-3, verbose = FALSE) {
  ot_method <- match.arg(ot_method)
  paths <- match.arg(paths)
  nthreads <- if (is.null(nthreads)) .default_threads() else max(1L, as.integer(nthreads))
  X <- .check_matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (p < 2L) stop("nlospca requires at least two variables.", call. = FALSE)
  n_components <- as.integer(n_components)
  if (length(n_components) != 1L || n_components < 1L)
    stop("n_components must be a positive integer.", call. = FALSE)
  n_components <- min(n_components, 2L * p)

  sz <- .grid_sizes(n, p, nR, nS)
  mR <- as.integer(if (is.null(mR)) sz$nR else mR)
  ## default: the smallest multiple of nS that is at least n, so that the
  ## direction grid is dense and contains the directions of the target grid
  mS <- as.integer(if (is.null(mS)) sz$nS * ceiling(n / sz$nS) else mS)
  n_cone <- as.integer(if (is.null(n_cone)) mS else n_cone)
  if (mR < 2L || mS < 2L || n_cone < 2L)
    stop("mR, mS and n_cone must be at least 2.", call. = FALSE)
  for (e in list(eps, eps_curves))
    if (length(e) != 1L || !is.finite(e) || e < 0)
      stop("eps and eps_curves must be single non-negative numbers.", call. = FALSE)
  if (length(delta) != 1L || !is.finite(delta) || delta < 0 || delta >= 1)
    stop("delta must be in [0, 1).", call. = FALSE)

  tf <- .ot_fit(X, sz$nR, sz$nS, sz$n0, kappa = kappa, method = ot_method,
                nthreads = nthreads, verbose = verbose)
  if (verbose) message("Principal directions ...")
  pd <- .nlos_directions(tf, 2L * p, mR, mS, eps, delta, tol_tie, verbose)
  ord <- order(pd$objective, decreasing = TRUE)[seq_len(n_components)]

  U <- t(pd$directions[ord, , drop = FALSE])
  objective <- pd$objective[ord]
  labels <- paste0("NLOSPC", seq_len(ncol(U)))
  colnames(U) <- labels
  rownames(U) <- colnames(X)
  curves <- lapply(seq_len(ncol(U)), function(j) {
    C <- rbind(tf$center, .nlos_curve_points(tf, U[, j], mR, eps_curves)[[1L]])
    colnames(C) <- colnames(X)
    C
  })
  names(curves) <- labels
  diag <- cbind(component = labels, pd$diagnostics[ord, , drop = FALSE])
  rownames(diag) <- NULL

  fit <- structure(
    list(directions = U, objective = objective, curves = curves, scores = NULL,
         center = tf$center, X = X, n = n, p = p, transport = tf,
         rays = .nlos_rays(tf),
         settings = list(nR = sz$nR, nS = sz$nS, n0 = sz$n0, mR = mR, mS = mS,
                         eps = eps, eps_curves = eps_curves, paths = paths,
                         delta = delta, kappa = tf$kappa, n_cone = n_cone,
                         ot_method = ot_method, nthreads = nthreads),
         diagnostics = diag, call = match.call()),
    class = "nlospca")
  if (verbose) message("Nonlinear principal components ...")
  Z <- vapply(seq_len(ncol(U)), function(j) .nlos_project(fit, t(U[, j]))$score,
              numeric(n))
  Z <- matrix(Z, nrow = n, dimnames = list(rownames(X), labels))
  fit$scores <- Z

  ## ---- tie / multiplicity warnings ----
  msgs <- character(0)
  if (any(diag$tie, na.rm = TRUE))
    msgs <- c(msgs, sprintf(
      "direction(s) %s are near-tied (relative criterion gap < tol_tie = %g): the maximiser may be non-unique; results can be unstable.",
      paste(diag$component[which(diag$tie)], collapse = ", "), tol_tie))
  if (length(objective) >= 2L) {
    rel <- abs(diff(objective)) / pmax(objective[-length(objective)], .Machine$double.eps)
    if (any(rel < tol_tie)) {
      jj <- which(rel < tol_tie)
      msgs <- c(msgs, sprintf(
        "retained directions %s have near-equal criterion (relative gap < %g): their ordering is essentially arbitrary.",
        paste(sprintf("%d/%d", jj, jj + 1L), collapse = ", "), tol_tie))
    }
  }
  if (length(msgs)) warning(paste(msgs, collapse = "\n"), call. = FALSE)
  fit
}

# Greedy search of the 2d nonlinear principal directions (selection order).
.nlos_directions <- function(tf, K, mR, mS, eps, delta, tol_tie = 1e-3,
                            verbose = FALSE, tol = 1e-8) {
  d <- ncol(tf$Y)
  md <- .trim_index(mR, delta)
  sel <- matrix(numeric(0), 0L, d)
  obj <- numeric(0)
  diag <- data.frame(step = integer(0), objective = numeric(0), gap = numeric(0),
                     tie = logical(0), candidates = integer(0))
  for (i in seq_len(K)) {
    if (i == 1L) {
      cand <- .sphere_grid(d, mS)
    } else {
      ## antipodes of selected directions that are not selected themselves
      neg <- -sel
      new <- apply(neg, 1L, function(v) all(colSums((t(sel) - v)^2) > tol))
      cand2 <- neg[new, , drop = FALSE]
      ## unit sphere of the orthogonal complement of the selected directions
      sv <- svd(t(sel), nu = d)
      rk <- sum(sv$d > tol)
      cand1 <- if (rk < d) {
        .sphere_grid(d - rk, mS) %*% t(sv$u[, (rk + 1L):d, drop = FALSE])
      } else matrix(numeric(0), 0L, d)
      cand <- rbind(cand1, cand2)
      if (nrow(cand) == 0L) break
    }
    M <- vapply(.nlos_curve_points(tf, cand, mR, eps), .nlos_crit, numeric(1), md = md)
    b <- which.max(M)
    far <- drop(cand %*% cand[b, ]) < cos(pi / 4)
    gap <- if (any(far)) max(0, (M[b] - max(M[far])) / max(M[b], .Machine$double.eps))
           else NA_real_
    sel <- rbind(sel, cand[b, ])
    obj <- c(obj, M[b])
    diag <- rbind(diag, data.frame(step = i, objective = M[b], gap = gap,
                                   tie = isTRUE(gap < tol_tie),
                                   candidates = nrow(cand)))
    if (verbose) message(sprintf("  step %d: criterion %.4g (%d candidates)",
                                 i, M[b], nrow(cand)))
  }
  list(directions = sel, objective = obj, diagnostics = diag)
}

# Empirical sign curves of the target grid: array nS x (nR + 1) x d with
# [k, 1, ] = centre and [k, l + 1, ] = Q_n(r_l u_k).
.nlos_rays <- function(tf) {
  d <- ncol(tf$Y); nR <- tf$nR; nS <- tf$nS
  R <- array(0, c(nS, nR + 1L, d))
  for (j in seq_len(d)) {
    R[, 1L, j] <- tf$center[j]
    R[, -1L, j] <- matrix(tf$Y[seq_len(nR * nS), j], nS, nR, byrow = TRUE)
  }
  R
}

# Candidate directions for the positive span of the directions in the columns of
# U (each pair orthogonal or antipodal): the directions themselves, the
# target-grid directions `dirs` lying in the span (with them, the reconstruction
# with all 2d directions and unsmoothed curves is exact), and a grid of n_cone
# points on the sphere of span(U) folded onto the cone (coordinates on axes
# without their antipode made non-negative); needed for d >= 3, where a cone of
# lower dimension contains (almost) no target-grid directions.
.cone_candidates <- function(U, n_cone, dirs = NULL, tol = 1e-6) {
  U <- t(U)
  axes <- matrix(numeric(0), 0L, ncol(U)); free <- logical(0)
  for (i in seq_len(nrow(U))) {
    s <- U[i, ]
    m <- if (nrow(axes)) which(abs(drop(axes %*% s) + 1) < tol) else integer(0)
    if (length(m)) { free[m[1L]] <- TRUE; next }
    axes <- rbind(axes, s); free <- c(free, FALSE)
  }
  W <- .sphere_grid(nrow(axes), n_cone)
  W[, !free] <- abs(W[, !free])
  V <- W %*% axes
  V <- V / sqrt(rowSums(V^2))
  if (!is.null(dirs)) {
    Cf <- dirs %*% t(axes)                      # coordinates on the axes
    inside <- rowSums(dirs^2) - rowSums(Cf^2) < tol &
      apply(Cf >= -tol | rep(free, each = nrow(Cf)), 1L, all)
    V <- rbind(dirs[inside, , drop = FALSE], V)
  }
  V <- rbind(U, V)
  V[!duplicated(round(V, 8)), , drop = FALSE]
}

# Point at arc length s (from the first point) along the polygon P; clamped to
# the end points.
.point_at_length <- function(P, s) {
  seg <- sqrt(rowSums(diff(P)^2))
  cl <- c(0, cumsum(seg))
  s <- pmin(pmax(s, 0), cl[length(cl)])
  m <- findInterval(s, cl, rightmost.closed = TRUE, all.inside = TRUE)
  w <- ifelse(seg[m] > 0, (s - cl[m]) / seg[m], 0)
  P[m, , drop = FALSE] + w * (P[m + 1L, , drop = FALSE] - P[m, , drop = FALSE])
}

# Nonlinear projection onto the manifold Qbar(cone spanned by the rows of V).
# An observation x with F(x) = r_l u_k has the sign curve C = S Q_k of its ray,
# Q_k = (Q_n(r_1 u_k), ..., Q_n(r_nR u_k)) with Q_n(r_l u_k) replaced by x, and
# S the Tikhonov filter with parameter eps (paths = "smooth") or the identity
# (paths = "raw"); C(r_0) = centre. Its score along the curve of direction v is
#   score(v) = sum_{j <= l} < t_j(v), C(r_j) - C(r_{j-1}) >^+,
# t_j(v) the unit tangents of the curve along v (same eps) at the target radii.
# The projection is the point at arc length max_v score(v), measured from the
# centre, on the curve of the maximising v; observations assigned to the origin
# project onto the centre. Replacing Q_n(r_l u_k) by x changes the increments by
# dS[j, l] (x - Q_n(r_l u_k)) with dS = diff(rbind(0, S)), so training
# observations (x = Q_n(r_l u_k)) use precomputed cumulative sums.
.nlos_project <- function(object, V, X = NULL, eps = object$settings$eps_curves,
                         paths = object$settings$paths) {
  tf <- object$transport; st <- object$settings
  nR <- st$nR; nS <- st$nS; d <- object$p
  if (is.null(X)) {
    X <- object$X
    g <- tf$obs_grid
  } else {
    g <- .fbar_index(X, tf$Gi, tf$Yi, tf$psi)
  }
  n <- nrow(X)
  orig <- g > nR * nS
  k <- (g - 1L) %/% nR + 1L
  l <- (g - 1L) %% nR + 1L
  k[orig] <- 1L; l[orig] <- 1L

  ## observations' sign curves along the target rays (smoothed or raw)
  S <- if (paths == "smooth" && eps > 0) .tikhonov(nR, eps) else diag(nR)
  rays <- object$rays
  for (h in seq_len(d)) rays[, -1L, h] <- rays[, -1L, h] %*% t(S)
  Delta <- rays[, -1L, , drop = FALSE] - rays[, -(nR + 1L), , drop = FALSE]

  ## new observations: deviation from the grid point they replace
  Yg <- matrix(object$rays[cbind(rep(k, d), rep(l + 1L, d), rep(seq_len(d), each = n))], n, d)
  E <- X - Yg
  E[orig, ] <- 0
  gen <- which(rowSums(abs(E)) > 0)
  if (length(gen)) {
    dS <- diff(rbind(0, S))
    dSl <- t(dS[, l[gen], drop = FALSE])                 # dS[j, l_i]
    upto <- outer(l[gen], seq_len(nR), ">=")             # j <= l_i
  }

  curves <- lapply(.nlos_curve_points(tf, V, nR, eps), function(C) rbind(tf$center, C))
  sc <- matrix(0, n, length(curves))
  for (j in seq_along(curves)) {
    Dc <- diff(curves[[j]])
    Tn <- Dc / pmax(sqrt(rowSums(Dc^2)), .Machine$double.xmin)
    IP <- matrix(0, nS, nR)
    for (h in seq_len(d)) IP <- IP + Delta[, , h] * rep(Tn[, h], each = nS)
    CS <- cbind(0, t(apply(pmax(IP, 0), 1L, cumsum)))
    s <- CS[cbind(k, l + 1L)]
    if (length(gen)) {
      M <- IP[k[gen], , drop = FALSE] + dSl * (E[gen, , drop = FALSE] %*% t(Tn))
      s[gen] <- rowSums(pmax(M, 0) * upto)
    }
    s[orig] <- 0
    sc[, j] <- s
  }
  best <- max.col(sc, ties.method = "first")
  score <- sc[cbind(seq_len(n), best)]
  P <- matrix(NA_real_, n, d)
  for (j in unique(best)) {
    ii <- which(best == j)
    P[ii, ] <- .point_at_length(curves[[j]], score[ii])
  }
  colnames(P) <- colnames(object$X)
  list(score = score, direction = best, points = P)
}
