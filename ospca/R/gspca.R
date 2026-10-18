#' Independent nonlinear component analysis (GS-PCA)
#'
#' A simple discrete implementation of the independent nonlinear component
#' analysis of Gunsilius and Schennach (2023), referred to here as GS-PCA. Let
#' \eqn{X} be a random vector in \eqn{R^d} with density \eqn{f}, and let
#' \eqn{T: R^d \to R^d} denote the Brenier map that pushes the distribution of
#' \eqn{X} forward to the \eqn{d}-variate standard normal distribution, whose
#' density is denoted by \eqn{\varphi}. Points in the space of the data are
#' denoted by \eqn{y} and points in the space of the normal reference
#' distribution by \eqn{x = T(y)}. For a map \eqn{g}, \eqn{J_g} denotes its
#' Jacobian matrix, the matrix of its first-order partial derivatives; as
#' \eqn{T} and its inverse \eqn{T^{-1}} are gradients of convex functions,
#' \eqn{J_T(y)} and \eqn{J_{T^{-1}}(x)} are symmetric, and
#' \eqn{J_{T^{-1}}(x) = J_T(T^{-1}(x))^{-1}}. The directions
#' \eqn{u_1, \ldots, u_d} are the eigenvectors of
#' \deqn{\bar J = -\int f(y) \ln J_T(y)\, dy = \int \varphi(x) \ln J_{T^{-1}}(x)\, dx,}
#' in decreasing order of the eigenvalues, where \eqn{\ln} denotes the matrix
#' logarithm (the logarithm of the eigenvalues in the spectral decomposition)
#' and the second equality follows from the change of variables
#' \eqn{x = T(y)}. The effective contribution of a unit vector \eqn{u} to the
#' entropy of \eqn{X} is \eqn{\frac{1}{2}\ln(2\pi e) + u^\top \bar J u}{(1/2) ln(2 pi e) + u' Jbar u}
#' (Gunsilius and Schennach, 2023, Lemma 1), so that the eigenvalue of
#' \eqn{u_j} is its contribution minus the entropy of a univariate standard
#' normal, and the first \eqn{k} directions maximize the total contribution of
#' \eqn{k} orthonormal directions (Theorem 2 of the same paper). The
#' principal curves are \eqn{T^{-1}(r u_j)}, \eqn{r \in R}, the scores of an
#' observation \eqn{X} are the coordinates \eqn{\langle T(X), u_j \rangle},
#' \eqn{j = 1, \ldots, k}, and its \eqn{k}-dimensional representation is
#' \eqn{T^{-1}(\sum_{j \le k} \langle T(X), u_j \rangle u_j)} (see
#' [reconstruct()]). For Gaussian data, \eqn{T} is linear and the directions
#' are the linear principal components (Gunsilius and Schennach, 2023,
#' Corollary 1). With covariance matrix \eqn{\Sigma},
#' \eqn{J_{T^{-1}} = \Sigma^{1/2}}{J_{T^{-1}} = Sigma^(1/2)}, so that the
#' eigenvalues of \eqn{\bar J} are half the logarithms of the eigenvalues of
#' \eqn{\Sigma}.
#'
#' The estimation consists of the following steps:
#' 1. Reference grid: the tensor product of the \eqn{n_g} = `n_grid` standard
#'    normal quantiles \eqn{\Phi^{-1}((i - 1/2) / n_g)},
#'    \eqn{i = 1, \ldots, n_g}, on each axis, where \eqn{\Phi} denotes the
#'    standard normal distribution function. This gives \eqn{m = n_g^d \le n}
#'    equally weighted grid points.
#' 2. Discrete Brenier map: the exact optimal coupling between the grid (mass
#'    \eqn{1/m} per grid point) and the sample \eqn{X_1, \ldots, X_n} (mass
#'    \eqn{1/n} per observation). \eqn{T^{-1}} at a grid point is the
#'    barycenter of the observations it is coupled with; this is a selection of
#'    the subdifferential of the convex potential, and hence cyclically
#'    monotone. \eqn{T(X_i)} is the (mass-weighted) grid point the observation
#'    \eqn{X_i} is coupled with.
#' 3. Jacobians of \eqn{T^{-1}} at the grid points by central differences over
#'    \eqn{\pm}`step` grid steps along each axis, symmetrized as in Theorem 4 of
#'    Gunsilius and Schennach (2023), since the Jacobian of a gradient map is
#'    symmetric. Their matrix logarithms are averaged over the grid, which
#'    approximates the expectation under the standard normal distribution. Grid
#'    points within `step` of the edge of the grid are treated according to
#'    `boundary`. Grid points whose Jacobian estimate is not positive definite
#'    (because of noise or ties) are left out of the average; their proportion
#'    is reported by `print()` and `summary()`, with a warning when it exceeds
#'    5 percent.
#' 4. \eqn{T^{-1}} at arbitrary points (for the curves and the reconstructions)
#'    by multilinear interpolation on the grid, clamped to the grid box, and
#'    \eqn{T} at new observations ([predict()]) by the piecewise constant
#'    cyclically monotone extension, i.e., a grid point. The latter needs the
#'    conjugate potentials, which are computed by Karp's algorithm when needed,
#'    in \eqn{O(m^3)} time for \eqn{m} grid points.
#'
#' Note that this is not the estimator used by Gunsilius and Schennach (2023,
#' Section 2.5), who compute the convex potential of the Brenier map
#' numerically, either on a grid with finite differences and a kernel density
#' estimate or by a series approximation, with a penalty against nonconvex
#' solutions. The difference step must be large relative to the error of
#' the discrete map: with `step = 1` and about one observation per grid point,
#' the differences are dominated by the random spacings of the data, and the
#' logarithm then biases the eigenvalues downwards. By default, the differences
#' are therefore taken over two grid steps (`step = 2`), the grid has about two
#' observations per grid point, and the differences are one-sided at the edge
#' of the grid (`boundary = "clamp"`). These choices were the most accurate in a
#' simulation with known \eqn{\bar J} (Gaussian and nonlinear Brenier maps,
#' \eqn{d = 2, 3}, \eqn{n} from 1000 to 4000), with eigenvalue errors of about
#' 0.02 to 0.05 and leading directions within about 1 to 2 degrees of the true
#' ones.
#'
#' @param X Numeric data matrix, observations in rows.
#' @param n_components Number of directions to keep (default `ncol(X)`).
#' @param n_grid Number of normal quantiles per axis, \eqn{n_g}. The grid has
#'   `n_grid^ncol(X)` points, at most `nrow(X)`. Defaults to
#'   `max(3, floor((nrow(X) / 2)^(1 / ncol(X))))`, that is, about two
#'   observations per grid point.
#' @param step Finite difference step, in grid indices (default 2).
#' @param boundary Treatment of the grid points within `step` of the edge of the
#'   grid: `"clamp"` (default) cuts the difference windows back there, which
#'   gives one-sided differences; `"nearest"` gives them the Jacobian of the
#'   nearest grid point with a full window; `"trim"` leaves them out, which
#'   truncates the normal tails and underestimates the eigenvalues of directions
#'   along which the map is more nonlinear in the tails.
#' @param eps_curves Tikhonov smoothing \eqn{(I + \varepsilon L'L)^{-1}} of the
#'   returned principal curves, as in [nlospca()], for display purposes (`0`, the
#'   default, for none). It does not affect the estimation.
#' @param tol Jacobian estimates with an eigenvalue `<= tol` are left out of the
#'   average.
#' @param ot_method Optimal transport algorithm of [transport::transport()],
#'   `"networkflow"` (default) or `"shortsimplex"`; both are exact.
#' @param nthreads Number of threads for the potentials computed by [predict()]
#'   for new observations; see [nlospca()].
#'
#' @return An object of class `gspca`, a list with the following components,
#'   where \eqn{n} denotes the number of observations, \eqn{p} the number of
#'   variables and \eqn{k} the number of retained directions:
#'   \item{directions}{\eqn{p \times k}{p x k} matrix of the directions
#'     \eqn{u_1, \ldots, u_k}, in its columns.}
#'   \item{eigenvalues}{All \eqn{p} eigenvalues of \eqn{\bar J}, in decreasing
#'     order.}
#'   \item{Jbar}{The estimated \eqn{p \times p}{p x p} matrix \eqn{\bar J}.}
#'   \item{scores}{\eqn{n \times k}{n x k} matrix of the scores
#'     \eqn{\langle T(X_i), u_j \rangle}.}
#'   \item{curves}{List of the \eqn{k} principal curves \eqn{T^{-1}(r u_j)},
#'     each a matrix of points in the space of the data.}
#'   \item{TX}{\eqn{n \times p}{n x p} matrix of the transformed observations
#'     \eqn{T(X_i)}.}
#'   \item{grid}{List describing the reference grid (`y`, the quantiles on each
#'     axis; `G`, the grid points; `n_grid`) and the discrete map \eqn{T^{-1}}
#'     at the grid points (`Yg`).}
#'   \item{X}{The data matrix.}
#'   \item{n, p}{Number of observations and of variables.}
#'   \item{settings}{List of the settings used (`n_grid`, `step`, `boundary`,
#'     `eps_curves`, `tol`, `ot_method`).}
#'   \item{diagnostics}{List with the proportion `dropped` of grid points left
#'     out of \eqn{\bar J}.}
#'   \item{nthreads}{Number of threads used by [predict()].}
#'   \item{call}{The matched call.}
#'
#' @references Gunsilius, F. and Schennach, S. (2023). Independent nonlinear
#'   component analysis. *Journal of the American Statistical Association*,
#'   118(542), 1305--1318. \doi{10.1080/01621459.2021.1990768}
#' @seealso [nlospca()] for the directional nonlinear PCA based on the spherical
#'   uniform distribution, and [ospca()] for the linear directional PCA.
#' @examples
#' ## Gaussian data: the directions are the linear principal components and the
#' ## eigenvalues are 0.5 * log(eigen(Sigma)$values)
#' set.seed(1)
#' Sigma <- matrix(c(4, 1, 1, 1), 2)
#' Xg <- matrix(rnorm(4000), ncol = 2) %*% chol(Sigma)
#' fit <- gspca(Xg)
#' fit$directions; eigen(Sigma)$vectors
#' fit$eigenvalues; 0.5 * log(eigen(Sigma)$values)
#'
#' ## a curved cloud
#' x <- rnorm(2000)
#' X <- cbind(x, x^2 + rnorm(2000, sd = 0.5))
#' fit <- gspca(X)
#' plot(fit)
#' screeplot(fit)
#' @export
gspca <- function(X, n_components = ncol(X), n_grid = NULL, step = 2L,
                  boundary = c("clamp", "nearest", "trim"), eps_curves = 0, tol = 1e-8,
                  ot_method = c("networkflow", "shortsimplex"), nthreads = NULL) {
  ot_method <- match.arg(ot_method)
  boundary <- match.arg(boundary)
  X <- .check_matrix(X)
  n <- nrow(X); p <- ncol(X)
  n_components <- as.integer(n_components)
  if (length(n_components) != 1L || n_components < 1L || n_components > p)
    stop("n_components must be an integer in [1, ncol(X)].", call. = FALSE)
  step <- as.integer(step)
  if (length(step) != 1L || is.na(step) || step < 1L)
    stop("step must be a positive integer.", call. = FALSE)
  if (length(eps_curves) != 1L || !is.finite(eps_curves) || eps_curves < 0)
    stop("eps_curves must be a single non-negative number.", call. = FALSE)
  if (length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("tol must be a single non-negative number.", call. = FALSE)
  if (is.null(n_grid)) n_grid <- max(3L, floor((n / 2)^(1 / p)))
  n_grid <- as.integer(n_grid)
  min_ng <- if (boundary == "clamp") 3L else 2L * step + 1L
  if (length(n_grid) != 1L || is.na(n_grid) || n_grid < min_ng)
    stop(sprintf("n_grid must be a single integer >= %d for step = %d and boundary = \"%s\"; more data are needed.",
                 min_ng, step, boundary), call. = FALSE)
  if (as.numeric(n_grid)^p > n)
    stop("the grid has n_grid^ncol(X) > nrow(X) points: decrease n_grid.", call. = FALSE)
  nthreads <- if (is.null(nthreads)) .default_threads() else max(1L, as.integer(nthreads))

  ## 1.-2. reference grid and discrete Brenier map
  grid <- .gs_grid(n_grid, p)
  cp <- .gs_coupling(grid$G, X, ot_method)

  ## 3. Jbar = average matrix logarithm of the symmetrised Jacobians of T^{-1}
  jb <- .gs_jbar(cp$Yg, grid, tol, step, boundary)
  eg <- eigen(jb$Jbar, symmetric = TRUE)
  idx <- seq_len(n_components)
  U <- eg$vectors[, idx, drop = FALSE]
  labels <- paste0("GSPC", idx)
  colnames(U) <- labels
  rownames(U) <- colnames(X)

  Z <- cp$TX %*% U
  dimnames(Z) <- list(rownames(X), labels)
  fit <- structure(
    list(directions = U, eigenvalues = eg$values, Jbar = jb$Jbar, scores = Z,
         curves = NULL, TX = cp$TX,
         grid = list(y = grid$y, G = grid$G, Yg = cp$Yg, n_grid = n_grid),
         X = X, n = n, p = p,
         settings = list(n_grid = n_grid, step = step, boundary = boundary,
                         eps_curves = eps_curves, tol = tol, ot_method = ot_method),
         diagnostics = list(dropped = jb$dropped), nthreads = nthreads,
         call = match.call()),
    class = "gspca")
  fit$curves <- lapply(idx, function(j) .gs_curve(fit, U[, j], eps_curves = eps_curves))
  names(fit$curves) <- labels

  if (jb$dropped > 0.05)
    warning(sprintf("%.1f%% of the grid points had a Jacobian estimate that is not positive definite and were left out of Jbar; consider a coarser grid or a larger step.",
                    100 * jb$dropped), call. = FALSE)
  fit
}

# Tensor grid of standard normal quantiles (expand.grid order: first axis
# fastest), with the per-axis quantiles y and index strides.
.gs_grid <- function(n_grid, p) {
  y1 <- stats::qnorm((seq_len(n_grid) - 0.5) / n_grid)
  y <- rep(list(y1), p)
  G <- as.matrix(do.call(expand.grid, y))
  dimnames(G) <- NULL
  list(y = y, G = G, n_grid = rep(n_grid, p),
       stride = as.integer(cumprod(c(1, rep(n_grid, p)))[seq_len(p)]))
}

# Exact optimal coupling between the grid (mass 1/m each) and the sample (mass
# 1/n each). Returns the barycentric images Yg of the grid points (discrete
# T^{-1}) and the barycentric images TX of the observations (discrete T).
.gs_coupling <- function(G, X, method) {
  m <- nrow(G); n <- nrow(X); p <- ncol(X)
  ## `transport` needs points of dimension >= 2: for p = 1 add a zero coordinate
  ## (costs unchanged)
  G2 <- if (p == 1L) cbind(G, 0) else G
  X2 <- if (p == 1L) cbind(X, 0) else X
  a <- transport::transport(transport::wpp(G2, rep(n, m)), transport::wpp(X2, rep(m, n)),
                            p = 2, method = method)
  fg <- factor(a$from, levels = seq_len(m)); fx <- factor(a$to, levels = seq_len(n))
  Yg <- matrix(0, m, p); TX <- matrix(0, n, p)
  for (h in seq_len(p)) {
    Yg[, h] <- tapply(a$mass * X[a$to, h], fg, sum) / n
    TX[, h] <- tapply(a$mass * G[a$from, h], fx, sum) / m
  }
  list(Yg = Yg, TX = TX)
}

# Average matrix logarithm of the symmetrised Jacobians of T^{-1} at the grid
# points, by central differences over +/- `step` grid steps along each axis.
# Grid points within `step` of the edge of the grid:
#   "nearest": take the Jacobian of the nearest grid point with a full window;
#   "clamp":   windows cut back at the edges (one-sided differences there);
#   "trim":    left out (this truncates the normal tails).
# Grid points whose symmetrised Jacobian has an eigenvalue <= tol (noise) are
# left out of the average; their proportion is returned as `dropped`.
.gs_jbar <- function(Yg, grid, tol, step = 1L, boundary = "nearest") {
  p <- ncol(Yg); ng <- grid$n_grid[1L]
  sub <- as.matrix(do.call(expand.grid, rep(list(seq_len(ng)), p)))   # multi-indices
  if (boundary == "trim")
    sub <- sub[apply(sub > step & sub <= ng - step, 1L, all), , drop = FALSE]
  at <- if (boundary == "nearest") pmin(pmax(sub, 1L + step), ng - step) else sub
  node <- 1L + drop((at - 1L) %*% grid$stride)
  D <- vector("list", p)               # D[[a]][r, ] = d T^{-1} / d x_a at point r
  for (a in seq_len(p)) {
    lo <- pmax(at[, a] - step, 1L); hi <- pmin(at[, a] + step, ng)
    D[[a]] <- (Yg[node + (hi - at[, a]) * grid$stride[a], , drop = FALSE] -
                 Yg[node + (lo - at[, a]) * grid$stride[a], , drop = FALSE]) /
      (grid$y[[a]][hi] - grid$y[[a]][lo])
  }
  N <- nrow(sub); Jsum <- matrix(0, p, p); used <- 0L
  for (r in seq_len(N)) {
    J <- matrix(vapply(D, function(Da) Da[r, ], numeric(p)), p, p)   # J[h, a]
    e <- eigen((J + t(J)) / 2, symmetric = TRUE)
    if (any(e$values <= tol)) next
    Jsum <- Jsum + e$vectors %*% (log(e$values) * t(e$vectors))
    used <- used + 1L
  }
  if (used == 0L)
    stop("no grid point with a positive-definite Jacobian estimate; use a coarser grid or a larger step.",
         call. = FALSE)
  list(Jbar = Jsum / used, dropped = 1 - used / N)
}

# T^{-1} at the rows of Z by multilinear interpolation on the grid (clamped to
# the grid box).
.gs_tinv <- function(object, Z) {
  gr <- object$grid; p <- object$p; ng <- gr$n_grid
  Z <- matrix(Z, ncol = p)
  y <- gr$y[[1L]]
  stride <- as.integer(cumprod(c(1, rep(ng, p)))[seq_len(p)])
  lo <- w <- matrix(0, nrow(Z), p)
  for (a in seq_len(p)) {
    z <- pmin(pmax(Z[, a], y[1L]), y[ng])
    i <- findInterval(z, y, rightmost.closed = TRUE, all.inside = TRUE)
    lo[, a] <- i
    w[, a] <- (z - y[i]) / (y[i + 1L] - y[i])
  }
  out <- matrix(0, nrow(Z), p)
  for (b in 0:(2^p - 1L)) {
    bit <- (b %/% 2^(seq_len(p) - 1L)) %% 2L
    wt <- rep(1, nrow(Z)); pos <- rep(1L, nrow(Z))
    for (a in seq_len(p)) {
      wt <- wt * (if (bit[a]) w[, a] else 1 - w[, a])
      pos <- pos + (lo[, a] - 1L + bit[a]) * stride[a]
    }
    out <- out + wt * gr$Yg[pos, , drop = FALSE]
  }
  colnames(out) <- colnames(object$X)
  out
}

# Principal curve T^{-1}(r u), r from -r_max to r_max, with r_max the largest
# radius keeping r u inside the grid box.
.gs_curve <- function(object, u, n_points = 101L, eps_curves = 0) {
  rmax <- max(abs(object$grid$y[[1L]])) / max(abs(u))
  r <- seq(-rmax, rmax, length.out = n_points)
  C <- .gs_tinv(object, outer(r, u))
  if (eps_curves > 0) C <- .tikhonov(n_points, eps_curves) %*% C   # optional smoothing (display)
  C
}

# T at new observations: grid point of the piecewise-constant cyclically
# monotone extension. The potentials (Karp, O(m^3) for m grid points) are only
# needed here and are therefore computed on demand.
.gs_t <- function(object, X) {
  gr <- object$grid
  psi <- .ot_potentials(gr$G, gr$Yg, object$nthreads)$psi
  gr$G[.fbar_index(X, gr$G, gr$Yg, psi), , drop = FALSE]
}

## ---------------------------------------------------------------- methods --

#' Methods for `gspca` objects
#'
#' Printing, fitted values and residuals of a fitted [gspca()] model. The fitted
#' values are the reconstructions with the first `n_components` directions (see
#' [reconstruct()]), and the residuals are the data minus their
#' reconstructions.
#'
#' @param x,object A `gspca` object.
#' @param n_components Number of directions to use where applicable.
#' @param ... Ignored.
#' @return The `print` method returns its argument invisibly; `fitted` and
#'   `residuals` return numeric matrices.
#' @examples
#' set.seed(1)
#' x <- rnorm(1000)
#' X <- cbind(x, x^2 + rnorm(1000, sd = 0.5))
#' fit <- gspca(X)
#' fit
#' head(fitted(fit, n_components = 1))
#' head(residuals(fit, n_components = 1))
#' @name gspca-methods
#' @rdname gspca-methods
#' @export
print.gspca <- function(x, ...) {
  cat("Independent nonlinear component analysis (GS-PCA) fit\n")
  cat("  Observations:        ", x$n, "\n", sep = "")
  cat("  Variables (d):       ", x$p, "\n", sep = "")
  cat("  Reference grid:      ", x$grid$n_grid, "^", x$p, " = ", nrow(x$grid$G),
      " normal quantiles\n", sep = "")
  if (!is.null(x$settings))
    cat("  Differences:         step ", x$settings$step, ", boundary \"",
        x$settings$boundary, "\"; ", sprintf("%.1f%%", 100 * x$diagnostics$dropped),
        " of the grid points left out\n", sep = "")
  cat("  Directions (k):      ", ncol(x$directions), "\n", sep = "")
  cat("  Eigenvalues of Jbar (decreasing):\n    ",
      paste(sprintf("%.4g", x$eigenvalues), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Summarize a GS-PCA fit
#'
#' @param object A `gspca` object.
#' @param ... Ignored.
#' @return An object of class `summary.gspca`.
#' @examples
#' set.seed(1)
#' x <- rnorm(1000)
#' X <- cbind(x, x^2 + rnorm(1000, sd = 0.5))
#' summary(gspca(X))
#' @export
summary.gspca <- function(object, ...) {
  out <- list(n = object$n, p = object$p, n_components = ncol(object$directions),
              eigenvalues = object$eigenvalues, directions = object$directions,
              n_grid = object$grid$n_grid, diagnostics = object$diagnostics)
  class(out) <- "summary.gspca"
  out
}

#' @rdname summary.gspca
#' @param x A `summary.gspca` object.
#' @export
print.summary.gspca <- function(x, ...) {
  cat("Independent nonlinear component analysis (GS-PCA) summary\n")
  cat("Dimensions: ", x$n, " x ", x$p, "   directions: ", x$n_components,
      "   grid: ", x$n_grid, "^", x$p, "\n", sep = "")
  cat("Grid points left out (Jacobian not positive definite): ",
      sprintf("%.1f%%", 100 * x$diagnostics$dropped), "\n\n", sep = "")
  tab <- rbind(eigenvalue = x$eigenvalues[seq_len(x$n_components)], x$directions)
  cat("Eigenvalues of Jbar and directions:\n"); print(round(tab, 4))
  invisible(x)
}

#' Predict GS-PCA scores
#'
#' Scores \eqn{\langle T(x), u_j \rangle} of new observations, where \eqn{T(x)}
#' is the grid point given by the piecewise constant cyclically monotone
#' extension of the discrete Brenier map. Note that the scores of the
#' observations used in the fit are based on the optimal coupling itself.
#'
#' @param object A `gspca` object.
#' @param newdata Optional numeric matrix with `object$p` columns.
#' @param n_components Number of directions to use. Defaults to all.
#' @param ... Ignored.
#' @return Matrix of scores.
#' @examples
#' set.seed(1)
#' x <- rnorm(1000)
#' X <- cbind(x, x^2 + rnorm(1000, sd = 0.5))
#' fit <- gspca(X)
#' predict(fit, newdata = cbind(c(-1, 0, 1), c(1, 0, 1)))
#' @export
predict.gspca <- function(object, newdata = NULL, n_components = NULL, ...) {
  idx <- .component_index(object, n_components)
  if (is.null(newdata)) return(object$scores[, idx, drop = FALSE])
  X <- .check_matrix(newdata)
  if (ncol(X) != object$p) stop("newdata must have the same number of columns.", call. = FALSE)
  Z <- .gs_t(object, X) %*% object$directions[, idx, drop = FALSE]
  dimnames(Z) <- list(rownames(X), colnames(object$directions)[idx])
  Z
}

#' @rdname scores
#' @export
scores.gspca <- function(object, n_components = NULL, ...)
  object$scores[, .component_index(object, n_components), drop = FALSE]

#' @rdname directions
#' @export
directions.gspca <- function(object, n_components = NULL, ...)
  object$directions[, .component_index(object, n_components), drop = FALSE]

#' @rdname gspca-methods
#' @export
fitted.gspca <- function(object, n_components = NULL, ...)
  reconstruct(object, n_components = n_components)

#' @rdname gspca-methods
#' @export
residuals.gspca <- function(object, n_components = NULL, ...)
  object$X - reconstruct(object, n_components = n_components)

#' @rdname reconstruct
#' @details For a `gspca` fit, the reconstruction is the representation
#'   \eqn{T^{-1}(\sum_{j \le k} \langle T(X), u_j \rangle u_j)} of Gunsilius and
#'   Schennach (2023) on the manifold
#'   \eqn{T^{-1}(\mathrm{span}(u_1, \ldots, u_k))}.
#' @export
reconstruct.gspca <- function(object, newdata = NULL, n_components = NULL,
                              add_center = TRUE, ...) {
  idx <- .component_index(object, n_components)
  U <- object$directions[, idx, drop = FALSE]
  Z <- predict(object, newdata, n_components = length(idx))
  Xhat <- .gs_tinv(object, Z %*% t(U))
  rownames(Xhat) <- rownames(if (is.null(newdata)) object$X else newdata)
  Xhat
}

#' @rdname screeplot.ospca
#' @export
screeplot.gspca <- function(x, type = c("both", "objective", "mse"), ...) {
  type <- match.arg(type)
  k <- ncol(x$directions)
  ev <- x$eigenvalues[seq_len(k)]
  labels <- colnames(x$directions)
  mse <- if (type == "objective") NULL else
    vapply(seq_len(k), function(j)
      reconstruction_mse(x$X, reconstruct(x, n_components = j)), numeric(1))
  ylab <- expression(paste("eigenvalue of ", bar(J)))
  if (type == "objective") {
    graphics::barplot(ev, names.arg = labels, col = "#3b5b92", border = NA,
                      ylab = ylab, main = "Eigenvalues", ...)
  } else if (type == "mse") {
    plot(seq_len(k), mse, type = "b", pch = 16, col = "#c0392b",
         xlab = "number of directions (k)", ylab = "reconstruction MSE",
         main = "Reconstruction error", ...)
  } else {
    op <- graphics::par(mar = c(5, 4, 4, 4) + 0.3); on.exit(graphics::par(op))
    bp <- graphics::barplot(ev, names.arg = labels, col = "#3b5b92", border = NA,
                            ylab = ylab,
                            main = "Eigenvalues and reconstruction error", ...)
    graphics::par(new = TRUE)
    plot(bp, mse, type = "b", pch = 16, col = "#c0392b", axes = FALSE,
         xlab = "", ylab = "", xlim = graphics::par("usr")[1:2], xaxs = "i",
         ylim = c(0, max(mse)))
    graphics::axis(4, col = "#c0392b", col.axis = "#c0392b")
    graphics::mtext("reconstruction MSE", side = 4, line = 2.5, col = "#c0392b")
  }
  invisible(x)
}

#' Plot a GS-PCA fit
#'
#' For bivariate data, draws the data together with the principal curves
#' \eqn{T^{-1}(r u_j)}. For trivariate data, draws the three pairwise
#' scatterplots with the projected curves, and for higher dimensions, plots two
#' chosen columns of the scores against each other.
#'
#' @param x A `gspca` object.
#' @param y Ignored.
#' @param choices For `d > 3`, the two score columns to plot.
#' @param col_curves Colors of the principal curves (recycled).
#' @param ... Passed to the underlying plotting call.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' x <- rnorm(1000)
#' X <- cbind(x, x^2 + rnorm(1000, sd = 0.5))
#' plot(gspca(X, eps_curves = 1))
#' @export
plot.gspca <- function(x, y = NULL, choices = 1:2,
                       col_curves = c("navy", "#c0392b", "#27ae60", "#8e44ad"), ...) {
  col_curves <- rep_len(col_curves, length(x$curves))
  nm <- colnames(x$X); if (is.null(nm)) nm <- paste0("X", seq_len(x$p))
  draw <- function(i, j, ...) {
    plot(x$X[, i], x$X[, j], xlab = nm[i], ylab = nm[j], pch = 16,
         col = "grey60", asp = 1, ...)
    for (h in seq_along(x$curves)) {
      C <- x$curves[[h]]
      graphics::lines(C[, i], C[, j], col = col_curves[h], lwd = 2)
      graphics::text(C[nrow(C), i], C[nrow(C), j], labels = h, col = col_curves[h], pos = 3)
    }
  }
  if (x$p == 2L) {
    draw(1L, 2L, ...)
  } else if (x$p == 3L) {
    op <- graphics::par(mfrow = c(1L, 3L)); on.exit(graphics::par(op))
    draw(1L, 2L, ...); draw(1L, 3L, ...); draw(2L, 3L, ...)
  } else {
    choices <- as.integer(choices)
    if (length(choices) != 2L) stop("choices must contain exactly two indices.", call. = FALSE)
    Z <- x$scores
    plot(Z[, choices[1]], Z[, choices[2]], xlab = colnames(Z)[choices[1]],
         ylab = colnames(Z)[choices[2]], pch = 16, col = "grey50", ...)
  }
  invisible(x)
}
