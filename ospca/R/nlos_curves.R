# Discrete (Tikhonov-regularised) sign curves and the nonlinear one-sided
# criterion M(u).

# (I + eps L'L)^{-1} with L the (m - 1) x m first-difference matrix.
.tikhonov <- function(m, eps) {
  if (eps <= 0 || m < 2L) return(diag(m))
  L <- diff(diag(m))
  solve(diag(m) + eps * crossprod(L))
}

# Regularised discrete curves Qbar(rho_j u), rho_j = j / (mR + 1), j = 1..mR,
# one mR x d matrix per row u of U.
.nlos_curve_points <- function(tf, U, mR, eps) {
  U <- matrix(U, ncol = ncol(tf$Yi))
  rho <- seq_len(mR) / (mR + 1)
  Z <- do.call(rbind, lapply(seq_len(nrow(U)), function(k) outer(rho, U[k, ])))
  QZ <- .qbar(Z, tf$Yi, tf$psi, tf$kappa, tf$nthreads)
  S <- .tikhonov(mR, eps)
  lapply(seq_len(nrow(U)), function(k)
    S %*% QZ[(k - 1L) * mR + seq_len(mR), , drop = FALSE])
}

# Trimming level m_delta = floor(mR (1 - delta)), kept within [2, mR].
.trim_index <- function(mR, delta) {
  as.integer(min(max(floor(mR * (1 - delta)), 2L), mR))
}

# Criterion (1 / m_delta) sum_{j = 2}^{m_delta} L_n^2(Qbar(r_j u)), where L_n is
# the polygonal length from r_1 to r_j of the (regularised) curve C.
.nlos_crit <- function(C, md) {
  ell <- c(0, cumsum(sqrt(rowSums(diff(C)^2))))
  sum(ell[seq_len(md)]^2) / md
}

.check_nlospca <- function(object) {
  if (!inherits(object, "nlospca")) stop("object must be an nlospca fit.", call. = FALSE)
  invisible(object)
}

#' Nonlinear one-sided PCA objective
#'
#' The empirical criterion maximized by [nlospca()], that is, the trimmed mean
#' of the squared polygonal lengths of the (Tikhonov regularized) sign curve
#' \eqn{r \mapsto \overline{Q}_\pm^{(n)}(r u)},
#' \deqn{M(u) = \frac{1}{m_\delta} \sum_{j=2}^{m_\delta}
#'   L_n^2\big(\overline{Q}_\pm^{(n)}(r_j u)\big), \qquad r_j = j / (m_R + 1),}
#' with \eqn{m_\delta = \lfloor m_R (1 - \delta) \rfloor}.
#'
#' @param u Direction (a numeric vector of length `object$p`) or a matrix with
#'   one direction per row. Directions are normalized internally.
#' @param object An `nlospca` fit, which provides the empirical quantile
#'   function.
#' @param mR Number of radii at which the curves are evaluated; the default is
#'   the value used in the fit.
#' @param eps Tikhonov regularization parameter; the default is the value used
#'   in the fit.
#' @param delta Trimming proportion; the default is the value used in the fit.
#' @return Numeric vector of objective values, one per direction.
#' @seealso [nlospca()], [nlos_curve()]
#' @examples
#' set.seed(1)
#' x <- rnorm(300); X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' fit <- nlospca(X)
#' nlos_objective(rbind(c(0, 1), c(0, -1)), fit)
#' @export
nlos_objective <- function(u, object, mR = object$settings$mR,
                          eps = object$settings$eps, delta = object$settings$delta) {
  .check_nlospca(object)
  U <- .check_directions(u, object$p)
  md <- .trim_index(mR, delta)
  vapply(.nlos_curve_points(object$transport, U, mR, eps), .nlos_crit, numeric(1), md = md)
}

#' Nonlinear one-sided sign curve
#'
#' Evaluates the (Tikhonov regularized) empirical sign curve
#' \eqn{\overline{Q}_\pm^{(n)}(r_j u)}, \eqn{r_j = j / (m_R + 1)}, along a
#' direction `u`, starting at the center \eqn{Q_\pm^{(n)}(0)}. For a principal
#' direction, this is the corresponding empirical principal curve.
#'
#' @inheritParams nlos_objective
#' @param u Direction, a numeric vector of length `object$p`.
#' @return An `(mR + 1)` by `p` matrix of curve points, starting at the center.
#' @seealso [nlospca()], [nlos_objective()]
#' @examples
#' set.seed(1)
#' x <- rnorm(300); X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' fit <- nlospca(X)
#' head(nlos_curve(directions(fit)[, 1], fit))
#' @export
nlos_curve <- function(u, object, mR = object$settings$mR, eps = object$settings$eps) {
  .check_nlospca(object)
  u <- .check_direction(u, object$p)
  C <- rbind(object$center, .nlos_curve_points(object$transport, u, mR, eps)[[1L]])
  colnames(C) <- colnames(object$X)
  C
}

#' Empirical center-outward quantile and distribution functions
#'
#' `nlos_quantile()` evaluates the interpolated empirical center-outward
#' quantile function \eqn{\overline{Q}_\pm^{(n)}} (piecewise constant for
#' `kappa = 0` and its Yosida regularization otherwise) at points of the open
#' unit ball. `nlos_distribution()` evaluates the piecewise constant cyclically
#' monotone extension of the empirical center-outward distribution function
#' \eqn{F_\pm^{(n)}}, which maps each point to a point of the target grid.
#'
#' @param object An `nlospca` fit.
#' @param z Matrix (or vector) of points in the unit ball, with `object$p`
#'   columns.
#' @param x Matrix (or vector) of points in the data space, with `object$p`
#'   columns.
#' @return A matrix with one row per input point.
#' @seealso [nlospca()]
#' @examples
#' set.seed(1)
#' x <- rnorm(300); X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' fit <- nlospca(X)
#' nlos_quantile(fit, rbind(c(0, 0.5), c(0.5, 0)))
#' nlos_distribution(fit, X[1:3, ])
#' @export
nlos_quantile <- function(object, z) {
  .check_nlospca(object)
  z <- .check_matrix(matrix(z, ncol = object$p))
  tf <- object$transport
  Q <- .qbar(z, tf$Yi, tf$psi, tf$kappa, tf$nthreads)
  colnames(Q) <- colnames(object$X)
  Q
}

#' @rdname nlos_quantile
#' @export
nlos_distribution <- function(object, x) {
  .check_nlospca(object)
  x <- .check_matrix(matrix(x, ncol = object$p))
  tf <- object$transport
  tf$Gi[.fbar_index(x, tf$Gi, tf$Yi, tf$psi), , drop = FALSE]
}
