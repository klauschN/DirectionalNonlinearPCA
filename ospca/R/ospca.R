#' One-sided (directional) principal component analysis
#'
#' In classical principal component analysis a direction and its negative are
#' equivalent, as the components are based on variances. One-sided
#' (directional) principal component analysis considers the two half-lines of a
#' direction separately.
#'
#' Let \eqn{X} be a centered random vector and \eqn{u} a unit vector. The
#' positive part of the projection of \eqn{X} onto \eqn{u} is
#' \deqn{\langle u, X \rangle^+ = \max(\langle u, X \rangle, 0),}{<u, X>^+ = max(<u, X>, 0),}
#' that is, the projection if \eqn{X} lies on the side of \eqn{u} and zero
#' otherwise. The importance of the direction \eqn{u} is measured by the
#' one-sided objective
#' \deqn{J(u) = E[(\langle u, X \rangle^+)^2],}{J(u) = E[(<u, X>^+)^2],}
#' the second moment of the positive part of the projection, which is estimated
#' by the average of \eqn{(\langle u, x_i \rangle^+)^2}{(<u, x_i>^+)^2} over
#' the centered observations \eqn{x_1, \ldots, x_n} (see [os_objective()]).
#' Note that
#' \deqn{\langle u, X \rangle^2 = (\langle u, X \rangle^+)^2 + (\langle -u, X \rangle^+)^2,}{<u, X>^2 = (<u, X>^+)^2 + (<-u, X>^+)^2,}
#' so that the variance criterion of classical PCA is \eqn{J(u) + J(-u)}. For a
#' symmetric distribution, \eqn{J(u) = J(-u)}, but for an asymmetric one the
#' two half-lines of a direction can have quite different objectives. For a
#' retained direction \eqn{u_j}, the one-sided scores (the directional
#' principal components) are the positive parts
#' \eqn{\langle u_j, x_i \rangle^+}{<u_j, x_i>^+}.
#'
#' The function first computes the greedy one-sided orthonormal axes
#' \eqn{a_1, \ldots, a_d} of the centered data (see [os_axes()]), then evaluates
#' the objective for the \eqn{2d} candidate directions
#' \eqn{\pm a_1, \ldots, \pm a_d}{+/- a_1, ..., +/- a_d} and keeps the `n_components` directions with
#' the largest objective, in decreasing order. By default, each axis is found
#' by a multi-start generalized power iteration; `method = "bfgs"` uses a
#' quasi-Newton optimizer instead.
#'
#' A warning is given if (i) an axis is a near-tie, that is, its objective is
#' within the relative tolerance `tol_tie` of the best clearly different
#' candidate, or (ii) two consecutive retained objectives differ by less than
#' `tol_tie` (relative). In both cases the solution may not be unique and should
#' be used with some caution.
#'
#' @param X Numeric data matrix, observations in rows.
#' @param n_components Number of directions to keep. Defaults to `2 * ncol(X)`,
#'   which is also the maximum.
#' @param center `TRUE`, `FALSE` or a numeric vector giving the center.
#' @param n_starts Number of random starting directions per axis, passed to
#'   [os_axes()].
#' @param n_polish Number of best starting directions that are refined, passed
#'   to [os_axes()].
#' @param tol_tie Relative tolerance below which two objectives are considered
#'   tied (see Details), passed to [os_axes()].
#' @param method Optimizer for the axes, `"power"` (default) or `"bfgs"`; see
#'   [os_axes()].
#' @param maxit Maximum number of power iteration steps per start, passed to
#'   [os_axes()].
#' @param tol Convergence tolerance of the power iteration, passed to
#'   [os_axes()].
#'
#' @return An object of class `ospca`, a list with the following components,
#'   where \eqn{n} denotes the number of observations, \eqn{p} the number of
#'   variables and \eqn{k} the number of retained directions:
#'   \item{directions}{\eqn{p \times k}{p x k} matrix of the one-sided
#'     directions, in its columns and in decreasing order of the objective.}
#'   \item{scores}{\eqn{n \times k}{n x k} matrix of the one-sided scores
#'     \eqn{\langle u_j, x_i \rangle^+}{<u_j, x_i>^+}, that is, the positive
#'     parts of the projections of the centered data onto the directions.}
#'   \item{X_centered}{The centered data matrix.}
#'   \item{center}{The center that was subtracted from the data.}
#'   \item{objective}{The one-sided objectives of the \eqn{k} directions.}
#'   \item{diagnostics}{Data frame with the axis diagnostics of [os_axes()]
#'     (`objective`, `gap` and `tie` for each axis).}
#'   \item{n, p}{Number of observations and of variables.}
#'   \item{call}{The matched call.}
#'
#' @references Hallin, M., Nordhausen, K. and Radojicic, U. Directional
#'   nonlinear principal components: a measure transportation approach.
#'
#' @seealso [os_axes()], [screeplot()], [reconstruct()], [reconstruction_mse()]
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(1000), ncol = 2) %*% diag(c(3, 1)),
#'            cbind(abs(rnorm(300, 5)), rnorm(300)))
#' fit <- ospca(X)
#' fit
#' screeplot(fit)
#' @export
ospca <- function(X, n_components = 2L * ncol(X), center = TRUE,
                   n_starts = 200L, n_polish = 10L, tol_tie = 1e-3,
                   method = c("power", "bfgs"), maxit = 200L, tol = 1e-10) {
  method <- match.arg(method)
  X <- .check_matrix(X)
  p <- ncol(X)
  n_components <- as.integer(n_components)
  if (length(n_components) != 1L || n_components < 1L)
    stop("n_components must be a positive integer.", call. = FALSE)
  n_components <- min(n_components, 2L * p)

  centered <- .center_matrix(X, center)
  Xc <- centered$X

  A    <- os_axes(Xc, n_starts = n_starts, n_polish = n_polish, tol_tie = tol_tie,
                   method = method, maxit = maxit, tol = tol)
  diag <- attr(A, "diagnostics")
  cand <- cbind(A, -A)
  vals <- vapply(seq_len(ncol(cand)), function(j) .os_obj(cand[, j], Xc), numeric(1))
  ord  <- order(vals, decreasing = TRUE)[seq_len(n_components)]

  U <- cand[, ord, drop = FALSE]
  objective <- vals[ord]
  Z <- pmax(Xc %*% U, 0)
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
  colnames(U) <- colnames(Z) <- paste0("OSPC", seq_len(ncol(U)))
  rownames(U) <- colnames(X)

  ## ---- tie / multiplicity warnings ----
  msgs <- character(0)
  if (any(diag$tie, na.rm = TRUE)) {
    kk <- diag$direction[which(diag$tie)]
    msgs <- c(msgs, sprintf(
      "axis/axes %s are near-tied (objective gap < tol_tie = %g): the one-sided maximiser may be non-unique (multiplicity); results can be unstable.",
      paste(kk, collapse = ", "), tol_tie))
  }
  if (length(objective) >= 2L) {
    rel <- abs(diff(objective)) / pmax(objective[-length(objective)], .Machine$double.eps)
    if (any(rel < tol_tie)) {
      jj <- which(rel < tol_tie)
      msgs <- c(msgs, sprintf(
        "retained directions %s have near-equal objective (relative gap < %g): their ordering is essentially arbitrary.",
        paste(sprintf("%d/%d", jj, jj + 1L), collapse = ", "), tol_tie))
    }
  }
  if (length(msgs)) warning(paste(msgs, collapse = "\n"), call. = FALSE)

  structure(
    list(directions = U, scores = Z, X_centered = Xc, center = centered$center,
         objective = objective, diagnostics = diag, n = nrow(X), p = p,
         call = match.call()),
    class = "ospca"
  )
}
