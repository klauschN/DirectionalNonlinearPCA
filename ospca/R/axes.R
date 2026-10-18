# Multi-start maximisation of the one-sided objective J over the unit sphere of the
# subspace spanned by the (orthonormal) columns of Q. The work is done in the
# m-dimensional subspace coordinates Y = X Q (an orthogonal change of basis that
# preserves inner products and angles), and the full direction is u = Q w.
#
# Two refiners are available for each random start:
#
#  * method = "power" (default): a structure-exploiting generalized power
#    iteration. Because J(w) = mean(pmax(Y w, 0)^2) is convex in w, the map
#    w <- normalize( t(Y) %*% pmax(Y w, 0) )   ( = normalize of grad J )
#    increases J at every step (Cauchy-Schwarz), converging to a stationary
#    point in a handful of cheap matrix-vector products -- no line search.
#  * method = "bfgs": the previous quasi-Newton polish via stats::optim, kept as
#    a fallback.
#
# Both are wrapped in a multi-start (best n_polish of n_starts random directions)
# to cope with the non-convexity (several local maxima) on the sphere.
#
# Multiplicity diagnostic: the relative gap between the attained objective and
# the best objective among sampled directions more than 45 degrees from the
# winner; a small gap flags a (near-)non-unique maximiser.
.optimize_direction <- function(X, Q, n_starts = 200L, n_polish = 10L,
                                tol_tie = 1e-3, method = c("power", "bfgs"),
                                maxit = 200L, tol = 1e-10) {
  method <- match.arg(method)
  m <- ncol(Q)
  if (m == 1L) {
    u <- as.numeric(Q[, 1])
    return(structure(u, objective = .os_obj(u, X), gap = NA_real_, tie = FALSE))
  }

  Y  <- X %*% Q                                  # n x m subspace coordinates
  Jm <- function(w) { s <- sqrt(sum(w^2)); mean(pmax(as.numeric(Y %*% (w / s)), 0)^2) }

  ## --- generalized power iteration (monotone ascent) ---
  run_power <- function(w) {
    w <- w / sqrt(sum(w^2))
    for (it in seq_len(maxit)) {
      g  <- as.numeric(crossprod(Y, pmax(as.numeric(Y %*% w), 0)))  # t(Y) relu(Y w)
      ng <- sqrt(sum(g^2))
      if (!is.finite(ng) || ng < 1e-14) break
      wn <- g / ng
      if (sqrt(sum((wn - w)^2)) < tol) { w <- wn; break }
      w <- wn
    }
    w
  }
  ## --- BFGS fallback ---
  run_bfgs <- function(w) {
    neg <- function(v) {
      nv <- sqrt(sum(v^2))
      if (!is.finite(nv) || nv < 1e-12) return(.Machine$double.xmax)
      -Jm(v)
    }
    opt <- stats::optim(w, neg, method = "BFGS", control = list(maxit = 500))
    opt$par / sqrt(sum(opt$par^2))
  }
  refine <- if (method == "power") run_power else run_bfgs

  starts     <- matrix(stats::rnorm(n_starts * m), n_starts, m)
  obj_starts <- apply(starts, 1L, Jm)
  ord        <- order(obj_starts, decreasing = TRUE)
  n_polish   <- max(1L, min(as.integer(n_polish), n_starts))

  cand <- vector("list", n_polish)
  cval <- numeric(n_polish)
  for (i in seq_len(n_polish)) {
    w <- refine(starts[ord[i], ])
    cand[[i]] <- w
    cval[i]   <- Jm(w)
  }
  best  <- which.max(cval)
  wbest <- cand[[best]]
  obest <- cval[best]
  ubest <- as.numeric(Q %*% wbest)
  ubest <- ubest / sqrt(sum(ubest^2))

  ## probe-based multiplicity gap (angles computed in subspace coords; Q is
  ## orthonormal so angles are preserved).
  dirs   <- starts / sqrt(rowSums(starts^2))
  far    <- abs(as.numeric(dirs %*% wbest)) < 0.707
  runner <- if (any(far)) max(obj_starts[far]) else NA_real_
  gap <- if (is.na(runner)) NA_real_
         else max(0, (obest - runner) / max(obest, .Machine$double.eps))
  tie <- isTRUE(gap < tol_tie)

  structure(ubest, objective = obest, gap = gap, tie = tie)
}

#' Greedy one-sided orthonormal axes
#'
#' Computes the greedy one-sided orthonormal axes \eqn{a_1, \ldots, a_d}{a_1, ..., a_d}: the
#' axis \eqn{a_1} maximizes the one-sided objective [os_objective()] over the
#' unit sphere, and each further axis \eqn{a_k} maximizes it over the unit
#' sphere of the orthogonal complement of the previous axes. Each maximization
#' uses a multi-start search in which the best `n_polish` of `n_starts` random
#' starting directions are refined (see Details).
#'
#' @param X Numeric matrix with observations in rows, assumed to be centered.
#' @param n_starts Number of random starting directions evaluated per axis.
#' @param n_polish Number of best starting directions that are refined; the
#'   best refined direction is kept.
#' @param tol_tie Relative objective gap below which an axis is declared a
#'   near-tie, that is, a possibly non-unique maximizer (see Value).
#' @param method Refinement of each start. `"power"` (default) is a generalized
#'   power iteration, `w <- normalize(t(Y) %*% pmax(Y w, 0))`, which increases
#'   the objective at every step and usually converges after a few iterations.
#'   `"bfgs"` uses `stats::optim(method = "BFGS")`.
#' @param maxit Maximum number of power iteration steps per start (ignored for
#'   `method = "bfgs"`).
#' @param tol Convergence tolerance for successive iterates of the power
#'   iteration.
#'
#' @return A `ncol(X)` by `ncol(X)` orthonormal matrix with the axes \eqn{a_i}
#'   in its columns and an attribute `"diagnostics"`, a data frame giving for
#'   each axis the attained `objective`, the relative `gap` and a logical `tie`
#'   flag. The `gap` is `(objective - best alternative) / objective`, where the
#'   best alternative is the largest objective among the sampled directions
#'   that are more than 45 degrees away from the chosen axis and its antipode. A
#'   small `gap` (below `tol_tie`) means that a clearly different direction does
#'   almost as well, i.e., the maximizer is effectively not unique.
#'
#' @details The one-sided objective
#'   \eqn{J(u) = E[(\langle u, X \rangle^+)^2]}{J(u) = E[(<u, X>^+)^2]}, with
#'   \eqn{\langle u, X \rangle^+ = \max(\langle u, X \rangle, 0)}{<u, X>^+ = max(<u, X>, 0)}
#'   (see [os_objective()]), is convex in \eqn{u}, but it is maximized over the sphere and can therefore
#'   have several local maxima. The power iteration uses the convexity to obtain
#'   a monotone ascent from each start, and the multi-start search addresses the
#'   local maxima. Refining the best of many starts makes the result much less
#'   dependent on the starting values than a single start.
#'
#' @seealso [ospca()], [os_objective()]
#' @examples
#' set.seed(1)
#' X <- scale(cbind(abs(rnorm(300, 4)), rnorm(300)), scale = FALSE)
#' A <- os_axes(X)
#' attr(A, "diagnostics")
#' A_bfgs <- os_axes(X, method = "bfgs")   # quasi-Newton optimizer
#' @export
os_axes <- function(X, n_starts = 200L, n_polish = 10L, tol_tie = 1e-3,
                     method = c("power", "bfgs"), maxit = 200L, tol = 1e-10) {
  method <- match.arg(method)
  X <- .check_matrix(X)
  p <- ncol(X)
  basis <- matrix(numeric(0), nrow = p, ncol = 0L)
  diag  <- data.frame(direction = integer(0), objective = numeric(0),
                      gap = numeric(0), tie = logical(0))
  for (k in seq_len(p)) {
    u <- .optimize_direction(X, .complement_basis(basis, p),
                             n_starts = n_starts, n_polish = n_polish,
                             tol_tie = tol_tie, method = method,
                             maxit = maxit, tol = tol)
    basis <- cbind(basis, as.numeric(u))
    diag <- rbind(diag, data.frame(direction = k, objective = attr(u, "objective"),
                                   gap = attr(u, "gap"), tie = attr(u, "tie")))
  }
  attr(basis, "diagnostics") <- diag
  basis
}
