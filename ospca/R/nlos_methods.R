#' Methods for `nlospca` objects
#'
#' Printing, fitted values and residuals of a fitted [nlospca()] model. The
#' fitted values are the reconstructions with the first `n_components`
#' directions (see [reconstruct()]), and the residuals are the data minus their
#' reconstructions.
#'
#' @param x,object An `nlospca` object.
#' @param n_components Number of directions to use where applicable.
#' @param ... Ignored.
#' @return The `print` method returns its argument invisibly; `fitted` and
#'   `residuals` return numeric matrices.
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' fit <- nlospca(X)
#' fit
#' head(fitted(fit, n_components = 2))
#' head(residuals(fit, n_components = 2))
#' @name nlospca-methods
#' @rdname nlospca-methods
#' @export
print.nlospca <- function(x, ...) {
  st <- x$settings
  cat("Nonlinear one-sided (directional) PCA fit\n")
  cat("  Observations:        ", x$n, "\n", sep = "")
  cat("  Variables (d):       ", x$p, "\n", sep = "")
  cat("  Target grid:         nR = ", st$nR, ", nS = ", st$nS, ", n0 = ", st$n0,
      "\n", sep = "")
  cat("  Directional dirs (k):", ncol(x$directions), " (max 2d = ", 2L * x$p, ")\n", sep = "")
  cat("  Criterion (decreasing):\n    ",
      paste(sprintf("%.4g", x$objective), collapse = ", "), "\n", sep = "")
  if (any(x$diagnostics$tie, na.rm = TRUE))
    cat("  Note: near-tied direction(s) detected (see $diagnostics).\n")
  invisible(x)
}

#' Summarize a nonlinear one-sided PCA fit
#'
#' @param object An `nlospca` object.
#' @param ... Ignored.
#' @return An object of class `summary.nlospca`.
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' summary(nlospca(X))
#' @export
summary.nlospca <- function(object, ...) {
  out <- list(n = object$n, p = object$p, n_components = ncol(object$directions),
              objective = object$objective, directions = object$directions,
              settings = object$settings, margin = object$transport$e0,
              diagnostics = object$diagnostics)
  class(out) <- "summary.nlospca"
  out
}

#' @rdname summary.nlospca
#' @param x A `summary.nlospca` object.
#' @export
print.summary.nlospca <- function(x, ...) {
  st <- x$settings
  cat("Nonlinear one-sided (directional) PCA summary\n")
  cat("Dimensions: ", x$n, " x ", x$p, "   directional directions: ",
      x$n_components, "\n", sep = "")
  cat("Target grid: nR = ", st$nR, ", nS = ", st$nS, ", n0 = ", st$n0,
      " (assignment: ", st$ot_method, ", potential margin ", signif(x$margin, 3),
      ")\n", sep = "")
  cat("Curves: mR = ", st$mR, ", mS = ", st$mS, ", eps = ", st$eps,
      ", eps_curves = ", st$eps_curves, ", paths = ", st$paths, ", delta = ",
      st$delta, ", kappa = ", signif(st$kappa, 3), "\n\n", sep = "")
  tab <- rbind(criterion = x$objective, x$directions)
  cat("Criterion values and directions:\n"); print(round(tab, 4))
  cat("\nDirection diagnostics:\n"); print(x$diagnostics, row.names = FALSE)
  invisible(x)
}

#' Predict nonlinear one-sided PCA scores
#'
#' New observations are mapped to the target grid by the cyclically monotone
#' extension of the empirical center-outward distribution function (see
#' [nlos_distribution()]). Their scores are then computed along the
#' corresponding empirical sign curve, with the observation in place of the grid
#' point; for smoothed paths, the regularized curve of the ray changes
#' accordingly. For the observations used in the fit, this reproduces the
#' fitted scores, with one exception: the \eqn{n_0} observations assigned to
#' the origin are represented in the fit by their mean, and as new data they
#' are mapped to the grid point of their own cell.
#'
#' @param object An `nlospca` object.
#' @param newdata Optional numeric matrix with `object$p` columns.
#' @param n_components Number of directions to use. Defaults to all.
#' @param ... Ignored.
#' @return Matrix of nonlinear one-sided scores (curvilinear projection lengths).
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' fit <- nlospca(X)
#' predict(fit, newdata = cbind(c(-1, 0, 1), c(1, 0, 1)), n_components = 2)
#' @export
predict.nlospca <- function(object, newdata = NULL, n_components = NULL, ...) {
  idx <- .component_index(object, n_components)
  if (is.null(newdata)) return(object$scores[, idx, drop = FALSE])
  X <- .check_matrix(newdata)
  if (ncol(X) != object$p) stop("newdata must have the same number of columns.", call. = FALSE)
  Z <- vapply(idx, function(j) .nlos_project(object, t(object$directions[, j]), X)$score,
              numeric(nrow(X)))
  Z <- matrix(Z, nrow = nrow(X), dimnames = list(rownames(X), colnames(object$directions)[idx]))
  Z
}

#' @rdname scores
#' @export
scores.nlospca <- function(object, n_components = NULL, ...)
  object$scores[, .component_index(object, n_components), drop = FALSE]

#' @rdname directions
#' @export
directions.nlospca <- function(object, n_components = NULL, ...)
  object$directions[, .component_index(object, n_components), drop = FALSE]

#' @rdname nlospca-methods
#' @export
fitted.nlospca <- function(object, n_components = NULL, ...)
  reconstruct(object, n_components = n_components)

#' @rdname nlospca-methods
#' @export
residuals.nlospca <- function(object, n_components = NULL, ...)
  object$X - reconstruct(object, n_components = n_components)

#' @rdname reconstruct
#' @param n_cone For `nlospca` fits: size of the sphere grid whose points in the
#'   positive span of the retained directions are used as candidate directions
#'   (default as in the fit).
#' @param smooth For `nlospca` fits: if `TRUE` (default), the observations are
#'   projected onto the regularized curves (parameter `eps_curves` of the fit),
#'   with correspondingly regularized observation paths; if `FALSE`, they are
#'   projected onto the unsmoothed curves \eqn{\overline{Q}_\pm^{(n)}(r_l u)}
#'   (\eqn{\varepsilon = 0}). With `kappa = 0` and all \eqn{2d} directions, the
#'   unsmoothed reconstruction is exact, whereas the smoothed one carries the
#'   smoothing bias.
#' @details For an `nlospca` fit, each observation is projected onto the
#'   nonlinear manifold given by the image under \eqn{\overline{Q}_\pm^{(n)}} of
#'   the positive span of the first `n_components` directions. The projection is
#'   the point on the curve of the direction in that span with the largest
#'   curvilinear projection score, at an arc length from the center equal to
#'   that score. The argument `add_center` is ignored.
#' @export
reconstruct.nlospca <- function(object, newdata = NULL, n_components = NULL,
                               add_center = TRUE, n_cone = object$settings$n_cone,
                               smooth = TRUE, ...) {
  idx <- .component_index(object, n_components)
  X <- NULL
  if (!is.null(newdata)) {
    X <- .check_matrix(newdata)
    if (ncol(X) != object$p) stop("newdata must have the same number of columns.", call. = FALSE)
  }
  V <- .cone_candidates(object$directions[, idx, drop = FALSE], n_cone,
                        dirs = object$transport$dirs)
  eps <- if (isTRUE(smooth)) object$settings$eps_curves else 0
  Xhat <- .nlos_project(object, V, X, eps = eps)$points
  rownames(Xhat) <- rownames(if (is.null(X)) object$X else X)
  Xhat
}

#' @rdname screeplot.ospca
#' @export
screeplot.nlospca <- function(x, type = c("both", "objective", "mse"), ...) {
  type <- match.arg(type)
  k <- ncol(x$directions)
  obj <- x$objective
  mse <- if (type == "objective") NULL else
    vapply(seq_len(k), function(j)
      reconstruction_mse(x$X, reconstruct(x, n_components = j)), numeric(1))
  labels <- colnames(x$directions)
  ylab <- expression(paste("criterion ", M(u[i])))

  if (type == "objective") {
    graphics::barplot(obj, names.arg = labels, col = "#3b5b92", border = NA,
                      ylab = ylab, main = "Criterion values", ...)
  } else if (type == "mse") {
    plot(seq_len(k), mse, type = "b", pch = 16, col = "#c0392b",
         xlab = "number of directions (k)",
         ylab = "reconstruction MSE", main = "Reconstruction error", ...)
  } else {
    op <- graphics::par(mar = c(5, 4, 4, 4) + 0.3); on.exit(graphics::par(op))
    bp <- graphics::barplot(obj, names.arg = labels, col = "#3b5b92", border = NA,
                            ylab = ylab,
                            main = "Criterion and reconstruction error", ...)
    graphics::par(new = TRUE)
    plot(bp, mse, type = "b", pch = 16, col = "#c0392b", axes = FALSE,
         xlab = "", ylab = "", xlim = graphics::par("usr")[1:2], xaxs = "i",
         ylim = c(0, max(mse)))
    graphics::axis(4, col = "#c0392b", col.axis = "#c0392b")
    graphics::mtext("reconstruction MSE", side = 4, line = 2.5, col = "#c0392b")
  }
  invisible(x)
}

#' Plot a nonlinear one-sided PCA fit
#'
#' For bivariate data, draws the data together with the principal curves, each
#' starting at the center \eqn{Q_\pm^{(n)}(0)} and labeled by its index. For
#' higher dimensions, plots two chosen columns of the scores against each
#' other.
#'
#' @param x An `nlospca` object.
#' @param y Ignored.
#' @param choices For `d > 2`, the two score columns to plot.
#' @param col_curves Colors of the principal curves (recycled).
#' @param ... Passed to the underlying plotting call.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' plot(nlospca(X))
#' @export
plot.nlospca <- function(x, y = NULL, choices = 1:2,
                        col_curves = c("navy", "#c0392b", "#27ae60", "#8e44ad"), ...) {
  if (x$p == 2L) {
    nm <- colnames(x$X); if (is.null(nm)) nm <- c("X1", "X2")
    plot(x$X[, 1], x$X[, 2], xlab = nm[1], ylab = nm[2], pch = 16,
         col = "grey60", asp = 1, ...)
    col_curves <- rep_len(col_curves, length(x$curves))
    for (i in seq_along(x$curves)) {
      C <- x$curves[[i]]
      graphics::lines(C[, 1], C[, 2], col = col_curves[i], lwd = 2)
      graphics::text(C[nrow(C), 1], C[nrow(C), 2], labels = i, col = col_curves[i], pos = 3)
    }
    graphics::points(x$center[1], x$center[2], pch = 19, cex = 1.3)
  } else {
    choices <- as.integer(choices)
    if (length(choices) != 2L) stop("choices must contain exactly two indices.", call. = FALSE)
    Z <- x$scores
    plot(Z[, choices[1]], Z[, choices[2]], xlab = colnames(Z)[choices[1]],
         ylab = colnames(Z)[choices[2]], pch = 16, col = "grey50", ...)
  }
  invisible(x)
}
