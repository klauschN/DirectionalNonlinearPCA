#' Methods for `ospca` objects
#'
#' Printing, fitted values and residuals of a fitted [ospca()] model. The fitted
#' values are the reconstructions with the first `n_components` directions (see
#' [reconstruct()]), and the residuals are the centered data minus their
#' reconstructions.
#'
#' @param x,object An `ospca` object.
#' @param n_components Number of directions to use where applicable.
#' @param ... Ignored.
#' @return The `print` method returns its argument invisibly; `fitted` and
#'   `residuals` return numeric matrices.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(400), ncol = 2) %*% diag(c(3, 1)),
#'            cbind(abs(rnorm(100, 5)), rnorm(100)))
#' fit <- ospca(X)
#' fit
#' head(fitted(fit, n_components = 2))
#' head(residuals(fit, n_components = 2))
#' @name ospca-methods
#' @rdname ospca-methods
#' @export
print.ospca <- function(x, ...) {
  cat("One-sided (directional) PCA fit\n")
  cat("  Observations:        ", x$n, "\n", sep = "")
  cat("  Variables (d):       ", x$p, "\n", sep = "")
  cat("  Directional dirs (k):", ncol(x$directions), " (max 2d = ", 2L * x$p, ")\n", sep = "")
  cat("  Objective (decreasing):\n    ",
      paste(sprintf("%.4g", x$objective), collapse = ", "), "\n", sep = "")
  if (!is.null(x$diagnostics) && any(x$diagnostics$tie, na.rm = TRUE))
    cat("  Note: near-tied axis/axes detected (see $diagnostics).\n")
  invisible(x)
}

#' Summarize a one-sided PCA fit
#'
#' @param object An `ospca` object.
#' @param ... Ignored.
#' @return An object of class `summary.ospca`.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(400), ncol = 2) %*% diag(c(3, 1)),
#'            cbind(abs(rnorm(100, 5)), rnorm(100)))
#' summary(ospca(X))
#' @export
summary.ospca <- function(object, ...) {
  out <- list(n = object$n, p = object$p, n_components = ncol(object$directions),
              objective = object$objective, directions = object$directions,
              diagnostics = object$diagnostics)
  class(out) <- "summary.ospca"
  out
}

#' @rdname summary.ospca
#' @param x A `summary.ospca` object.
#' @export
print.summary.ospca <- function(x, ...) {
  cat("One-sided (directional) PCA summary\n")
  cat("Dimensions: ", x$n, " x ", x$p, "   directional directions: ",
      x$n_components, "\n\n", sep = "")
  tab <- rbind(objective = x$objective); colnames(tab) <- colnames(x$directions)
  cat("Objective values:\n"); print(round(tab, 4))
  cat("\nAxis diagnostics:\n"); print(x$diagnostics, row.names = FALSE)
  invisible(x)
}

.component_index <- function(object, n_components = NULL) {
  k <- ncol(object$directions)
  if (is.null(n_components)) return(seq_len(k))
  seq_len(min(as.integer(n_components), k))
}

#' Predict one-sided PCA scores
#'
#' One-sided scores of new observations, that is, the positive parts of the
#' projections of the observations, centered with the center of the fit, onto
#' the directions.
#'
#' @param object An `ospca` object.
#' @param newdata Optional numeric matrix with `object$p` columns.
#' @param n_components Number of directions to use. Defaults to all.
#' @param ... Ignored.
#' @return Matrix of one-sided scores.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(400), ncol = 2) %*% diag(c(3, 1)),
#'            cbind(abs(rnorm(100, 5)), rnorm(100)))
#' fit <- ospca(X[1:250, ])
#' predict(fit, newdata = X[251:255, ], n_components = 2)
#' @export
predict.ospca <- function(object, newdata = NULL, n_components = NULL, ...) {
  idx <- .component_index(object, n_components)
  if (is.null(newdata)) return(object$scores[, idx, drop = FALSE])
  X <- .check_matrix(newdata)
  if (ncol(X) != object$p) stop("newdata must have the same number of columns.", call. = FALSE)
  Z <- pmax(sweep(X, 2L, object$center, "-") %*% object$directions[, idx, drop = FALSE], 0)
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
  colnames(Z) <- colnames(object$directions)[idx]
  Z
}

#' @rdname scores
#' @export
scores.ospca <- function(object, n_components = NULL, ...)
  object$scores[, .component_index(object, n_components), drop = FALSE]

#' @rdname directions
#' @export
directions.ospca <- function(object, n_components = NULL, ...)
  object$directions[, .component_index(object, n_components), drop = FALSE]

#' @rdname ospca-methods
#' @export
fitted.ospca <- function(object, n_components = NULL, ...)
  reconstruct(object, n_components = n_components)

#' @rdname ospca-methods
#' @export
residuals.ospca <- function(object, n_components = NULL, ...)
  object$X_centered - reconstruct(object, n_components = n_components, add_center = FALSE)

#' @rdname reconstruct
#' @export
reconstruct.ospca <- function(object, newdata = NULL, n_components = NULL,
                               add_center = TRUE, ...) {
  idx <- .component_index(object, n_components)
  U <- object$directions[, idx, drop = FALSE]
  Z <- if (is.null(newdata)) object$scores[, idx, drop = FALSE]
       else predict(object, newdata, n_components = length(idx))
  Xhat <- Z %*% t(U)
  if (isTRUE(add_center)) Xhat <- sweep(Xhat, 2L, object$center, "+")
  colnames(Xhat) <- rownames(U)
  Xhat
}

#' Mean squared reconstruction error
#'
#' The mean, over the observations, of the squared Euclidean distance between
#' an observation and its reconstruction, `sum((Y - X)^2) / nrow(X)`.
#'
#' @param X Original matrix.
#' @param Y Reconstructed matrix of the same dimension.
#' @return A single numeric value.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(400), ncol = 2) %*% diag(c(3, 1)),
#'            cbind(abs(rnorm(100, 5)), rnorm(100)))
#' fit <- ospca(X)
#' sapply(1:4, function(k) reconstruction_mse(X, reconstruct(fit, n_components = k)))
#' @export
reconstruction_mse <- function(X, Y) {
  X <- .check_matrix(X); Y <- .check_matrix(Y)
  if (!all(dim(X) == dim(Y))) stop("X and Y must have the same dimensions.", call. = FALSE)
  sum((Y - X)^2) / nrow(X)
}

#' @rdname screeplot.ospca
#' @export
screeplot.ospca <- function(x, type = c("both", "objective", "mse"), ...) {
  type <- match.arg(type)
  k <- ncol(x$directions)
  obj <- x$objective
  mse <- vapply(seq_len(k), function(j)
    reconstruction_mse(x$X_centered, reconstruct(x, n_components = j, add_center = FALSE)),
    numeric(1))
  labels <- colnames(x$directions)

  if (type == "objective") {
    graphics::barplot(obj, names.arg = labels, col = "#3b5b92", border = NA,
                      ylab = expression(paste("one-sided objective ", J(u[i]))),
                      main = "One-sided objectives", ...)
  } else if (type == "mse") {
    plot(seq_len(k), mse, type = "b", pch = 16, col = "#c0392b",
         xlab = "number of directions (k)",
         ylab = "reconstruction MSE", main = "Reconstruction error", ...)
  } else {
    op <- graphics::par(mar = c(5, 4, 4, 4) + 0.3); on.exit(graphics::par(op))
    bp <- graphics::barplot(obj, names.arg = labels, col = "#3b5b92", border = NA,
                            ylab = expression(paste("one-sided objective ", J(u[i]))),
                            main = "Objectives and reconstruction error", ...)
    graphics::par(new = TRUE)
    plot(bp, mse, type = "b", pch = 16, col = "#c0392b", axes = FALSE,
         xlab = "", ylab = "", xlim = graphics::par("usr")[1:2], xaxs = "i",
         ylim = c(0, max(mse)))
    graphics::axis(4, col = "#c0392b", col.axis = "#c0392b")
    graphics::mtext("reconstruction MSE", side = 4, line = 2.5, col = "#c0392b")
  }
  invisible(x)
}

#' Plot a one-sided PCA fit
#'
#' For bivariate data, draws the centered data with an arrow for each direction,
#' its length proportional to the square root of the objective. For higher
#' dimensions, plots two chosen columns of the one-sided scores against each
#' other.
#'
#' @param x An `ospca` object.
#' @param y Ignored.
#' @param choices For `d > 2`, the two score columns to plot.
#' @param scale Arrow scaling factor (2D only).
#' @param ... Passed to the underlying plotting call.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(400), ncol = 2) %*% diag(c(3, 1)),
#'            cbind(abs(rnorm(100, 5)), rnorm(100)))
#' plot(ospca(X))
#' @export
plot.ospca <- function(x, y = NULL, choices = 1:2, scale = 1, ...) {
  if (x$p == 2L) {
    Xc <- x$X_centered
    plot(Xc[, 1], Xc[, 2], xlab = "X1 (centered)", ylab = "X2 (centered)",
         pch = 16, col = "grey60", asp = 1, ...)
    graphics::abline(h = 0, v = 0, col = "grey85")
    len <- sqrt(x$objective); len <- len / max(len) * scale * 0.6 * max(abs(Xc))
    for (i in seq_len(ncol(x$directions))) {
      u <- x$directions[, i] * len[i]
      graphics::arrows(0, 0, u[1], u[2], col = "navy", lwd = 2, length = 0.1)
      graphics::text(u[1], u[2], labels = i, col = "navy", pos = 3)
    }
  } else {
    choices <- as.integer(choices)
    if (length(choices) != 2L) stop("choices must contain exactly two indices.", call. = FALSE)
    Z <- x$scores
    plot(Z[, choices[1]], Z[, choices[2]], xlab = colnames(Z)[choices[1]],
         ylab = colnames(Z)[choices[2]], pch = 16, col = "grey50", ...)
  }
  invisible(x)
}
