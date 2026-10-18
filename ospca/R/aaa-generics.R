# Generics and help topics shared by ospca, nlospca and gspca fits. The file name
# makes roxygen read these blocks first, so that the shared help topics take
# their names and texts from here (the methods add their usage with @rdname).

#' Extract scores
#'
#' Scores of a fitted [ospca()] (positive-part projections), [nlospca()]
#' (curvilinear projection lengths) or [gspca()] (coordinates
#' \eqn{\langle T(X), u_j \rangle}) model.
#'
#' @param object A fitted `ospca`, `nlospca` or `gspca` object.
#' @param n_components Number of score columns. Defaults to all.
#' @param ... Ignored.
#' @return Matrix of scores.
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' head(scores(ospca(X), n_components = 2))
#' head(scores(nlospca(X), n_components = 2))
#' head(scores(gspca(X)))
#' @name scores
#' @export
scores <- function(object, n_components = NULL, ...) UseMethod("scores")

#' Extract principal directions
#'
#' Principal directions (unit vectors, in columns) of a fitted [ospca()],
#' [nlospca()] or [gspca()] model, in decreasing order of the criterion.
#'
#' @param object A fitted `ospca`, `nlospca` or `gspca` object.
#' @param n_components Number of directions. Defaults to all.
#' @param ... Ignored.
#' @return Matrix of directions, one per column.
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' directions(ospca(X))
#' directions(nlospca(X), n_components = 2)
#' directions(gspca(X))
#' @name directions
#' @export
directions <- function(object, n_components = NULL, ...) UseMethod("directions")

#' Reconstruct data from a fitted model
#'
#' Dimension reduction with the first `n_components` directions. For an `ospca`
#' fit, this is the linear directional reconstruction onto the positive span of
#' the directions, \eqn{\sum_{i=1}^k P_i^+ U_i^+}; for `nlospca` and `gspca`
#' fits, see Details.
#'
#' @param object A fitted `ospca`, `nlospca` or `gspca` object.
#' @param newdata Optional data matrix to reconstruct.
#' @param n_components Number of directions to use. Defaults to all.
#' @param add_center For `ospca` fits: whether to add the fitted center back
#'   (ignored otherwise).
#' @param ... Ignored.
#' @return Matrix of reconstructed observations.
#' @examples
#' set.seed(1)
#' x <- rnorm(300)
#' X <- cbind(x, x^2 + rnorm(300, sd = 0.5))
#' fit_lin <- ospca(X)
#' fit_nl <- nlospca(X)
#' fit_gs <- gspca(X)
#' ## reconstruction errors with one direction
#' reconstruction_mse(X, reconstruct(fit_lin, n_components = 1))
#' reconstruction_mse(X, reconstruct(fit_nl, n_components = 1))
#' reconstruction_mse(X, reconstruct(fit_gs, n_components = 1))
#' ## nlospca: unsmoothed curves reproduce the data with all 2d directions,
#' ## apart from the observations assigned to the origin
#' reconstruction_mse(X, reconstruct(fit_nl, smooth = FALSE))
#' @name reconstruct
#' @export
reconstruct <- function(object, newdata = NULL, n_components = NULL,
                        add_center = TRUE, ...) UseMethod("reconstruct")

#' Scree-type plot
#'
#' Shows how much each direction contributes. With `type = "objective"` the
#' criterion values are drawn as bars (the one-sided objective \eqn{J(u_i)}
#' for `ospca`, the expected squared curve length for `nlospca` and the
#' eigenvalues of \eqn{\bar J} for `gspca`). With `type = "mse"` the
#' reconstruction error is drawn as a function of the number of retained
#' directions, and `type = "both"` (default) shows the two together, with the
#' reconstruction error on a second axis.
#'
#' @param x A fitted `ospca`, `nlospca` or `gspca` object.
#' @param type One of `"both"`, `"objective"`, `"mse"`.
#' @param ... Passed to the underlying plotting calls.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' fit <- ospca(rbind(matrix(rnorm(800), ncol = 2) %*% diag(c(3, 1)),
#'                     cbind(abs(rnorm(250, 5)), rnorm(250))))
#' screeplot(fit, type = "both")
#' @name screeplot.ospca
NULL
