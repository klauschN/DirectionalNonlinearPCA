#' One-sided projection scores
#'
#' Positive parts
#' \deqn{\langle u, x_i \rangle^+ = \max(\langle u, x_i \rangle, 0)}{<u, x_i>^+ = max(<u, x_i>, 0)}
#' of the projections of the rows \eqn{x_i} of `X` onto a direction \eqn{u},
#' computed as `pmax(X %*% u, 0)`. An observation on the side of \eqn{u} thus
#' gets its projection and an observation on the opposite side gets zero. The
#' direction is normalized internally. Note that `X` is not centered by this
#' function; if needed, center it beforehand.
#'
#' @param u Numeric direction vector of length `ncol(X)`.
#' @param X Numeric matrix with observations in rows.
#' @return Numeric vector of one-sided projection scores.
#' @seealso [os_objective()], [ospca()]
#' @examples
#' os_project(c(1, 0), rbind(c(1, 0), c(-1, 0), c(0, 2)))
#' @export
os_project <- function(u, X) {
  X <- .check_matrix(X)
  u <- .check_direction(u, ncol(X))
  as.numeric(pmax(drop(X %*% u), 0))
}

#' One-sided PCA objective
#'
#' The one-sided objective
#' \deqn{J(u) = E[(\langle u, X \rangle^+)^2],}{J(u) = E[(<u, X>^+)^2],}
#' where \eqn{\langle u, X \rangle^+ = \max(\langle u, X \rangle, 0)}{<u, X>^+ = max(<u, X>, 0)}
#' is the positive part of the projection (see [os_project()]), estimated by
#' the average
#' \deqn{J_n(u) = \frac{1}{n} \sum_{i=1}^n (\langle u, x_i \rangle^+)^2}{J_n(u) = (1/n) sum_i (<u, x_i>^+)^2}
#' over the rows \eqn{x_i} of `X`, that is, by `mean(pmax(X %*% u, 0)^2)`. In
#' [ospca()] the data are centered first; this function uses `X` as it is.
#'
#' @inheritParams os_project
#' @return A single numeric value.
#' @seealso [os_project()], [os_axes()], [ospca()]
#' @examples
#' os_objective(c(1, 0), rbind(c(1, 0), c(-1, 0), c(0, 2)))
#' @export
os_objective <- function(u, X) {
  z <- os_project(u, X)
  mean(z^2)
}

# Internal: objective on data assumed already centred.
.os_obj <- function(u, X) {
  u <- u / sqrt(sum(u^2))
  z <- pmax(as.numeric(X %*% u), 0)
  mean(z^2)
}
