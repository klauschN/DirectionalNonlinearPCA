# Internal helpers (not exported).

.check_matrix <- function(X) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.", call. = FALSE)
  if (nrow(X) < 1L || ncol(X) < 1L)
    stop("X must have at least one row and one column.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values.", call. = FALSE)
  X
}

.check_direction <- function(u, p) {
  u <- as.numeric(u)
  if (length(u) != p) stop("u must have length ncol(X).", call. = FALSE)
  if (anyNA(u)) stop("u must not contain missing values.", call. = FALSE)
  nrm <- sqrt(sum(u^2))
  if (!is.finite(nrm) || nrm <= .Machine$double.eps)
    stop("u must be nonzero.", call. = FALSE)
  u / nrm
}

# Directions given as a vector or as a matrix with one direction per row;
# returns the row-normalised matrix.
.check_directions <- function(u, p) {
  U <- if (is.matrix(u)) u else matrix(u, nrow = 1L)
  if (ncol(U) != p) stop("directions must have length ncol(X).", call. = FALSE)
  matrix(t(apply(U, 1L, .check_direction, p = p)), ncol = p)
}

.center_matrix <- function(X, center = TRUE) {
  if (isTRUE(center)) {
    mu <- colMeans(X); Xc <- sweep(X, 2L, mu, FUN = "-")
  } else if (isFALSE(center)) {
    mu <- rep(0, ncol(X)); Xc <- X
  } else {
    mu <- as.numeric(center)
    if (length(mu) != ncol(X))
      stop("center must be TRUE, FALSE, or a numeric vector of length ncol(X).",
           call. = FALSE)
    Xc <- sweep(X, 2L, mu, FUN = "-")
  }
  list(X = Xc, center = mu)
}

# Orthonormal basis of the orthogonal complement of span(basis).
.complement_basis <- function(basis, p) {
  if (is.null(basis) || ncol(basis) == 0L) return(diag(p))
  qr.Q(qr(basis), complete = TRUE)[, (ncol(basis) + 1L):p, drop = FALSE]
}
