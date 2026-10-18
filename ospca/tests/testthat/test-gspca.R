options(ospca.threads = 2L)

## data with a known Brenier map: T^{-1}(z) = R g(R'z), g separable and monotone,
## so that Jbar = R diag(E log g'(Z)) R' (eigenvalues 0.686 and 0)
gs_rot <- cbind(c(cos(pi / 6), sin(pi / 6)), c(-sin(pi / 6), cos(pi / 6)))
gs_gen <- function(n, seed = 1) {
  set.seed(seed)
  W <- matrix(rnorm(2 * n), ncol = 2) %*% gs_rot
  cbind(W[, 1] + 0.5 * W[, 1]^3, 2 * exp(W[, 2] / 2)) %*% t(gs_rot)
}
gs_truth <- c(integrate(function(z) dnorm(z) * log(1 + 1.5 * z^2), -Inf, Inf)$value, 0)
angle <- function(u, v) acos(min(1, abs(sum(u * v)))) * 180 / pi

test_that("Gaussian data: linear PCs and half log-eigenvalues (Corollary 1)", {
  set.seed(1)
  Sigma <- matrix(c(4, 1, 1, 1), 2)
  X <- matrix(rnorm(8000), ncol = 2) %*% chol(Sigma)
  fit <- gspca(X)
  ev <- eigen(Sigma)
  expect_lt(angle(fit$directions[, 1], ev$vectors[, 1]), 3)
  expect_lt(max(abs(fit$eigenvalues - 0.5 * log(ev$values))), 0.06)
})

test_that("known nonlinear Brenier map: directions and eigenvalues of Jbar", {
  fit <- gspca(gs_gen(8000))
  expect_lt(angle(fit$directions[, 1], gs_rot[, 1]), 3)
  expect_lt(max(abs(fit$eigenvalues - gs_truth)), 0.06)
})

test_that("gspca object structure", {
  X <- gs_gen(2000)
  fit <- gspca(X, n_components = 1)
  expect_s3_class(fit, "gspca")
  expect_equal(dim(fit$directions), c(2L, 1L))
  expect_length(fit$eigenvalues, 2L)
  expect_true(all(diff(fit$eigenvalues) <= 0))
  expect_equal(unname(fit$scores), unname(fit$TX %*% fit$directions))
  expect_equal(unname(scores(fit)), unname(predict(fit)))
  expect_length(fit$curves, 1L)
})

test_that("T^{-1} interpolates the discrete map, and F inverts it on the nodes", {
  fit <- gspca(gs_gen(2000))
  gr <- fit$grid
  expect_equal(unname(ospca:::.gs_tinv(fit, gr$G)), gr$Yg)
  expect_equal(unname(predict(fit, gr$Yg)), unname(gr$G %*% fit$directions))
})

test_that("principal curves pass through T^{-1}(0)", {
  fit <- gspca(gs_gen(2000))
  C <- fit$curves[[1]]
  expect_equal(unname(C[(nrow(C) + 1) / 2, ]),
               unname(ospca:::.gs_tinv(fit, matrix(0, 1, 2))[1, ]))
})

test_that("reconstruction error decreases with the number of directions", {
  X <- gs_gen(2000)
  fit <- gspca(X)
  m <- sapply(1:2, function(k) reconstruction_mse(X, reconstruct(fit, n_components = k)))
  expect_lt(m[2], m[1])
  expect_equal(residuals(fit, n_components = 1), X - fitted(fit, n_components = 1))
})

test_that("step and boundary options, and a single variable", {
  X <- gs_gen(2000)
  for (b in c("nearest", "clamp", "trim")) {
    fit <- gspca(X, n_grid = 20, step = 2, boundary = b)
    expect_equal(fit$settings$boundary, b)
    expect_lt(angle(fit$directions[, 1], gs_rot[, 1]), 5)
    expect_true(fit$diagnostics$dropped >= 0 && fit$diagnostics$dropped < 0.1)
  }
  expect_error(gspca(X, n_grid = 4, step = 2, boundary = "nearest"), "n_grid")
  set.seed(4)
  z <- rnorm(500)
  f1 <- gspca(cbind(2 * z + 1), n_grid = 25)                # linear map: log 2
  expect_equal(f1$eigenvalues, log(2), tolerance = 0.05)
})

test_that("duplicated observations: degenerate local Jacobians are left out", {
  set.seed(6)
  X0 <- gs_gen(300)
  X <- X0[sample(300, 2000, replace = TRUE), ]
  fit <- suppressWarnings(gspca(X))
  expect_gt(fit$diagnostics$dropped, 0)
  expect_true(all(is.finite(fit$eigenvalues)))
  expect_lt(max(abs(fit$eigenvalues - gs_truth)), 0.3)      # ties bias, but no breakdown
})

test_that("argument names shared with nlospca", {
  X <- gs_gen(1000)
  f0 <- gspca(X)
  f1 <- gspca(X, eps_curves = 2, tol = 1e-6, ot_method = "networkflow")
  expect_equal(f1$settings[c("eps_curves", "tol", "ot_method")],
               list(eps_curves = 2, tol = 1e-6, ot_method = "networkflow"))
  ## smoothing changes the displayed curves only, not the estimate
  expect_equal(f1$eigenvalues, f0$eigenvalues)
  expect_false(isTRUE(all.equal(f1$curves[[1]], f0$curves[[1]])))
  expect_error(gspca(X, eps_curves = -1), "eps_curves")
})

test_that("the default number of threads respects OMP_THREAD_LIMIT", {
  old_opt <- options(ospca.threads = NULL)
  old_env <- Sys.getenv(c("OMP_THREAD_LIMIT", "_R_CHECK_LIMIT_CORES_"), unset = NA)
  on.exit({
    options(old_opt)
    for (v in names(old_env))
      if (is.na(old_env[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old_env[v]))
  })
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  Sys.setenv(OMP_THREAD_LIMIT = "1")
  expect_equal(ospca:::.default_threads(), 1L)
  Sys.unsetenv("OMP_THREAD_LIMIT")
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_lte(ospca:::.default_threads(), 2L)
  options(ospca.threads = 3L)
  Sys.setenv(OMP_THREAD_LIMIT = "2")
  expect_equal(ospca:::.default_threads(), 2L)
})

test_that("input checks", {
  expect_error(gspca(gs_gen(200), n_grid = 20), "n_grid")
  expect_error(gspca(gs_gen(200), n_grid = 2), "n_grid")
  expect_error(gspca(gs_gen(200), n_components = 3), "n_components")
})
