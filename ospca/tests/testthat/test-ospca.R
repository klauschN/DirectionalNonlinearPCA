asym2d <- function(seed = 1) {
  set.seed(seed)
  rbind(matrix(rnorm(2000), ncol = 2) %*% diag(c(3, 1)),
        cbind(abs(rnorm(600, 5)), rnorm(600)))
}

test_that("power and bfgs optimisers agree on the optimum", {
  X <- asym2d(11)
  fp <- ospca(X, n_starts = 120, n_polish = 10, method = "power")
  fb <- ospca(X, n_starts = 120, n_polish = 10, method = "bfgs")
  # top one-sided objective should match closely
  expect_lt(abs(fp$objective[1] - fb$objective[1]), 1e-4 * fp$objective[1])
  # leading directions equal up to sign
  expect_lt(min(sqrt(colSums((fp$directions[, 1] - outer(fb$directions[, 1], c(1, -1)))^2))), 1e-3)
})

test_that("power iteration is the default method", {
  expect_identical(formals(ospca)$method, formals(os_axes)$method)
  expect_equal(eval(formals(os_axes)$method)[1], "power")
})

test_that("ospca returns 2d ordered directions and diagnostics", {
  fit <- ospca(asym2d(), n_starts = 80, n_polish = 5)
  expect_s3_class(fit, "ospca")
  expect_equal(ncol(fit$directions), 4L)
  expect_true(all(diff(fit$objective) <= 1e-8))
  expect_true(all(c("direction", "objective", "gap", "tie") %in% names(fit$diagnostics)))
})

test_that("os_axes is orthonormal and carries diagnostics", {
  X <- scale(asym2d(2), scale = FALSE)
  A <- os_axes(X, n_starts = 80, n_polish = 5)
  expect_lt(max(abs(crossprod(A) - diag(2))), 1e-6)
  expect_s3_class(attr(A, "diagnostics"), "data.frame")
})

test_that("antipodal directions differ in objective (one-sidedness)", {
  fit <- ospca(asym2d(3), n_starts = 80, n_polish = 5)
  j <- which(apply(fit$directions, 2, function(u) sum(u * fit$directions[, 1])) < -0.99)
  expect_length(j, 1L)
  expect_gt(abs(fit$objective[1] - fit$objective[j]), 1e-3)
})

test_that("reconstruction MSE is monotone and zero at k = 2d", {
  fit <- ospca(asym2d(4), n_starts = 80, n_polish = 5)
  m <- sapply(1:4, function(k)
    reconstruction_mse(fit$X_centered, reconstruct(fit, n_components = k, add_center = FALSE)))
  expect_true(all(diff(m) <= 1e-8))
  expect_lt(m[4], 1e-8)
})

test_that("a genuine rotational multiplicity triggers a warning", {
  set.seed(10)
  # Build data with exact 4-fold rotational symmetry: the empirical one-sided
  # objective is then invariant under 90-degree rotation, so the maximiser is
  # non-unique (a and its 90-degree rotation tie).
  B <- matrix(rnorm(2000), ncol = 2) %*% diag(c(3, 1))
  R <- matrix(c(0, 1, -1, 0), 2)            # 90-degree rotation
  X4 <- rbind(B, B %*% t(R), B %*% t(R %*% R), B %*% t(R %*% R %*% R))
  expect_warning(ospca(X4, n_starts = 150, n_polish = 10),
                 "tie|non-unique|multiplicity|arbitrary")
})

test_that("asymmetric data is decisive (large gap for the first direction)", {
  fit <- suppressWarnings(ospca(asym2d(7), n_starts = 120, n_polish = 10))
  expect_false(isTRUE(fit$diagnostics$tie[1]))
  expect_gt(fit$diagnostics$gap[1], 0.05)
})

test_that("robust optimiser beats a careless single start on objective", {
  X <- scale(asym2d(5), scale = FALSE)
  fit <- ospca(X, n_starts = 120, n_polish = 12)
  # the top objective should match the best axis-aligned guess or exceed it
  cand <- max(os_objective(c(1, 0), X), os_objective(c(0, 1), X),
              os_objective(c(-1, 0), X), os_objective(c(0, -1), X))
  expect_gte(fit$objective[1] + 1e-8, cand)
})

test_that("screeplot runs for all types", {
  fit <- ospca(asym2d(6), n_starts = 60, n_polish = 4)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_invisible(screeplot(fit, type = "objective"))
  expect_invisible(screeplot(fit, type = "mse"))
  expect_invisible(screeplot(fit, type = "both"))
  dev.off(); unlink(tmp)
})
