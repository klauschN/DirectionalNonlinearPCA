options(ospca.threads = 2L)

parab <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  cbind(x1 = x, x2 = x^2 + rnorm(n, sd = 0.5))
}

test_that("grid sizes respect n = nR * nS + n0 with 0 <= n0 <= min(nR, nS)", {
  for (n in c(50, 301, 1000)) for (d in 2:3) {
    s <- ospca:::.grid_sizes(n, d)
    expect_equal(s$nR * s$nS + s$n0, n)
    expect_true(s$n0 >= 1 && s$n0 <= min(s$nR, s$nS))
    if (d > 2) expect_equal(s$nS %% 2L, 0L)
  }
  expect_error(ospca:::.grid_sizes(100, 2, nR = 3, nS = 20), "n0")
})

test_that("potentials reproduce the coupling and F inverts Q, also with n0 > 1", {
  X <- parab(395)
  fit <- nlospca(X, nR = 15, nS = 26, n_components = 1)   # n0 = 5 origin copies
  tf <- fit$transport
  expect_equal(fit$settings$n0, 5L)
  expect_gt(tf$e0, 0)
  N <- nrow(tf$Gi)
  expect_equal(N, 15L * 26L + 1L)                          # origin copies merged
  expect_equal(tf$center, colMeans(tf$Y[391:395, ]), ignore_attr = TRUE)
  S <- tf$Gi %*% t(tf$Yi) - rep(tf$psi, each = N)
  expect_equal(max.col(S, ties.method = "first"), seq_len(N))
  expect_equal(ospca:::.fbar_index(tf$Yi, tf$Gi, tf$Yi, tf$psi), seq_len(N))
  expect_equal(nlos_distribution(fit, tf$Yi), tf$Gi, ignore_attr = TRUE)
})

test_that("Karp margin is the maximal margin of the potentials", {
  fit <- nlospca(parab(200), n_components = 1)
  tf <- fit$transport
  S <- tf$Gi %*% t(tf$Yi) - rep(tf$psi, each = nrow(tf$Gi))
  m <- diag(S); diag(S) <- -Inf
  expect_equal(min(m - apply(S, 1, max)), tf$e0, tolerance = 1e-6)
})

test_that("Yosida interpolation: exact at the nodes for kappa = 'auto', optimal for larger kappa", {
  fit <- nlospca(parab(200), kappa = "auto", n_components = 1)
  tf <- fit$transport
  expect_gt(fit$settings$kappa, 0)
  expect_equal(nlos_quantile(fit, tf$Gi), tf$Yi, ignore_attr = TRUE)
  ## optimality: w = Qbar(z) is the gradient of the Moreau envelope, i.e. the
  ## prox point p = z - kappa w minimises phi(p) + ||z - p||^2 / (2 kappa)
  kappa <- 0.05
  set.seed(3); Z <- matrix(runif(40, -0.7, 0.7), ncol = 2)
  W <- ospca:::.qbar(Z, tf$Yi, tf$psi, kappa, 2L)
  phi <- function(p) max(drop(tf$Yi %*% p) - tf$psi)
  f <- function(p, z) phi(p) + sum((z - p)^2) / (2 * kappa)
  for (i in seq_len(nrow(Z))) {
    p <- Z[i, ] - kappa * W[i, ]
    pert <- matrix(rnorm(200, sd = 1e-3), ncol = 2)
    expect_true(all(apply(pert, 1, function(e) f(p + e, Z[i, ])) >= f(p, Z[i, ]) - 1e-12))
  }
})

test_that("nlospca returns 2d ordered, orthogonal / antipodal directions", {
  fit <- nlospca(parab())
  expect_s3_class(fit, "nlospca")
  expect_equal(ncol(fit$directions), 4L)
  expect_true(all(diff(fit$objective) <= 1e-12))
  ip <- crossprod(fit$directions)
  expect_true(all(abs(ip) < 1e-8 | abs(abs(ip) - 1) < 1e-8))
  expect_equal(unname(fit$objective), unname(nlos_objective(t(fit$directions), fit)))
  expect_equal(dim(fit$scores), c(nrow(parab()), 4L))
  expect_true(all(fit$scores >= 0))
  # the first direction of a parabola opening upwards points upwards
  expect_gt(fit$directions[2, 1], 0)
})

test_that("principal curves start at the centre and use eps_curves", {
  fit <- nlospca(parab(), eps_curves = 0.1)
  C <- nlos_curve(fit$directions[, 1], fit, eps = 0.1)
  expect_equal(C[1, ], fit$center, ignore_attr = TRUE)
  expect_equal(fit$curves[[1]], C)
})

test_that("predict / reconstruct on the training data reproduce the fit", {
  X <- parab()
  for (paths in c("smooth", "raw")) {
    fit <- nlospca(X, paths = paths)
    # observations assigned to the origin are represented by their mean, so
    # only the others are mapped back to their own grid point by F
    keep <- fit$transport$obs_grid <= fit$settings$nR * fit$settings$nS
    expect_equal(unname(predict(fit, X)[keep, ]), unname(scores(fit)[keep, ]))
    expect_equal(unname(reconstruct(fit, X, n_components = 2)[keep, ]),
                 unname(reconstruct(fit, n_components = 2)[keep, ]))
    # scores of new data are continuous in the observation
    expect_equal(unname(predict(fit, X + 1e-8)[keep, ]), unname(scores(fit)[keep, ]),
                 tolerance = 1e-5)
  }
  fit <- nlospca(X)
  expect_equal(residuals(fit, n_components = 1), X - fitted(fit, n_components = 1))
})

test_that("unsmoothed reconstruction with all 2d directions is exact (kappa = 0)", {
  X <- parab(301)
  fit <- nlospca(X, nR = 10, nS = 30)                       # n0 = 1
  expect_lt(reconstruction_mse(X, reconstruct(fit, smooth = FALSE)), 1e-20)
  m <- sapply(1:4, function(k) reconstruction_mse(X, reconstruct(fit, n_components = k)))
  expect_true(all(diff(m) <= 1e-8))
})

test_that("smoothed observation paths remove spurious scores opposite to U1", {
  X <- parab(400)
  rel_score <- sapply(c("raw", "smooth"), function(paths) {
    fit <- nlospca(X, paths = paths, n_components = 1)
    tf <- fit$transport; m <- tf$nR * tf$nS
    g <- tf$obs_grid; ok <- g <= m
    opp <- drop(tf$dirs[(g[ok] - 1) %/% tf$nR + 1, ] %*% fit$directions[, 1]) < -0.5
    mean(fit$scores[ok, 1][opp]) / sum(sqrt(rowSums(diff(fit$curves[[1]])^2)))
  })
  expect_lt(rel_score[["smooth"]], rel_score[["raw"]] / 2)
})

test_that("projection cone is non-degenerate in d = 3", {
  set.seed(2)
  X3 <- cbind(rnorm(250), rnorm(250)); X3 <- cbind(X3, X3[, 1]^2 + rnorm(250, sd = 0.3))
  fit <- nlospca(X3, n_components = 3)
  V <- ospca:::.cone_candidates(fit$directions[, 1:2], 40, dirs = fit$transport$dirs)
  expect_gt(nrow(V), 2L)
  # every candidate lies in the positive span of the two directions
  coef <- V %*% fit$directions[, 1:2]
  expect_true(all(coef > -1e-6))
  expect_equal(unname(rowSums(coef^2)), rep(1, nrow(V)), tolerance = 1e-6)
})

test_that("input checks", {
  expect_error(nlospca(matrix(rnorm(50), ncol = 1)), "two variables")
  expect_error(nlospca(parab(), eps = -1), "eps")
  expect_error(nlospca(parab(), delta = 1), "delta")
})

test_that("duplicated observations give a zero-margin warning, not a non-optimality one", {
  set.seed(5)
  X <- parab(150)
  X <- X[c(seq_len(150), seq_len(150)), ]                 # every observation twice
  w <- character(0)
  withCallingHandlers(nlospca(X, n_components = 1),
                      warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("zero margin", w)))
  expect_false(any(grepl("not optimal", w)))
})
