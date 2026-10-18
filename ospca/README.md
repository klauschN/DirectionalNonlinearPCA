# ospca

Principal component analysis (PCA) is a bidirectional method, as a direction and
its negative describe the same component. Many data sets, however, are
asymmetric, and the two half-lines of a direction can carry quite different
information. One-sided (directional) PCA therefore considers each half-line
separately, which gives up to `2d` directions in dimension `d`. The package
`ospca` contains a linear and a nonlinear version of one-sided PCA and, for
comparison, the independent nonlinear component analysis of Gunsilius and
Schennach (2023):

- `ospca()`: linear one-sided PCA, based on the positive parts of the
  projections;
- `nlospca()`: nonlinear one-sided PCA based on measure transportation, in which
  the half-lines are replaced by the curves traced by the empirical
  center-outward quantile function;
- `gspca()`: the (bidirectional) independent nonlinear component analysis of
  Gunsilius and Schennach (2023).

For all three methods, `directions()`, `scores()`, `predict()`, `reconstruct()`,
`reconstruction_mse()`, `screeplot()` and `plot()` give the directions, the
scores of the observed and of new data, the reconstructions and their errors,
and graphical summaries. The vignettes `vignette("ospca")`,
`vignette("nlospca")` and `vignette("gspca")` illustrate the three methods and
their helper functions on simple data.

## Installation

The package contains compiled C++ code, so installing it from source needs a
C++ compiler (Rtools on Windows).

```r
install.packages("ospca")
```

## Linear one-sided PCA

The importance of a direction `u` is measured by the one-sided objective
`J(u) = E[(<u, X>^+)^2]`, the second moment of the positive part
`<u, X>^+ = max(<u, X>, 0)` of the projection of the centered data. Since
`<u, X>^2 = (<u, X>^+)^2 + (<-u, X>^+)^2`, the variance criterion of classical
PCA is `J(u) + J(-u)`, which treats the two half-lines of a direction alike. `ospca()` first computes greedy one-sided orthonormal axes of the
centered data and then keeps, among the `2d` directions `+/- a_1, ..., +/- a_d`,
those with the largest objective.

Since `J(u)` is convex in `u` but maximized over the unit sphere, it can have
several local maxima. Each axis is therefore found by a multi-start search, in
which the best starting directions are refined by a generalized power iteration,

```
w <- normalize(t(Y) %*% pmax(Y %*% w, 0))      # Y = X projected to the search subspace
```

which increases the objective at every step. A quasi-Newton optimizer is
available as an alternative (`method = "bfgs"`). When the best and the best
clearly different direction (more than 45 degrees away) have nearly the same
objective, the maximizer is effectively not unique; such near-ties are recorded
in `fit$diagnostics` and reported with a warning.

```r
library(ospca)

set.seed(1)
X <- rbind(matrix(rnorm(1000), ncol = 2) %*% diag(c(3, 1)),
           cbind(abs(rnorm(300, 5)), rnorm(300)))

fit <- ospca(X)
fit
summary(fit)                        # objectives and axis diagnostics
screeplot(fit)                      # objectives and reconstruction errors
reconstruction_mse(X, reconstruct(fit, n_components = 2))
```

The arguments `n_starts` and `n_polish` give the number of random starting
directions and the number of them that are refined, and `maxit` and `tol`
control the power iteration.

## Nonlinear one-sided PCA

In `nlospca()`, the sample is optimally assigned to a regular grid in the unit
ball, which gives the empirical center-outward quantile function, and this
function is interpolated in a cyclically monotone way. Then

- the principal directions maximize the expected squared length of the
  (Tikhonov regularized) sign curves `r -> Qbar(r u)`, greedily over the
  orthogonal complement of the previous directions and their antipodes;
- the principal curves are the images of the corresponding radii;
- the principal components are curvilinear projections, that is, the sums of
  the positive parts of the increments of the observation's own (by default
  regularized) sign curve, projected onto the tangents of the principal curve;
- `reconstruct(fit, n_components = k)` projects onto the manifold spanned by the
  first `k` principal curves. With `smooth = FALSE`, it uses the unsmoothed
  curves and, with all `2d` directions and `kappa = 0`, reproduces the data.

```r
set.seed(1)
x <- rnorm(400)
X <- cbind(x, x^2 + rnorm(400, sd = 0.5))

fit <- nlospca(X)                   # the grid is chosen automatically
summary(fit)
plot(fit)                           # data with the principal curves
screeplot(fit)                      # criterion and reconstruction errors
predict(fit, X[1:5, ])              # scores of (new) observations
nlos_curve(directions(fit)[, 1], fit)   # the first principal curve
```

The main tuning parameters are `nR` and `nS` (the target grid), `mS` (the size
of the direction grids, by default the smallest multiple of `nS` that is at
least `n`), `eps` (the Tikhonov smoothing of
the curves that define the criterion), `eps_curves` (the smoothing of the
returned curves and projections, by default `eps`), `paths` (smoothed or raw
observation paths), `delta` (trimming) and `kappa` (Yosida regularization of
the interpolation, with `0` for a piecewise constant interpolation).

The optimal assignment is computed with the exact network simplex algorithm of
the package `transport`. The cyclically monotone potentials (Karp's
minimum-mean-cycle algorithm, with `O(n^3)` time and, up to `n = 6000`, `O(n^2)`
memory) and the Yosida interpolation are compiled C++ code, parallelized with
OpenMP. The number of threads is set by `nthreads` or by
`options(ospca.threads = )`. The computing time is dominated by the potentials
and thus grows with the cube of `n`: on a laptop with 15 threads, a complete fit
took between 0.5 and 1.2 seconds for `n = 1000` and between 15 and 60 seconds
for `n = 4000` (`d = 2`), depending on the state of the machine.

## GS-PCA (Gunsilius and Schennach, 2023)

For comparison, `gspca()` implements the independent nonlinear component
analysis of Gunsilius and Schennach (2023). With `T` the Brenier map that pushes
the data forward to the standard normal distribution, the directions are the
eigenvectors of `Jbar = E_Phi[log J_{T^{-1}}(x)]`, the average matrix logarithm
of the Jacobian of `T^{-1}`. The curves are `T^{-1}(r u_j)`, the scores are
`<T(X), u_j>`, and the `k`-dimensional representation is
`T^{-1}(sum_{j <= k} <T(X), u_j> u_j)` (see `reconstruct()`). For Gaussian data,
the method reduces to linear PCA.

The implementation is a simple discrete estimator and not the smooth potential
fit of the paper: an exact optimal coupling between a tensor grid of normal
quantiles and the sample, barycentric images of the grid points, and
symmetrized central difference Jacobians, where grid points whose Jacobian
estimate is not positive definite are left out. By default, the
differences are taken over two grid steps (`step = 2`), the grid has about two
observations per grid point, and the differences are one-sided at the edge of
the grid (`boundary = "clamp"`); these choices were the most accurate in a small
simulation with known `Jbar`.

```r
fit <- gspca(X)                     # about two observations per grid point
fit
plot(fit)
screeplot(fit)
reconstruct(fit, n_components = 1)
gspca(X, step = 1, boundary = "nearest", n_grid = 20)   # other estimation choices
```

## References

Gunsilius, F. and Schennach, S. (2023). Independent nonlinear component
analysis. *Journal of the American Statistical Association*, 118(542),
1305–1318. <doi:10.1080/01621459.2021.1990768>

Hallin, M., del Barrio, E., Cuesta-Albertos, J. and Matrán, C. (2021).
Distribution and quantile functions, ranks and signs in dimension d: a measure
transportation approach. *Annals of Statistics* 49, 1139–1165.
<doi:10.1214/20-AOS1996>

Hallin, M., Nordhausen, K. and Radojicic, U. Directional nonlinear principal
components: a measure transportation approach.
