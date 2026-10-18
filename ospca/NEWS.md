# ospca 1.0.0

First CRAN release.

* `gspca()`: arguments renamed to match `nlospca()`: `alpha` is now `eps_curves`
  (smoothing of the returned curves), `eps` is now `tol` (eigenvalue threshold
  for the local Jacobians), and `method` is now `ot_method`.
* `nlospca()`: `method` is now `ot_method` (the optimal transport algorithm);
  `method` in `ospca()` still selects the optimizer.
* `nlospca()`: the default size `mS` of the direction grids is now the smallest
  multiple of `nS` that is at least `n` (was `4 * nS`). The finer grid finds
  narrow maxima of the rough empirical criterion that the coarser grid missed,
  at almost no extra cost. The directions, and thus the results, can therefore
  differ from those of earlier versions; `mS = 4 * nS` gives the earlier grid.
* `successive_optimizers()` has been removed; use `directions(ospca(X))`.
* `screeplot()` for `ospca` and `nlospca` fits: the reconstruction errors are now
  drawn above the bar centers (they were shifted).
* The default number of threads also respects the environment variable
  `OMP_THREAD_LIMIT`.
* `screeplot()`: shorter titles and axis labels.
* Examples for all exported functions; references with DOIs in `DESCRIPTION`.
* Three vignettes, `vignette("ospca")`, `vignette("nlospca")` and
  `vignette("gspca")`, illustrate the methods and their helper functions.
* Help pages: the Value sections list the components of the fitted objects,
  each argument is documented separately, and the notation is defined (e.g.,
  the positive part `<u, X>^+ = max(<u, X>, 0)` in `ospca()`), with plain-text
  versions of the main formulas for the help in the console.

# ospca 0.5.0

* **Renamed** (breaking): the nonlinear one-sided PCA is now `nlospca()` (was
  `onlpca()`), its helpers `nlos_objective()`, `nlos_curve()`, `nlos_quantile()`
  and `nlos_distribution()` (were `onl_*()`), its class `"nlospca"` and its
  component labels `NLOSPC1, ...`. The linear and nonlinear one-sided methods are
  now `ospca()` / `nlospca()` with helpers `os_*()` / `nlos_*()`.
* `gspca()`: new estimation choices, with differences over two grid steps and
  a finer grid. New arguments `step` (finite-difference step, default 2) and `boundary`
  (`"clamp"`, default, `"nearest"`, `"trim"`); the default grid has about two
  observations per grid point (was eight). Grid points whose Jacobian estimate
  is not positive definite are left out of `Jbar` (they were floored at `eps`,
  which made `Jbar` meaningless with tied observations); a warning is given
  when more than 5% are left out. New `alpha` for smoothing the returned
  curves. Single-variable data are supported. In a simulation with known `Jbar`
  the new defaults reduced the mean eigenvalue error from 0.048 to 0.036, most
  in three dimensions.

# ospca 0.4.1

* `nlospca()`: a margin of the potentials that is negative only by rounding
  (ties, e.g., duplicated observations) is reported as a zero margin, not as a
  non-optimal assignment.

# ospca 0.4.0

* New `gspca()`: independent nonlinear component analysis of Gunsilius and
  Schennach (2023), with the usual methods.

# ospca 0.3.0

* Compiled potentials (Karp's algorithm) and Yosida interpolation (C++, OpenMP);
  `quadprog` no longer needed.
* Smoothed observation paths (`paths`), separate smoothing of the returned
  curves (`eps_curves`), `reconstruct(smooth = FALSE)`, origin observations
  merged into their mean.

# ospca 0.2.0

* New nonlinear one-sided PCA (measure-transportation-based).
