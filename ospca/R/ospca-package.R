#' ospca: One-Sided Linear and Nonlinear Principal Component Analysis
#'
#' Principal component analysis (PCA) is a bidirectional method, as a direction
#' and its negative describe the same component. Many data sets, however, are
#' asymmetric, and the two half-lines of a direction can carry quite different
#' information. One-sided (directional) PCA therefore considers each half-line
#' separately, which gives up to \eqn{2d} directions in dimension \eqn{d}. The
#' package contains a linear and a nonlinear version of one-sided PCA and, for
#' comparison, the independent nonlinear component analysis of Gunsilius and
#' Schennach (2023).
#'
#' @section Main functions:
#' * [ospca()] fits the linear one-sided PCA, based on the positive parts of
#'   the projections. The underlying greedy orthonormal axes are computed by
#'   [os_axes()], and [os_objective()] and [os_project()] give the objective and
#'   the one-sided projections for a given direction.
#' * [nlospca()] fits the nonlinear one-sided PCA based on measure
#'   transportation. Here the half-lines are replaced by the curves traced by
#'   the empirical center-outward quantile function, the directions maximize the
#'   expected squared length of these curves, and the components are
#'   curvilinear projections onto them. The functions [nlos_objective()],
#'   [nlos_curve()], [nlos_quantile()] and [nlos_distribution()] give access to
#'   the criterion, the curves and the fitted quantile and distribution
#'   functions.
#' * [gspca()] fits the (bidirectional) independent nonlinear component
#'   analysis of Gunsilius and Schennach (2023), which is based on the optimal
#'   transport to the standard normal distribution.
#'
#' For all three methods, [directions()], [scores()], [predict()] and
#' [reconstruct()] give the directions, the scores of the observed and of new
#' data, and the reconstructions with a chosen number of directions;
#' [reconstruction_mse()] gives the corresponding reconstruction error, and
#' [screeplot()] and [plot()] give graphical summaries. The three methods are
#' illustrated on simple data in the vignettes `vignette("ospca")`,
#' `vignette("nlospca")` and `vignette("gspca")`.
#'
#' @section Ties:
#' When the best direction and the best clearly different direction have
#' nearly the same criterion value, the maximizer is effectively not unique,
#' similar to tied eigenvalues in classical PCA. Such near-ties are recorded in
#' the fit and reported with a warning.
#'
#' @section Compiled code:
#' The potentials of the cyclically monotone interpolation (Karp's
#' minimum-mean-cycle algorithm) and its Yosida regularization are computed in
#' C++ and parallelized with OpenMP. The number of threads is set by the
#' `nthreads` argument of [nlospca()] and [gspca()] or by the option
#' `ospca.threads`; the environment variable `OMP_THREAD_LIMIT` is respected.
#'
#' @references Gunsilius, F. and Schennach, S. (2023). Independent nonlinear
#'   component analysis. *Journal of the American Statistical Association*,
#'   118(542), 1305--1318. \doi{10.1080/01621459.2021.1990768}
#'
#'   Hallin, M., del Barrio, E., Cuesta-Albertos, J. and Matrán, C. (2021).
#'   Distribution and quantile functions, ranks and signs in dimension d: a
#'   measure transportation approach. *Annals of Statistics*, 49, 1139--1165.
#'   \doi{10.1214/20-AOS1996}
#'
#'   Hallin, M., Nordhausen, K. and Radojicic, U. Directional nonlinear
#'   principal components: a measure transportation approach.
#'
#' @keywords internal
#' @useDynLib ospca, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats optim rnorm predict fitted residuals screeplot
#' @importFrom graphics abline arrows text points barplot lines axis mtext
"_PACKAGE"
