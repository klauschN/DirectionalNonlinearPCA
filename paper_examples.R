###############################################################################
## Examples of Hallin, Nordhausen and Radojicic, "Directional nonlinear
## principal components: a measure transportation approach", computed with the
## R package ospca (version 1.0.0).
##
## The script reproduces the examples of the manuscript: the linear directional
## principal components of Section 2.2 (three bivariate examples and a
## trivariate one, Figures 1 to 4), the first illustration of nonlinear
## directional PCA in Section 2.4 (the banana), the six simulated settings of
## Section 5.1, the Old Faithful data of Section 5.2 and the Gvessel data of
## Section 5.3. Only the figure of Section 3.1 (the grid and the empirical
## quantile contours of a mixture of three bivariate Gaussians) is not included,
## as the mixture is not specified in the manuscript. The statistics reported in
## the text (the MSE values) are printed in the console. In
## the examples of Section 5, the nonlinear directional principal curves are
## compared with linear PCA, the principal curve of Hastie and Stuetzle (1989),
## an autoencoder (Scholz et al., 2005) and GS-PCA (Gunsilius and Schennach,
## 2023).
##
## Besides ospca, the bivariate examples of Section 2.2 need the package mvtnorm
## (and MASS, which comes with R), and the comparisons of Section 5 need the
## packages princurve (CRAN) and pcaMethods (Bioconductor, installed with
## BiocManager::install("pcaMethods")). The Gvessel data are contained in the
## package LPCM, and Setting 6 needs the package png together with the image of
## the bird silhouette. Missing packages or files are reported, and the
## corresponding parts are skipped.
##
## Run the script with the working directory set to its folder (in RStudio, for
## example, with Session > Set Working Directory > To Source File Location), as
## the image for Setting 6 is read from there and the figures are written to the
## subfolder paper_figures. The data and the fitted models of all examples are
## saved in paper_figures/results. A second run reuses them and only redraws the
## figures (see recompute below), so that the figures can be modified without
## fitting the models again. The saved results can also be loaded directly, e.g.,
## readRDS("paper_figures/results/section5_2_old_faithful.rds").
##
## Note that the smoothing parameters reported in the manuscript refer to the
## earlier implementation, which applied the Tikhonov regularization twice and
## to the two-sided curve. The values used below were tuned by U. Radojicic for
## the current implementation, which is also the one in ospca.
###############################################################################

library(ospca)

## ===========================================================================
## Settings
## ===========================================================================
out_dir   <- "paper_figures"          # the figures are saved as PDF files in this folder
bird_file <- "bird_silhouette.png"   # image for Setting 6: the silhouette from the source given
                                      # in the manuscript, placed on a white background (a dark
                                      # bird on a light or transparent background also works)
run_ae    <- TRUE                     # fit the autoencoder of pcaMethods
ae_steps  <- NULL                     # NULL uses the default number of training steps of
                                      # pcaMethods::nlpca(), as in the manuscript (on a laptop
                                      # about 2 minutes for each simulated data set); a small value,
                                      # e.g., 20, gives a quick but untrained autoencoder.
nthreads  <- NULL                     # threads of the compiled code (NULL: all cores but one)
make_html <- TRUE                     # also save the 3D figures as interactive web pages (needs
                                      # plotly and htmlwidgets)
pt_cex    <- 0.6                      # size of the plotting symbols of the data in all figures
text_width <- 7                       # text width of the manuscript in inches (geometry package:
                                      # total = {7in, 9in}); the figures are drawn at the size at
                                      # which they are included, see fig_pdf() below
fig_pointsize <- 9                    # font size of the figures in points
html_dir  <- "docs"                   # folder of the web pages, one page per figure and an index
                                      # page; in a GitHub repository, GitHub Pages can publish it
res_dir   <- file.path(out_dir, "results")
                                      # the data and the fitted models are saved in this folder
recompute <- FALSE                    # FALSE: results saved by an earlier run are reused, so that
                                      # a second run only redraws the figures (in seconds); a result
                                      # is recomputed automatically when its code or settings have
                                      # changed. TRUE recomputes everything.

## Simulated settings: n = nR * nS + n0 = 2561 with nR = 20, nS = 128 and n0 = 1.
n_sim  <- 2561
nR_sim <- 20
nS_sim <- 128
eps_sim <- c("1" = 0.5, "2" = 1.5, "3" = 1.5, "4" = 1.5, "5" = 1.5, "6" = 1)
                                      # Tikhonov parameter of the criterion, i.e., for the choice
                                      # of the directions (manuscript: 0.5, 0.5, 0.5, 2, 2 and 1)
eps_curves_sim <- 1                   # Tikhonov parameter of the curves of the fits
eps_display_sim <- c("1" = 1, "2" = 2, "3" = 1, "4" = 1, "5" = 1, "6" = 1)
                                      # Tikhonov parameter of the curves drawn in Figure 7; it only
                                      # affects the drawing, so it can be changed without refitting
delta_sim <- 0.1                      # trimming proportion of the criterion

dir.create(out_dir, showWarnings = FALSE)
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
if (!has_pkg("princurve")) message("Package 'princurve' is not installed: the principal curves are skipped.")
if (run_ae && !has_pkg("pcaMethods"))
  message("Package 'pcaMethods' is not installed (see BiocManager::install(\"pcaMethods\")): the autoencoder is skipped.")

## ===========================================================================
## Helper functions
## ===========================================================================

## Each figure is drawn at the size at which it is included in the manuscript, so
## that its text has about the size of the text of the paper (fig_pointsize). A
## figure drawn with width_frac = w is meant to be included with
##   \includegraphics[width=w\textwidth]{...}
## (the width is given for each figure below).
fig_pdf <- function(file, width_frac, height) {
  pdf(file.path(out_dir, file), width = width_frac * text_width, height = height,
      pointsize = fig_pointsize)
  par(mgp = c(1.8, 0.6, 0), tcl = -0.3)
}

## The data and the fits of each example are saved in res_dir and reused in later
## runs, unless recompute is TRUE (see the settings above). The expression is
## only evaluated when no saved result is used.
## Each result is saved together with the code that computes it and the
## settings and functions it depends on (deps), and it is recomputed
## automatically when one of these has changed since it was saved.
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
cached <- function(name, expr, deps = list()) {
  f <- file.path(res_dir, paste0(name, ".rds"))
  signature <- paste(c(deparse(substitute(expr)), deparse(deps)), collapse = "\n")
  if (!recompute && file.exists(f)) {
    saved <- readRDS(f)
    if (is.list(saved) && identical(saved$signature, signature)) {
      message("  (using the saved results in ", f, ")")
      return(saved$value)
    }
    message("  (the code or the settings have changed since the results were saved: recomputing)")
  }
  value <- expr
  saveRDS(list(value = value, signature = signature), f)
  value
}

## Symmetric inverse square root of a covariance matrix.
inv_sqrt <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% diag(1 / sqrt(e$values)) %*% t(e$vectors)
}

## Multivariate normal random numbers.
rmvn <- function(n, mu, Sigma) {
  sweep(matrix(rnorm(n * length(mu)), n) %*% chol(Sigma), 2, mu, "+")
}

## Resampling with replacement produces exact ties, and the optimal assignment is
## then not unique (the potentials have a zero margin). The ties are therefore
## broken by a small jitter of 1/1000 of the standard deviation of each variable.
break_ties <- function(X) {
  X + matrix(rnorm(length(X)), nrow(X)) %*% diag(1e-3 * apply(X, 2, sd), ncol(X))
}

## Curves of the competing methods: the segment of the first linear principal
## component over the range of its scores, the principal curve of Hastie and
## Stuetzle, the curve of the autoencoder (the reconstructed points, ordered by
## their score) and the first GS-PCA curve. Methods whose package is missing are
## left out. With gs_scale = TRUE, GS-PCA is fitted to the standardized data and
## its curve is transformed back to the original units (see the Gvessel example).
comparison_curves <- function(X, gs_eps_curves = 1, gs_scale = FALSE) {
  out <- list()
  pc <- stats::prcomp(X)
  s <- range(pc$x[, 1])
  out$pca <- rbind(pc$center + s[1] * pc$rotation[, 1], pc$center + s[2] * pc$rotation[, 1])
  if (has_pkg("princurve")) {
    hs <- princurve::principal_curve(X)
    out$hs <- hs$s[hs$ord, , drop = FALSE]
  }
  if (run_ae && has_pkg("pcaMethods")) {
    args <- list(as.matrix(X), nPcs = 1, verbose = FALSE)
    if (!is.null(ae_steps)) args$maxSteps <- ae_steps
    ae <- do.call(pcaMethods::nlpca, args)
    rec <- stats::fitted(ae, as.matrix(X))
    out$ae <- rec[order(pcaMethods::scores(ae)[, 1]), , drop = FALSE]
  }
  ## the GS-PCA curves are smoothed slightly for display; this does not affect the
  ## estimated directions
  mu <- if (gs_scale) colMeans(X) else rep(0, ncol(X))
  sdv <- if (gs_scale) apply(X, 2, sd) else rep(1, ncol(X))
  gs <- gspca(sweep(sweep(X, 2, mu), 2, sdv, "/"), eps_curves = gs_eps_curves,
              nthreads = nthreads)
  out$gs <- sweep(sweep(gs$curves[[1]], 2, sdv, "*"), 2, mu, "+")
  out
}

## Each method has the same color in all figures: directional PCA (linear and
## nonlinear) blue, classical PCA black, GS-PCA red, the principal curve magenta
## and the autoencoder orange. The directions of directional PCA are numbered and
## distinguished by the line types in lty_dir (solid, dashed, dotted, ...).
cols <- c(nlos = "blue", gs = "red", pca = "black", hs = "magenta", ae = "orange")
labs <- c(nlos = "directional PCA", gs = "GS-PCA", pca = "PCA",
          hs = "principal curve", ae = "autoencoder")
lty_dir <- c(1, 2, 3, 4, 5, 6)

## Data together with the first two directional principal curves (solid and
## dashed) and the curves of the competing methods, for the coordinates i and j.
## With eps_display, the directional curves are drawn with this Tikhonov
## parameter instead of the one of the fit.
draw_curves <- function(X, fit, cmp, main, i = 1, j = 2, xlim = NULL, ylim = NULL,
                        xlab = colnames(X)[i], ylab = colnames(X)[j], eps_display = NULL) {
  if (is.null(xlim)) xlim <- range(X[, i])
  if (is.null(ylim)) ylim <- range(X[, j])
  plot(X[, i], X[, j], pch = 16, cex = pt_cex, col = "grey70", xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main)
  for (m in c("pca", "hs", "ae", "gs"))
    if (!is.null(cmp[[m]])) lines(cmp[[m]][, i], cmp[[m]][, j], col = cols[m], lwd = 2)
  for (k in 1:2) {
    C <- if (is.null(eps_display)) fit$curves[[k]]
         else nlos_curve(directions(fit)[, k], fit, eps = eps_display)
    lines(C[, i], C[, j], col = cols["nlos"], lwd = 2.5, lty = lty_dir[k])
  }
}

## The directional principal curves `which` of a bivariate fit, numbered at
## their ends, with the center as a black dot.
draw_dir_curves <- function(fit, which = seq_along(fit$curves)) {
  for (k in which) {
    C <- fit$curves[[k]]
    lines(C, col = cols["nlos"], lwd = 2.5, lty = lty_dir[k])
    text(C[nrow(C), 1], C[nrow(C), 2], labels = k, col = cols["nlos"], pos = 3)
  }
  points(fit$center[1], fit$center[2], pch = 19)
}

## Legend below a figure with several panels, in `nrow` rows.
bottom_legend <- function(nrow = 1) {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("bottom", legend = labs, col = cols, lwd = 2, lty = 1,
         ncol = ceiling(length(labs) / nrow), bty = "n", cex = 0.9)
}

## ===========================================================================
## Section 2.2: linear directional principal components
## ===========================================================================

## The data are generated as in the scripts for the figures of the manuscript
## (with the same random number generators, mvtnorm::rmvnorm() and
## MASS::mvrnorm()), so that the figures can be compared directly.

## Arrows of the linear one-sided PCA (blue: the directions, with lengths
## proportional to the one-sided objectives) and of classical PCA (black: the
## eigenvectors, with lengths proportional to the eigenvalues), as in Figure 1.
draw_arrows_2d <- function(Y, fit, os_scale = 1, pca_scale = 2, main = "") {
  ev <- eigen(stats::cov(Y), symmetric = TRUE)
  A_os <- t(directions(fit)) * (os_scale * fit$objective)
  A_pca <- t(ev$vectors) * (pca_scale * ev$values / 2)
  lims <- range(Y, 1.15 * A_os, 1.15 * A_pca)
  plot(Y, pch = 16, cex = pt_cex, col = "grey80", asp = 1, xlim = lims, ylim = lims,
       xlab = "", ylab = "", main = main)
  arrows(0, 0, A_pca[, 1], A_pca[, 2], col = cols["pca"], lwd = 1.3, length = 0.06)
  text(1.12 * A_pca, labels = seq_len(nrow(A_pca)), col = cols["pca"], cex = 0.9)
  arrows(0, 0, A_os[, 1], A_os[, 2], col = cols["nlos"], lwd = 1.3, length = 0.06)
  text(1.12 * A_os, labels = seq_len(nrow(A_os)), col = cols["nlos"], cex = 0.9)
}

## Reconstructions with k = 1, ..., 2d directions (one row of Figure 3).
draw_reconstructions_2d <- function(Y, fit) {
  for (k in seq_len(2 * ncol(Y))) {
    RC <- reconstruct(fit, n_components = k, add_center = FALSE)
    mse <- reconstruction_mse(Y, RC)
    plot(RC, pch = 16, cex = pt_cex, col = "darkgray", asp = 1, xlim = range(Y[, 1]),
         ylim = range(Y[, 2]), xlab = "", ylab = "",
         main = paste0("k = ", k, ", MSE = ", if (round(mse, 4) == 0) "0" else round(mse, 3)))
  }
}

## Reconstruction errors of classical PCA, for comparison.
pca_mse <- function(Y, ks) {
  V <- eigen(stats::cov(Y), symmetric = TRUE)$vectors
  sapply(ks, function(k) reconstruction_mse(Y, tcrossprod(Y %*% V[, 1:k, drop = FALSE],
                                                          V[, 1:k, drop = FALSE])))
}

if (has_pkg("mvtnorm")) {
  message("Section 2.2: bivariate examples ...")
  lin2d <- cached("section2_2_bivariate", {
    ## (i) an elliptical (Gaussian) distribution
    set.seed(12345)
    Y_ex1 <- mvtnorm::rmvnorm(n = 2000, mean = c(4, 0), sigma = matrix(c(3, 1, 1, 2), 2))
    ## (ii) an asymmetric Gaussian mixture
    set.seed(123)
    Y_ex2 <- rbind(MASS::mvrnorm(1000, c(4, 0), matrix(c(8, 0, 0, 0.1), 2)),
                   MASS::mvrnorm(3000, c(0, 0), matrix(c(0.1, 0, 0, 10), 2)),
                   MASS::mvrnorm(750, c(0, 0), matrix(c(5, 0, 0, 5), 2)))
    ## (iii) a second asymmetric Gaussian mixture, standardized
    set.seed(123)
    Y_ex3 <- rbind(MASS::mvrnorm(2000, c(12, 1), matrix(c(80, 0, 0, 1), 2)),
                   MASS::mvrnorm(2000, c(-4, 1), matrix(c(0.2, 0, 0, 40), 2)),
                   MASS::mvrnorm(750, c(0, 0), matrix(c(40, 0, 0, 40), 2)))
    lapply(list(scale(Y_ex1, scale = FALSE), scale(Y_ex2, scale = FALSE), scale(Y_ex3)),
           function(Y) { Y <- unname(Y); list(Y = Y, fit = ospca(Y, center = FALSE)) })
  })
  ## lengths of the arrows in Figure 1 (one-sided, classical), per example
  arrow_scales <- list(c(1, 2), c(1, 2), c(2, 4))

  ## Figure 1 (width = \textwidth): one-sided and classical principal directions,
  ## the three examples side by side (no panel titles: the examples are referred
  ## to as the left, middle and right panel in the caption)
  fig_pdf("linear_directions_2d.pdf", 1, 2.55)
  par(mfrow = c(1, 3), mar = c(2.2, 2.2, 0.5, 0.5))
  for (i in seq_along(lin2d))
    draw_arrows_2d(lin2d[[i]]$Y, lin2d[[i]]$fit, arrow_scales[[i]][1], arrow_scales[[i]][2])
  dev.off()

  ## Figure 3 (width = \textwidth): reconstructions with k = 1, ..., 4 one-sided
  ## directions, one row per example
  fig_pdf("linear_reconstructions_2d.pdf", 1, 5.6)
  par(mfrow = c(3, 4), mar = c(2.2, 2.2, 1.6, 0.5))
  for (i in seq_along(lin2d)) with(lin2d[[i]], draw_reconstructions_2d(Y, fit))
  dev.off()

  for (i in seq_along(lin2d))
    cat(sprintf("Example (%s): MSE one-sided k = 1, ..., 4: %s; PCA k = 1, 2: %s\n", tolower(as.roman(i)),
                paste(round(sapply(1:4, function(k) reconstruction_mse(lin2d[[i]]$Y,
                  reconstruct(lin2d[[i]]$fit, n_components = k, add_center = FALSE))), 3), collapse = ", "),
                paste(round(pca_mse(lin2d[[i]]$Y, 1:2), 3), collapse = ", ")))
} else {
  message("The bivariate examples of Section 2.2 are skipped: they need the package 'mvtnorm'.")
}

## The trivariate example: two intersecting half-planes together with a ridge
## along the second axis (Figures 2 and 4).
message("Section 2.2: trivariate example ...")
ex_3d <- cached("section2_2_trivariate", {
  set.seed(1234)
  plane1 <- cbind(runif(1000, 0, 20)^2, runif(1000, -100, 100), rnorm(1000, 0, 2)^2)
  plane2 <- cbind(rnorm(1000, 0, 2), runif(1000, -100, 100), runif(1000, 0, 20)^2)
  ridge <- cbind(rnorm(4000, 0, 0.03), runif(4000, -100, 100), rnorm(4000, 0, 0.03))
  Y <- scale(rbind(plane1, plane2, ridge), scale = FALSE)
  dimnames(Y) <- list(NULL, c("X1", "X2", "X3"))
  list(Y = Y, fit = ospca(Y, center = FALSE))
})
Y_3d <- ex_3d$Y
fit_3d <- ex_3d$fit
print(fit_3d)

## Static view of points in three dimensions inside the box lim^3, with arrows
## from the origin; the arrows in A (one-sided directions) are drawn in blue and
## those in B (classical PCA) in black. Points outside the box are not drawn, and
## arrows are shortened, if necessary, so that they end within the fraction
## `clip` of the box.
## The labels of the arrows in A are placed at the side label_pos of the tips
## (one value per arrow, recycled; 1 below, 2 left, 3 above, 4 right).
draw_3d <- function(P, lim, A = NULL, B = NULL, clip = 1, main = "", theta = 45, phi = 20,
                    label_cex = 0.75, label_font = 1, label_pos = 3) {
  pm <- persp(lim, lim, matrix(NA, 2, 2), zlim = lim, theta = theta, phi = phi,
              xlab = "X1", ylab = "X2", zlab = "X3", main = main, border = NA,
              ticktype = "simple", cex.main = 0.9)
  inside <- rowSums(P < lim[1] | P > lim[2]) == 0
  points(trans3d(P[inside, 1], P[inside, 2], P[inside, 3], pm), pch = 16,
         cex = 0.6 * pt_cex, col = "grey60")
  o <- trans3d(0, 0, 0, pm)
  for (M in list(list(A, cols[["nlos"]], label_pos), list(B, cols[["pca"]], 3))) {
    if (is.null(M[[1]])) next
    pos <- rep_len(M[[3]], nrow(M[[1]]))
    for (i in seq_len(nrow(M[[1]]))) {
      a <- M[[1]][i, ]
      s <- ifelse(a > 0, clip * lim[2] / a, ifelse(a < 0, clip * lim[1] / a, Inf))
      a <- a * min(1, s)
      e <- trans3d(a[1], a[2], a[3], pm)
      arrows(o$x, o$y, e$x, e$y, col = M[[2]], lwd = 1.5, length = 0.07)
      text(e$x, e$y, labels = i, col = M[[2]], cex = label_cex, font = label_font, pos = pos[i])
    }
  }
}

ev_3d <- eigen(stats::cov(Y_3d), symmetric = TRUE)
A_3d <- t(directions(fit_3d)) * (fit_3d$objective / 6)     # one-sided directions
B_3d <- t(ev_3d$vectors) * (ev_3d$values / 2 / 12)         # classical PCA

## Figure 2 (width = 0.7\textwidth): the data with the one-sided (blue) and
## classical (black) directions; for the manuscript, a snapshot of the interactive
## version (see the end of the script) may be preferable
fig_pdf("linear_directions_3d.pdf", 0.7, 4.6)
par(mar = c(0.5, 0.5, 0.5, 0.5))
draw_3d(Y_3d, range(Y_3d, 1.1 * A_3d, 1.1 * B_3d), A = A_3d, B = B_3d)
dev.off()

## Figure 4 (width = 0.9\textwidth): reconstructions with k = 1, ..., 6
## one-sided directions
lim_3d <- c(min(Y_3d), max(Y_3d)) / 2
fig_pdf("linear_reconstructions_3d.pdf", 0.9, 4.3)
par(mfrow = c(2, 3), mar = c(0.5, 0.3, 1.6, 0.3))
for (k in 1:6) {
  RC <- reconstruct(fit_3d, n_components = k, add_center = FALSE)
  mse <- reconstruction_mse(Y_3d, RC)
  draw_3d(RC, lim_3d, A = A_3d[seq_len(k), , drop = FALSE], clip = 0.8,
          main = paste0("k = ", k, ", MSE = ", if (round(mse, 4) == 0) "0" else round(mse, 3)),
          label_cex = 1, label_font = 2,       # bold numbers of the directions
          label_pos = c(3, 3, 3, 3, 3, 1))     # the label of direction 6 below its tip
}
dev.off()
cat("Trivariate example: MSE one-sided k = 1, ..., 6:",
    round(sapply(1:6, function(k) reconstruction_mse(Y_3d,
      reconstruct(fit_3d, n_components = k, add_center = FALSE))), 3),
    "; PCA k = 1, 2, 3:", round(pca_mse(Y_3d, 1:3), 3), "\n")

## ===========================================================================
## Section 2.4: a first illustration of nonlinear directional PCA
## ===========================================================================

## The banana-shaped distribution of Setting 1 of Section 5.1. The manuscript
## does not state the sample size; we use n = 2561 = 20 * 128 + 1 and the same
## smoothing as in Section 5.1.
message("Section 2.4: banana ...")
banana <- cached("section2_4_banana", {
  set.seed(1)
  x_b <- runif(n_sim, -4, 4)
  X <- cbind(X = x_b, Y = (0.4 * x_b)^2 + rnorm(n_sim, sd = 0.5))
  list(X = X, fit = nlospca(X, nR = nR_sim, nS = nS_sim, eps = eps_sim[["1"]],
                            eps_curves = eps_curves_sim, delta = delta_sim, nthreads = nthreads))
}, deps = list(n_sim, nR_sim, nS_sim, eps_sim[["1"]], eps_curves_sim, delta_sim))
X_b <- banana$X
fit_b <- banana$fit
print(fit_b)
mse_b <- sapply(1:4, function(k) reconstruction_mse(X_b, reconstruct(fit_b, n_components = k)))

## Projections of the bivariate data X onto the manifold spanned by the first k
## directional principal curves (medium grey), together with the data (light
## grey) and the curves.
## With newdata = TRUE, X are not the observations of the fit, and their
## projections are computed as for new observations.
draw_reduction <- function(X, fit, k, main, newdata = FALSE, ...) {
  P <- if (newdata) reconstruct(fit, newdata = X, n_components = k)
       else reconstruct(fit, n_components = k)
  plot(X, pch = 16, cex = pt_cex, col = "grey85", main = main, ...)
  points(P, pch = 16, cex = pt_cex, col = "grey55")
  for (j in seq_len(k))
    lines(fit$curves[[j]], col = cols["nlos"], lwd = 2.5, lty = lty_dir[j])
}

## Figure 5, top (width = 0.4\textwidth): the data with the four directional
## principal curves
fig_pdf("nonlinear_banana_curves.pdf", 0.4, 2.5)
par(mar = c(2.8, 2.8, 0.5, 0.5))
plot(X_b, pch = 16, cex = pt_cex, col = "grey70")
draw_dir_curves(fit_b)
dev.off()

## Figure 5, bottom (width = 0.65\textwidth): the projections onto the manifolds
## spanned by the first k principal curves
fig_pdf("nonlinear_banana_reduction.pdf", 0.65, 3.9)
par(mfrow = c(2, 2), mar = c(2.8, 2.8, 1.6, 0.5))
for (k in 1:4)
  draw_reduction(X_b, fit_b, k, main = paste0("k = ", k, ", MSE = ", round(mse_b[k], 3)))
dev.off()
cat("Section 2.4 (banana): MSE k = 1, ..., 4:", round(mse_b, 3), "\n")

## ===========================================================================
## Section 5.1: simulated examples
## ===========================================================================

## Data-generating mechanisms of the six settings. Settings 1, 3 and 6 are as
## described in the manuscript, Settings 2, 4 and 5 as in the code of the
## published figure (and in U. Radojicic's scripts), for which the smoothing
## parameters eps_sim were tuned; the text of the manuscript should describe
## these. Note that N(0, 0.25) denotes the variance.
sim_data <- function(setting, n) {
  X <- switch(setting,
    ## Setting 1: banana-shaped distribution
    "1" = {
      x1 <- runif(n, -4, 4)
      cbind(x1, (0.4 * x1)^2 + rnorm(n, sd = 0.5))
    },
    ## Setting 2: noisy banana, as in the code of the figure (and in U. Radojicic's
    ## scripts): 2001 points of a banana with X1 ~ U(-6, 6) and 560 points of
    ## uniform noise on [-7, 7] times the range of the banana, standardized and
    ## whitened (multiplied by the inverse square root of the covariance matrix)
    "2" = {
      stopifnot(n == 2561)
      set.seed(1)
      x1 <- runif(2001, -6, 6)
      x2 <- (0.4 * x1)^2 + 0.5 * rnorm(2001)
      noise <- cbind(runif(560, -7, 7), runif(560, round(min(x2)), round(max(x2))))
      W <- scale(rbind(cbind(x1, x2), noise))
      W %*% inv_sqrt(stats::cov(W))
    },
    ## Setting 3: exponential transformation of a Gaussian vector
    "3" = {
      x1 <- rnorm(n); x2 <- rnorm(n)
      y2 <- exp(x1 + x2)
      cbind(y2 + x1, y2)
    },
    ## Setting 4: sigmoid-shaped distribution, as in the code of the figure: the
    ## logit of U(0, 1) variables plus 0.4 times t_5 noise, centered
    "4" = {
      x1 <- runif(n, 0, 1)
      x2 <- log(abs(x1 / (x1 - 1))) + 0.4 * rt(n, 5)
      sweep(cbind(x1, x2), 2, colMeans(cbind(x1, x2)))
    },
    ## Setting 5: two-component Gaussian mixture, as in the code of the figure:
    ## 1200 observations from N((3, 3), Sigma_1), squared componentwise, and 1361
    ## from N((1, -1), Sigma_2) multiplied by 5, standardized
    "5" = {
      rnorm2 <- if (has_pkg("mvtnorm")) mvtnorm::rmvnorm else rmvn
      Z1 <- rnorm2(1200, c(3, 3), matrix(c(1.5, 1, 1, 1.5), 2))^2
      Z2 <- 5 * rnorm2(n - 1200, c(1, -1), matrix(2 * c(0.2, -0.2, -0.2, 0.9), 2))
      scale(rbind(Z1, Z2))
    },
    ## Setting 6: noisy bird-shaped distribution. The foreground pixels of the
    ## thresholded image are sampled uniformly, and each sampled point is placed
    ## uniformly within its pixel, which avoids ties. Note that the original image
    ## of the source has a transparent background (all pixels are black, and only
    ## the alpha channel separates the bird from the background). The file
    ## bird_silhouette.png is this image placed on a white background; an image
    ## with an alpha channel is placed on a white background here as well.
    "6" = {
      img <- png::readPNG(bird_file)
      g <- if (length(dim(img)) == 3) apply(img[, , 1:3, drop = FALSE], c(1, 2), mean) else img
      if (length(dim(img)) == 3 && dim(img)[3] == 4) g <- img[, , 4] * g + (1 - img[, , 4])
      fg <- which(g < 0.5, arr.ind = TRUE)
      b <- rbinom(n, 1, 2000 / 2561)
      P <- fg[sample(nrow(fg), n, replace = TRUE), , drop = FALSE]
      bird <- cbind((P[, 2] - 1 + runif(n)) / ncol(g), 1 - (P[, 1] - 1 + runif(n)) / nrow(g))
      b * bird + (1 - b) * cbind(runif(n), runif(n))
    })
  X <- unname(as.matrix(X))
  colnames(X) <- c("X1", "X2")
  X
}

setting_names <- c("1" = "Setting 1: banana", "2" = "Setting 2: noisy banana",
                   "3" = "Setting 3: exponential transformation", "4" = "Setting 4: sigmoid",
                   "5" = "Setting 5: Gaussian mixture", "6" = "Setting 6: noisy bird")
## plot windows for the heavy-tailed Settings 3 and 4, as in the published figure
## (the other settings use the range of the data)
windows <- list("3" = list(xlim = c(-2, 4), ylim = c(-0.5, 5)), "4" = list(xlim = c(-0.6, 0.6)))

settings <- c("1", "2", "3", "4", "5", "6")
if (!(file.exists(bird_file) && has_pkg("png"))) {
  message("Setting 6 is skipped: it needs the package 'png' and the file '", bird_file, "'.")
  settings <- setdiff(settings, "6")
}

sim <- list()
for (s in settings) {
  message("Simulated ", setting_names[s], " ...")
  sim[[s]] <- cached(paste0("section5_1_setting", s), {
    set.seed(123)
    X <- sim_data(s, n_sim)
    fit <- nlospca(X, nR = nR_sim, nS = nS_sim, eps = eps_sim[[s]],
                   eps_curves = eps_curves_sim, delta = delta_sim, nthreads = nthreads)
    list(X = X, fit = fit, cmp = comparison_curves(X))
  }, deps = list(n_sim, nR_sim, nS_sim, eps_sim[[s]], eps_curves_sim, delta_sim,
                sim_data, comparison_curves, run_ae, ae_steps,
                if (s == "6") bird_file))
}

## Figure 7 (width = \textwidth): the six settings
fig_pdf("simulated_settings.pdf", 1, 8.2)
par(mfrow = c(3, 2), oma = c(2, 0, 0, 0), mar = c(2.8, 2.8, 1.6, 0.5))
for (s in names(sim))
  with(sim[[s]], draw_curves(X, fit, cmp, main = setting_names[s],
                             xlim = windows[[s]]$xlim, ylim = windows[[s]]$ylim,
                             eps_display = eps_display_sim[[s]]))
bottom_legend()
dev.off()

## ===========================================================================
## Section 5.2: Old Faithful
## ===========================================================================

## For the directional method, the 272 observations of datasets::faithful are
## resampled with replacement to n = 2561 = 20 * 128 + 1, and the ties are broken
## by a small jitter. Both variables are measured in minutes and are used
## without standardization. The competing methods are fitted to the original
## observations, since they do not need a sample size of this form, and the
## ties created by the resampling make the discrete Jacobians of GS-PCA
## degenerate (for most grid points with the resampled data).
message("Old Faithful ...")
old_faithful <- cached("section5_2_old_faithful", {
  X0 <- as.matrix(faithful[, c("eruptions", "waiting")])
  set.seed(1)
  X <- break_ties(X0[sample(nrow(X0), 2561, replace = TRUE), ])
  rownames(X) <- NULL
  ## The criterion evaluates the curves at mR = 60 radii, three times the radii of
  ## the target grid. With the default mR = nR = 20, the second direction was
  ## orthogonal to the first for all six resamples we tried, so that the first two
  ## curves both lead into the cluster of short eruptions. With mR = 60 and eps = 2
  ## (the value in the manuscript), the second direction was the antipode of the
  ## first for all six resamples, and the first two curves lead into the two
  ## clusters, as in the figure of the manuscript (eps = 4 gave this for five of
  ## the six).
  list(X0 = X0, X = X,
       fit = nlospca(X, nR = 20, nS = 128, mR = 60, eps = 2, delta = 0, nthreads = nthreads),
       cmp = comparison_curves(X0))
}, deps = list(break_ties, comparison_curves, run_ae, ae_steps))
X_of <- old_faithful$X
fit_of <- old_faithful$fit
cmp_of <- old_faithful$cmp
print(fit_of)

## The figures show the 272 original observations (the resampled data contain
## mostly copies of them), and the reconstruction errors are computed for these.
X_of0 <- old_faithful$X0
lab_of <- c("eruption duration (min)", "waiting time (min)")

## Figure 8 (width = 0.7\textwidth): comparison of the methods (with the legend
## in two rows)
fig_pdf("old_faithful_comparison.pdf", 0.7, 4.6)
par(oma = c(3, 0, 0, 0), mar = c(2.8, 2.8, 0.5, 0.5))
draw_curves(X_of0, fit_of, cmp_of, main = "", xlab = lab_of[1], ylab = lab_of[2])
bottom_legend(nrow = 2)
dev.off()

## Figure 9, left (width = 0.48\textwidth): the four directional principal
## curves (plot(fit_of) would use equal scales on both axes, which is not
## suitable here, as the two variables have quite different ranges)
fig_pdf("old_faithful_curves.pdf", 0.48, 3.3)
par(mar = c(2.8, 2.8, 0.5, 0.5))
## room above the highest curve for its label
ylim_of <- range(X_of0[, 2], sapply(fit_of$curves, function(C) C[, 2]))
plot(X_of0, pch = 16, cex = pt_cex, col = "grey70", xlab = lab_of[1], ylab = lab_of[2],
     ylim = ylim_of + c(0, 0.05) * diff(ylim_of))
draw_dir_curves(fit_of)
dev.off()

## Figure 9 (2): the projections onto the manifolds spanned by the first k
## principal curves, k = 1, ..., 4
mse_of <- sapply(1:4, function(k)
  reconstruction_mse(X_of0, reconstruct(fit_of, newdata = X_of0, n_components = k)))
## Figure 9, right (width = 0.48\textwidth); the panels are small, so the axis
## labels are left out (they are those of the left panel)
fig_pdf("old_faithful_reduction.pdf", 0.48, 3.3)
par(mfrow = c(2, 2), mar = c(1.8, 1.8, 1.4, 0.3))
for (k in 1:4)
  draw_reduction(X_of0, fit_of, k, newdata = TRUE, xlab = "", ylab = "",
                 main = paste0("k = ", k, ", MSE = ", round(mse_of[k], 2)))
dev.off()
cat("Old Faithful, reconstruction errors (original observations) for k = 1, ..., 4:",
    round(mse_of, 3), "\n")

## ===========================================================================
## Section 5.3: Gvessel
## ===========================================================================

## Salinity, temperature and depth of LPCM::gvessel, used in their original
## units. As for Old Faithful, the 643 observations are resampled with
## replacement to n = 4001 = 40 * 100 + 1 for the directional method only, and
## the competing methods are fitted to the original observations. GS-PCA,
## however, is fitted to the standardized data: in the original units the
## standard deviations range from 0.29 (salinity) to 1319 (depth), and the
## discrete Jacobians are then not positive definite at almost all grid points.
## Even after standardization about one fifth of the grid points are left out,
## as the data lie close to a curve, and gspca() gives a warning. The warning can
## be ignored here: with other grid sizes or difference steps, the first GS-PCA
## direction changes by at most 3 degrees. Note also that the manuscript uses a
## 20 x 20 x 10 grid for GS-PCA, whereas gspca() uses the same number of grid
## points on each axis (here its default grid).
if (has_pkg("LPCM")) {
  message("Gvessel ...")
  gv <- cached("section5_3_gvessel", {
    utils::data("gvessel", package = "LPCM", envir = environment())
    X0 <- as.matrix(gvessel[, c("salg", "tempg", "depthg")])
    colnames(X0) <- c("salinity", "temperature", "depth")
    rownames(X0) <- NULL
    set.seed(1)
    X <- break_ties(X0[sample(nrow(X0), 4001, replace = TRUE), ])
    list(X0 = X0, X = X,
         fit = nlospca(X, nR = 40, nS = 100, eps = 2, delta = 0, nthreads = nthreads),
                                      # manuscript: alpha = 1 in the earlier implementation
         cmp = comparison_curves(X0, gs_scale = TRUE))
  }, deps = list(break_ties, comparison_curves, run_ae, ae_steps))
  X_gv <- gv$X
  fit_gv <- gv$fit
  cmp_gv <- gv$cmp
  print(fit_gv)

  ## As for Old Faithful, the figures show the 643 original observations.
  X_gv0 <- gv$X0
  pairs_3 <- list(c(1, 2), c(1, 3), c(2, 3))
  ## Figure 10 (width = \textwidth): the three pairs of variables
  fig_pdf("gvessel_comparison.pdf", 1, 2.9)
  par(mfrow = c(1, 3), oma = c(2, 0, 0, 0), mar = c(2.8, 2.8, 0.5, 0.5))
  for (p in pairs_3)
    draw_curves(X_gv0, fit_gv, cmp_gv, main = "", i = p[1], j = p[2])
  bottom_legend()
  dev.off()

  ## projections onto the manifold spanned by the first two principal curves
  ## (medium grey), with the data (light grey) and the two curves
  P2_gv <- reconstruct(fit_gv, newdata = X_gv0, n_components = 2)
  ## (width = \textwidth)
  fig_pdf("gvessel_reduction.pdf", 1, 2.5)
  par(mfrow = c(1, 3), mar = c(2.8, 2.8, 0.5, 0.5))
  for (p in pairs_3) {
    plot(X_gv0[, p], pch = 16, cex = pt_cex, col = "grey85")
    points(P2_gv[, p], pch = 16, cex = pt_cex, col = "grey55")
    for (j in 1:2)
      lines(fit_gv$curves[[j]][, p], col = cols["nlos"], lwd = 2, lty = lty_dir[j])
  }
  dev.off()
} else {
  message("The Gvessel example is skipped: it needs the package 'LPCM' for the data.")
}

## ===========================================================================
## Interactive 3D figures (HTML)
## ===========================================================================

## The three-dimensional figures are also saved as interactive plotly graphics,
## which can be rotated and zoomed in a web browser: one web page per figure in
## the folder html_dir, together with an index page and the shared JavaScript
## libraries (subfolder lib). The whole folder is needed to view the pages. In
## a GitHub repository, GitHub Pages can publish it (Settings > Pages, deploy
## from the branch, folder /docs), and each page then has its own address, e.g.,
##   https://<user>.github.io/<repository>/linear_directions_3d.html
## which can be given in the caption of the static figure. A page can also be
## included in another web page, e.g., with
##   <iframe src="linear_directions_3d.html" width="800" height="650"></iframe>

## Arrows from the origin to the rows of A, with heads of four short barbs (so
## that they look like arrow heads from every angle) of length `head`, and labels.
## The shafts and heads form one trace, and the labels a second one.
cross3 <- function(a, b) c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
                           a[1] * b[2] - a[2] * b[1])
add_arrows_3d <- function(p, A, color, name, head, visible = TRUE, showlegend = TRUE) {
  seg <- NULL
  for (i in seq_len(nrow(A))) {
    a <- A[i, ]
    u <- a / sqrt(sum(a^2))
    v1 <- cross3(u, if (abs(u[3]) < 0.9) c(0, 0, 1) else c(0, 1, 0))
    v1 <- v1 / sqrt(sum(v1^2))
    v2 <- cross3(u, v1)
    base <- a - head * u
    for (b in list(c(0, 0, 0), base + 0.5 * head * v1, base - 0.5 * head * v1,
                   base + 0.5 * head * v2, base - 0.5 * head * v2))
      seg <- rbind(seg, a, b, NA)
  }
  p <- plotly::add_trace(p, type = "scatter3d", mode = "lines", x = seg[, 1], y = seg[, 2],
                         z = seg[, 3], line = list(color = color, width = 6), name = name,
                         showlegend = showlegend, hoverinfo = "name", visible = visible)
  plotly::add_trace(p, type = "scatter3d", mode = "text", x = 1.08 * A[, 1],
                    y = 1.08 * A[, 2], z = 1.08 * A[, 3], text = seq_len(nrow(A)),
                    textfont = list(color = color, size = 14), showlegend = FALSE,
                    hoverinfo = "skip", visible = visible)
}
## Equal scales on the three axes within the box lim^3. The camera gives the
## initial view, which can also be used for the static figure in the paper.
scene_3d <- function(lim, titles = c("X1", "X2", "X3"))
  list(xaxis = list(title = titles[1], range = lim), yaxis = list(title = titles[2], range = lim),
       zaxis = list(title = titles[3], range = lim), aspectmode = "cube",
       camera = list(eye = list(x = 1.5, y = 1.5, z = 0.9)))
## One web page per figure. The pages share the JavaScript libraries in the
## subfolder lib of html_dir, so that each page contains only its own data.
save_html <- function(p, file, title)
  htmlwidgets::saveWidget(p, file.path(html_dir, file), selfcontained = FALSE,
                          libdir = "lib", title = title)

if (make_html && has_pkg("plotly") && has_pkg("htmlwidgets")) {
  message("Interactive 3D figures ...")
  dir.create(html_dir, showWarnings = FALSE)
  html_dir <- normalizePath(html_dir)
  pages <- character(0)                 # titles of the pages written, named by file
  lim_html <- range(Y_3d, 1.15 * A_3d, 1.15 * B_3d)
  head_html <- 0.04 * diff(lim_html)      # length of the arrow heads

  ## Figure 2: the data with the one-sided (blue) and classical (black) directions
  p2 <- plotly::plot_ly()
  p2 <- plotly::add_trace(p2, type = "scatter3d", mode = "markers", x = Y_3d[, 1],
                          y = Y_3d[, 2], z = Y_3d[, 3], name = "data", hoverinfo = "skip",
                          marker = list(size = 2.5 * pt_cex, color = "grey"))
  p2 <- add_arrows_3d(p2, A_3d, cols[["nlos"]], "one-sided PCA", head_html)
  p2 <- add_arrows_3d(p2, B_3d, cols[["pca"]], "classical PCA", head_html)
  p2 <- plotly::layout(p2, title = "Directional and classical PCA", scene = scene_3d(lim_html),
                       margin = list(l = 0, r = 0, b = 0, t = 40))
  save_html(p2, "linear_directions_3d.html", "Directional and classical PCA")
  pages["linear_directions_3d.html"] <- paste("Trivariate example: the data with the",
    "one-sided (blue) and classical (black) principal directions")

  ## Figure 4: the reconstructions with k = 1, ..., 6 directions; the buttons
  ## switch between the values of k
  p4 <- plotly::plot_ly()
  group <- integer(0)
  mse_html <- numeric(6)
  for (k in 1:6) {
    RC <- reconstruct(fit_3d, n_components = k, add_center = FALSE)
    mse_html[k] <- reconstruction_mse(Y_3d, RC)
    p4 <- plotly::add_trace(p4, type = "scatter3d", mode = "markers", x = RC[, 1],
                            y = RC[, 2], z = RC[, 3], hoverinfo = "skip", visible = k == 1,
                            marker = list(size = 2.5 * pt_cex, color = "grey"), showlegend = FALSE)
    p4 <- add_arrows_3d(p4, A_3d[seq_len(k), , drop = FALSE], cols[["nlos"]], "one-sided PCA",
                        head_html, visible = k == 1, showlegend = FALSE)
    group <- c(group, rep(k, 3))            # points, arrows, labels
  }
  title_k <- function(k) paste0("Reconstruction with k = ", k, " directions, MSE = ",
                                format(round(mse_html[k], 2), big.mark = ","))
  buttons <- lapply(1:6, function(k)
    list(method = "update", label = paste0("k = ", k),
         args = list(list(visible = as.list(group == k)), list(title = title_k(k)))))
  p4 <- plotly::layout(p4, title = title_k(1), scene = scene_3d(lim_html),
                       margin = list(l = 0, r = 0, b = 0, t = 40),
                       updatemenus = list(list(type = "buttons", direction = "right",
                                               x = 0, y = 0, xanchor = "left",
                                               yanchor = "bottom", buttons = buttons)))
  save_html(p4, "linear_reconstructions_3d.html", "Directional reconstructions")
  pages["linear_reconstructions_3d.html"] <- paste("Trivariate example: directional dimension",
    "reduction with k = 1, ..., 6 directions")

  ## Gvessel: the data with the first two directional principal curves (solid
  ## and dashed) and the curves of the competing methods
  if (exists("fit_gv")) {
    X0 <- gv$X0
    pg <- plotly::plot_ly()
    pg <- plotly::add_trace(pg, type = "scatter3d", mode = "markers", x = X0[, 1],
                            y = X0[, 2], z = X0[, 3], name = "data", hoverinfo = "skip",
                            marker = list(size = 3 * pt_cex, color = "grey"))
    for (m in c("pca", "hs", "ae", "gs"))
      if (!is.null(cmp_gv[[m]]))
        pg <- plotly::add_trace(pg, type = "scatter3d", mode = "lines", x = cmp_gv[[m]][, 1],
                                y = cmp_gv[[m]][, 2], z = cmp_gv[[m]][, 3], name = labs[[m]],
                                line = list(color = cols[[m]], width = 5))
    for (j in 1:2)
      pg <- plotly::add_trace(pg, type = "scatter3d", mode = "lines",
                              x = fit_gv$curves[[j]][, 1], y = fit_gv$curves[[j]][, 2],
                              z = fit_gv$curves[[j]][, 3],
                              name = paste("directional PCA, curve", j),
                              line = list(color = cols[["nlos"]], width = 7,
                                          dash = if (j == 1) "solid" else "dash"))
    pg <- plotly::layout(pg, title = "Gvessel", margin = list(l = 0, r = 0, b = 0, t = 40),
                         scene = list(xaxis = list(title = colnames(X0)[1]),
                                      yaxis = list(title = colnames(X0)[2]),
                                      zaxis = list(title = colnames(X0)[3]),
                                      aspectmode = "cube"))
    save_html(pg, "gvessel_3d.html", "Gvessel")
    pages["gvessel_3d.html"] <- paste("Gvessel data: the first two directional principal",
      "curves, compared with PCA, the principal curve, the autoencoder and GS-PCA")
  }

  ## index page with links to the figures; the empty file .nojekyll tells GitHub
  ## Pages to publish the folder as it is
  writeLines(c(
    "<!DOCTYPE html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Directional nonlinear principal components: interactive figures</title>",
    "<style>body { font-family: Georgia, serif; max-width: 46em; margin: 2em auto;",
    "padding: 0 1em; line-height: 1.5; color: #222; } a { color: #1a4f8a; }</style>",
    "</head>", "<body>",
    "<h1>Directional nonlinear principal components: interactive figures</h1>",
    paste("<p>Interactive versions of the three-dimensional figures of M. Hallin,",
          "K. Nordhausen and U. Radojicic, <em>Directional nonlinear principal components:",
          "a measure transportation approach</em>. The figures can be rotated with the",
          "mouse and zoomed with the mouse wheel. They were computed with the R package",
          "ospca and the script paper_examples.R.</p>"),
    "<ul>",
    sprintf("<li><a href=\"%s\">%s</a></li>", names(pages), pages),
    "</ul>", "</body>", "</html>"), file.path(html_dir, "index.html"))
  file.create(file.path(html_dir, ".nojekyll"))
} else if (make_html) {
  message("The interactive figures are skipped: they need the packages 'plotly' and 'htmlwidgets'.")
}

message("Done. The figures are in the folder '", out_dir, "'.")
