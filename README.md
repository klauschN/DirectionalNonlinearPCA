# **Directional Nonlinear Principal Components: A Measure Transportation Approach**

This project contains code to reproduce almost all of the results of the paper with the same name by Hallin, M., Nordhausen, K. and Radojicic, U.

The goal is to introduce a novel directional nonlinear principal component analysis (NLPCA) based on measure transportation. All methods of the paper are implemented in the R package ospca, which is contained in this repository. In addition, the repository contains a script that reproduces most of the figures of the paper and interactive versions of its three-dimensional figures.

**The R package ospca**

The folder ospca contains the source of the R package ospca. Its main functions are:

-   ospca() Linear directional (one-sided) PCA
-   nlospca() Directional nonlinear PCA, based on the empirical center-outward quantile function
-   gspca() A basic implementation of the independent nonlinear component analysis (GS-PCA) of Gunsilius and Schennach (2023), for comparison

For all three methods, the functions directions(), scores(), predict(), reconstruct(), reconstruction_mse(), screeplot() and plot() give the directions, the scores, the reconstructions and graphical summaries. The package has three vignettes, which illustrate the methods on simple data. Its authors are Nordhausen, K., del Barrio Tellado, E., Hallin, M. and Radojicic, U.

The package contains C++ code, so that installing it from source needs a C++ compiler (Rtools on Windows). It can be installed directly from this repository with

```r
remotes::install_github("klauschN/DirectionalNonlinearPCA", subdir = "ospca", build_vignettes = TRUE)
```

or, from a local copy of the repository, with R CMD INSTALL ospca.

**Reproducing the figures**

The script paper_examples.R reproduces most of the figures and the reported reconstruction errors of the paper with the package ospca: the linear directional principal components of Section 2.2, the first illustration of nonlinear directional PCA in Section 2.4, the simulated settings of Section 5.1 and the real data examples of Section 5.2 (Old Faithful, the Hertzsprung-Russell diagram of the Hipparcos stars, and Gvessel). Besides ospca, it needs the packages mvtnorm, princurve, pcaMethods (from Bioconductor), LPCM, astrodatR and png, and, for the interactive figures, plotly and htmlwidgets. The package astrodatR, which contains the Hipparcos data, is archived on CRAN and can be installed with

```r
install.packages("https://cran.r-project.org/src/contrib/Archive/astrodatR/astrodatR_0.1.tar.gz",
                 repos = NULL, type = "source")
```

The image bird_silhouette.png is needed for the bird of Setting 6; it is the silhouette from the source given in the paper, placed on a white background.

The script should be run with the working directory set to its folder. The figures are saved in the folder paper_figures, and the data and the fitted models in paper_figures/results, so that a second run only redraws the figures. The comments in the script give more details.

**Interactive figures**

The folder docs contains interactive versions of the three-dimensional figures, which can be rotated and zoomed in a web browser. They are published at https://klauschn.github.io/DirectionalNonlinearPCA/.

**Authors**

Hallin, M., Nordhausen, K. and Radojicic, U.

**License**

GNU GPLv3

**References**

Gunsilius, F. and Schennach, S. (2023): Independent Nonlinear Component Analysis. Journal of the American Statistical Association, 118(542), 1305–1318. https://doi.org/10.1080/01621459.2021.1990768

Hallin, M., del Barrio, E., Cuesta-Albertos, J. and Matrán, C. (2021): Distribution and Quantile Functions, Ranks and Signs in Dimension d: A Measure Transportation Approach. Annals of Statistics, 49, 1139–1165. https://doi.org/10.1214/20-AOS1996
