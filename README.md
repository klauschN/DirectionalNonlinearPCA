# **Directional Nonlinear Principal Components: A Measure Transportation Approach**

This project contains code to reproduce almost all of the results of the paper with the same name by Hallin, M. Radojicic, U., and Nordhausen, K.

The goal is to introduce a novel directional nonlinear principal component analysis (NLPCA) based on optimal transport. The code below is once for dimension 2 and once for dimension 3. In addition there 
is some code for linear directional PCA as well as some basic implmentation of GSPCA (Gunsilius and Schennach 2021). Code to reproduce most of the figures of the paper is also included.

**Main functions**

All code here is written in R and requires the packages transport, SphereOptimize and mvtnorm.

The main functions are:

-   dirNLPCA_2D.R Functions for directional NLPCA in the 2D setting
-   dirNLPCA_3D.R Functions for directional NLPCA in the 3D setting
-   dirPCA.R Functions for directional (linear) PCA
-   GSPCA.R Functions for a basic implementation of GSPCA

**Authors**

Hallin, M., Radojicic, U. and Nordhausen, K.

**License**

GNU GPLv3

**References**

Gunsilius, F. and Schennach, S. (2021): Independent Nonlinear Component Analysis. Journal of the American Statistical Association, 118, 1305–1318.\
