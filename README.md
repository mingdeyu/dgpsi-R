# dgpsi
<!-- badges: start -->
  [![R-CMD-check](https://github.com/mingdeyu/dgpsi_R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mingdeyu/dgpsi-R/actions/workflows/R-CMD-check.yaml)
  ![GitHub](https://img.shields.io/github/license/mingdeyu/DGP)
  [![DOI](https://img.shields.io/badge/DOI-10.1137%2F20M1323771-informational)](https://epubs.siam.org/doi/abs/10.1137/20M1323771)
<!-- badges: end -->
  
The R package `dgpsi` provides R interface to Python package [`dgpsi`](https://github.com/mingdeyu/DGP) for deep and linked Gaussian process emulations. The package currently has following features:

* Deep Gaussian process emulation with flexible architecture construction: 
    - multiple layers;
    - multiple GP nodes;
    - separable or non-separable squared exponential and Mat&eacute;rn-2.5 kernels;
    - global input connections;
    - non-Gaussian likelihoods (Poisson, Negative-Binomial, heteroskedastic Gaussian, and more to come);
* Linked emulation of feed-forward systems of computer models:
    - linking GP emulators of deterministic individual computer models;
    - linking GP and DGP emulators of deterministic individual computer models;
* Multi-core predictions from GP, DGP, and Linked (D)GP emulators.
* Fast Leave-One-Out (LOO) and Out-Of-Sample (OOS) validations for GP, DGP, and linked (D)GP emulators.

## Documentation
See [https://mingdeyu.github.io/dgpsi-R](https://mingdeyu.github.io/dgpsi-R/) to learn more about the package.

## Installation
You can install development version of the package from the GitHub repo:

1. In your RStudio Console, type:
```r
devtools::install_github('mingdeyu/dgpsi-R')
```

2. Restart your RStudio.

3. Load the package and initialize the required Python environment:
```r
library(dgpsi)
init_py()
```

## References
> [Ming, D., Williamson, D., and Guillas, S. (2022) Deep Gaussian process emulation using stochastic imputation. <i>Technometrics</i> (to appear).](https://arxiv.org/abs/2107.01590)

> [Ming, D. and Guillas, S. (2021) Linked Gaussian process emulation for systems of computer models using Mat&eacute;rn kernels and adaptive design, <i>SIAM/ASA Journal on Uncertainty Quantification</i>. 9(4), 1615-1642.](https://epubs.siam.org/doi/abs/10.1137/20M1323771)
