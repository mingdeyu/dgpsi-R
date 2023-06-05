# dgpsi
  [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/dgpsi)](https://CRAN.R-project.org/package=dgpsi)
  [![Download](https://cranlogs.r-pkg.org/badges/grand-total/dgpsi?color=brightgreen)](https://CRAN.R-project.org/package=dgpsi)
  [![R-CMD-check](https://github.com/mingdeyu/dgpsi_R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mingdeyu/dgpsi-R/actions/workflows/R-CMD-check.yaml)
  [![DOC](https://img.shields.io/badge/DOC-release-brightgreen)](https://mingdeyu.github.io/dgpsi-R/)
  [![REF](https://img.shields.io/badge/REF-Linked%20GP-informational)](https://doi.org/10.1137/20M1323771)
  [![REF](https://img.shields.io/badge/REF-Deep%20GP-informational)](https://doi.org/10.1080/00401706.2022.2124311)
  [![REF](https://img.shields.io/badge/REF-Linked%20DGP-informational)](https://arxiv.org/abs/2306.01212)
  [![python](https://img.shields.io/badge/Python-dgpsi%20v2.2.0-informational)](https://github.com/mingdeyu/DGP)
  ![CRAN_License](https://img.shields.io/cran/l/dgpsi?color=green)
  
The R package `dgpsi` provides R interface to Python package [`dgpsi`](https://github.com/mingdeyu/DGP) for deep and linked Gaussian process emulations. 

> **Hassle-free Python Setup**  
> You don't need prior knowledge of Python to start using the package, all you need is a single click in R (see [Installation](#installation) section below) that automatically installs and activates the required Python environment for you!

## Features
`dgpsi` currently has following features:

* Gaussian process emulations with separable or non-separable squared exponential and Mat&eacute;rn-2.5 kernels.
* Deep Gaussian process emulations with flexible structures including: 
    - multiple layers;
    - multiple GP nodes;
    - separable or non-separable squared exponential and Mat&eacute;rn-2.5 kernels;
    - global input connections;
    - non-Gaussian likelihoods (Poisson, Negative-Binomial, and heteroskedastic Gaussian).
* Linked emulations of feed-forward systems of computer models by linking (D)GP emulators of deterministic individual computer models.
* Fast Leave-One-Out (LOO) and Out-Of-Sample (OOS) validations for GP, DGP, and linked (D)GP emulators.
* Multi-core predictions and validations for GP, DGP, and Linked (D)GP emulators.
* Sequential designs for (D)GP emulators and bundles of (D)GP emulators.

## Getting started
* Check [A Quick Guide to dgpsi](https://mingdeyu.github.io/dgpsi-R/articles/dgpsi.html) to get started with the package.
* For experimental features, check out our [website](https://mingdeyu.github.io/dgpsi-R/dev/) for the development version.

## Installation
You can install the package from CRAN:

```r
install.packages('dgpsi')
```

or its development version from GitHub:

```r
devtools::install_github('mingdeyu/dgpsi-R')
```

After the installation, run 

```r
library(dgpsi)
init_py()
```

to install and activate the required Python environment. That's it, the package is now ready to use!

> **Note**  
> Always run `init_py()` after `library(dgpsi)`, telling R to invoke the required Python environment.
> 
> If you experience issues while running `init_py()`, please try to reinstall the Python environment:    
> 
> ```r
> dgpsi::init_py(reinstall = T)
> ```
> 
> or uninstall completely the Python environment:
> 
> ```r
> dgpsi::init_py(uninstall = T)
> ```
> 
> And then restart the R and rerun:
>
> ```r
> library(dgpsi)
> init_py()
> ```

## References
> [Ming, D. and Williamson, D. (2023) Linked deep Gaussian process emulation for model networks. arXiv:2306.01212](https://arxiv.org/abs/2306.01212)

> [Ming, D., Williamson, D., and Guillas, S. (2023) Deep Gaussian process emulation using stochastic imputation. <i>Technometrics</i>. 65(2), 150-161.](https://doi.org/10.1080/00401706.2022.2124311)

> [Ming, D. and Guillas, S. (2021) Linked Gaussian process emulation for systems of computer models using Mat&eacute;rn kernels and adaptive design, <i>SIAM/ASA Journal on Uncertainty Quantification</i>. 9(4), 1615-1642.](https://doi.org/10.1137/20M1323771)
