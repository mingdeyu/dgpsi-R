# dgpsi
  [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/dgpsi)](https://CRAN.R-project.org/package=dgpsi)
  [![Download](https://cranlogs.r-pkg.org/badges/grand-total/dgpsi?color=brightgreen)](https://CRAN.R-project.org/package=dgpsi)
  [![R-CMD-check](https://github.com/mingdeyu/dgpsi_R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mingdeyu/dgpsi-R/actions/workflows/R-CMD-check.yaml)
  [![DOC](https://img.shields.io/badge/DOC-release-brightgreen)](https://mingdeyu.github.io/dgpsi-R/)
  [![REF](https://img.shields.io/badge/REF-Linked%20GP-informational)](https://doi.org/10.1137/20M1323771)
  [![REF](https://img.shields.io/badge/REF-Deep%20GP-informational)](https://doi.org/10.1080/00401706.2022.2124311)
  [![REF](https://img.shields.io/badge/REF-Linked%20DGP-informational)](https://arxiv.org/abs/2306.01212)
  [![python](https://img.shields.io/badge/Python-dgpsi%20v2.4.0-informational)](https://github.com/mingdeyu/DGP)
  ![CRAN_License](https://img.shields.io/cran/l/dgpsi?color=green)
  
The R package `dgpsi` provides R interface to Python package [`dgpsi`](https://github.com/mingdeyu/DGP) for deep and linked Gaussian process emulations using stochastic imputation (SI). 

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
* Automatic pruning of DGP emulators, both statically and dynamically.

## Getting started
* Check [A Quick Guide to dgpsi](https://mingdeyu.github.io/dgpsi-R/articles/dgpsi.html) to get started with the package.
* For experimental features, check out our [website for the development version](https://mingdeyu.github.io/dgpsi-R/dev/).

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
```

to load the package. To install or activate the required Python environment automatically, simply run a function from the package. That's it, the package is now ready to use!

> **Note**  
> After loading `dgpsi`, the package may take some time to compile and initiate the underlying Python environment the first
> time a function from `dgpsi` is executed. Any subsequent function calls won't require re-compiling or re-activation of the 
> Python environment, and will be faster.
>
> If you experience Python related issues while using the package, please try to reinstall the Python environment:    
> 
> ```r
> dgpsi::init_py(reinstall = T)
> ```
> 
> Or uninstall completely the Python environment:
> 
> ```r
> dgpsi::init_py(uninstall = T)
> ```
> 
> and then reinstall:
>
> ```r
> dgpsi::init_py()
> ```

## Research Notice
This package is part of an ongoing research initiative. For detailed information about the research aspects and guidelines for use, please refer to our [Research Notice](https://github.com/mingdeyu/dgpsi-R/blob/master/RESEARCH-NOTICE.md).

## References
> [Ming, D. and Williamson, D. (2023) Linked deep Gaussian process emulation for model networks. arXiv:2306.01212](https://arxiv.org/abs/2306.01212)

> [Ming, D., Williamson, D., and Guillas, S. (2023) Deep Gaussian process emulation using stochastic imputation. <i>Technometrics</i>. 65(2), 150-161.](https://doi.org/10.1080/00401706.2022.2124311)

> [Ming, D. and Guillas, S. (2021) Linked Gaussian process emulation for systems of computer models using Mat&eacute;rn kernels and adaptive design, <i>SIAM/ASA Journal on Uncertainty Quantification</i>. 9(4), 1615-1642.](https://doi.org/10.1137/20M1323771)
