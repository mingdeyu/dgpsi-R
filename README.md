# dgpsi
  [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/dgpsi)](https://CRAN.R-project.org/package=dgpsi)
  [![Download](https://cranlogs.r-pkg.org/badges/grand-total/dgpsi?color=brightgreen)](https://CRAN.R-project.org/package=dgpsi)
  [![R-CMD-check](https://github.com/mingdeyu/dgpsi_R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mingdeyu/dgpsi-R/actions/workflows/R-CMD-check.yaml)
  ![CRAN_License](https://img.shields.io/cran/l/dgpsi?color=green)
  [![DOC](https://img.shields.io/badge/DOC-Linked%20GP-informational)](https://epubs.siam.org/doi/abs/10.1137/20M1323771)
  [![DOC](https://img.shields.io/badge/DOC-Deep%20GP-informational)](https://doi.org/10.1080/00401706.2022.2124311)
  
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

## Getting started
Check [A Quick Guide to dgpsi](https://mingdeyu.github.io/dgpsi-R/articles/dgpsi.html) to get started with the package.

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

> **Warning**  
> If you are Linux users and encountered importing errors similar to
>
> ```
> /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'GLIBCXX_3.4.30' not found
> ```
>
> during the execution of `init_py()`, please try the following procedure to resolve the issue:
> 
> 1. Open your terminal.
> 2. Activate the conda environment
>
>    ```
>    conda activate dgp_si_R_X_Y_Z
>    ```
>    if you installed the release version `X.Y.Z`, or
>    
>    ```
>    conda activate dgp_si_R_X_Y_Z
>    ```
>    if you installed the development version `X.Y.Z-9000`.
> 3. Install the `libstdcxx-ng` package:
>
>    ```
>    conda install -c conda-forge libstdcxx-ng
>    ```
> 4. Search for the path of the latest version of `libstdc++.so` (e.g., `libstdc++.so.6.0.30`):
>
>    ```
>    find / -name "libstdc++.so*"
>    ```
> 5. Remove the current version of `libstdc++.so` installed in the system:
>
>    ```
>    sudo rm /usr/lib/x86_64-linux-gnu/libstdc++.so.6
>    ```
> 6. Create a symlink in the system lib to the latest version of `libstdc++.so` found in Step 4:
>
>    ```
>    sudo ln -s PATH-FOUND-IN-STEP4 /usr/lib/x86_64-linux-gnu/libstdc++.so.6
>    ```
> 7. Restart R and reload the package.

## References
> [Ming, D., Williamson, D., and Guillas, S. (2022) Deep Gaussian process emulation using stochastic imputation. <i>Technometrics</i>. 0(0), 1-12.](https://doi.org/10.1080/00401706.2022.2124311)

> [Ming, D. and Guillas, S. (2021) Linked Gaussian process emulation for systems of computer models using Mat&eacute;rn kernels and adaptive design, <i>SIAM/ASA Journal on Uncertainty Quantification</i>. 9(4), 1615-1642.](https://epubs.siam.org/doi/abs/10.1137/20M1323771)
