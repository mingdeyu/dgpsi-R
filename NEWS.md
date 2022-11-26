# dgpsi 2.1.5-9000 (development version)

- A bug is found in multi-core predictions in `predict()` when `object` is an instance of `lgp` class and `x` is a list. This bug has been fixed in the development version `2.1.5-9000` and will be fixed in the next `dgpsi` release version `2.1.6`. The single-core predictions in the current release version `2.1.5` are not affected by the bug.  
- An issue (`/usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'GLIBCXX_3.4.30' not found`) encountered in Linux machines is fixed automatically during the execution of `init_py()`.
- `gp()` and `dgp()` allow users to specify the value of scale parameters and whether to estimate the parameters.
- Training data are now contained in the S3 classes `gp` and `dgp`.
- The RMSEs (without the min-max normalization) of emulators are now contained in the S3 classes `gp`, `dgp`, and `lgp` after the execution of `validate()`.
- The min-max normalization can now be switched off in `plot()` by setting the value of `min_max`.
- The default number of imputations `B` for `dgp()` is changed from `50` to `30` to better balance the uncertainty and the speed of DGP emulator predictions. A new function `set_imp()` is made available to change the number of imputations of a trained DGP emulator so one can either achieve faster predictions by further reducing the number of imputations, or account for more imputation uncertainties by increasing the number of imputations, without re-training the emulator.
- The speed of predictions from DGP emulators with squared exponential kernels is significantly improved and is roughly 3x faster than the implementations in version `2.1.5`. 
- The implementation of sequential designs (with a vignette) is made available.

# dgpsi 2.1.5

- Initial release of the R interface to the Python package `dgpsi v2.1.5`.
