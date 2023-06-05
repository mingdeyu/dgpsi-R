# dgpsi 2.2.0

- The efficiency and speed of imputations involved in the training and predictions of DGP emulators are significantly improved (achieving roughly 3x faster training and imputations) by utilizing blocked Gibbs sampling that imputes latent variables layer-wise rather than node-wise. The blocked Gibbs sampling is now the default method for DGP emulator inference and can be changed back to the old node-wise approach by setting `blocked_gibbs = FALSE` in `dgp()`.
- One can now optimize GP components that are contained in the same layer of a DGP emulator in parallel during the DGP emulator training, using multiple cores by setting the new argument `cores` in `dgp()`. This option is useful and can accelerate the training speed when the input dimension is moderately large (in which case there is a large number of GP components to be optimized) and the optimization of GP components is computationally expensive, e.g., when `share = FALSE` in which case input dimensions to individual GP components have different lengthscales.
- Thanks to @tjmckinley, a bug in `update()` when the `object` is an instance of the `dgp` class (that has been trimmed by `window()`) is fixed.
- Thanks to @tjmckinley, some R memory issues due to the underlying Python implementations are rectified.
- `set_seed()` function is added to ensure reproducible results from the package.
- A bug is fixed when candidate sets `x_cand` and `y_cand` are provided to `design()`.
- One can choose different color palettes using the new argument `color` in `plot()` when `style = 2`.
- `set_linked_idx()` allows constructions of different (D)GP emulators (in terms of different connections to the feeding layers) from a same (D)GP emulator. 

# dgpsi 2.1.6

- A bug is found in multi-core predictions in `predict()` when `object` is an instance of `lgp` class and `x` is a list. This bug has been fixed in this version.  
- Thanks to @tjmckinley, an issue (`/usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'GLIBCXX_3.4.30' not found`) encountered in Linux machines is fixed automatically during the execution of `init_py()`.
- `gp()` and `dgp()` allow users to specify the value of scale parameters and whether to estimate the parameters.
- `gp()` and `dgp()` allow users to specify the bounds of lengthscales.
- The jointly robust prior (Gu, 2019) is implemented as the default inference approach for GP emulators in `gp()`.
- The default value of `lengthscale` in `gp()` is changed from `0.2` to `0.1`, and the default value for `nugget` in `gp()` is changed from `1e-6` to `1e-8` if `nugget_est = FALSE`.
- One can now specify the number of GP nodes in each layer (except for the final layer) of a DGP emulator via the `node` argument in `dgp()`.
- Training data are now contained in the S3 classes `gp` and `dgp`.
- The RMSEs (without the min-max normalization) of emulators are now contained in the S3 classes `gp`, `dgp`, and `lgp` after the execution of `validate()`.
- `window()` function is added to trim the traces and obtain new point estimates of DGP model parameters for predictions.
- The min-max normalization can now be switched off in `plot()` by setting the value of `min_max`.
- The default number of imputations `B` for `dgp()` is changed from `50` to `30` to better balance the uncertainty and the speed of DGP emulator predictions. A new function `set_imp()` is made available to change the number of imputations of a trained DGP emulator so one can either achieve faster predictions by further reducing the number of imputations, or account for more imputation uncertainties by increasing the number of imputations, without re-training the emulator.
- The default number of imputations `B` for `continue()` is set to `NULL`, in which case the same number of imputations used in `object` will be applied.
- `nugget` argument of `dgp()` now specifies the nugget values for GP nodes in different layers rather than GP nodes in the final layer.
- The speed of predictions from DGP emulators with squared exponential kernels is significantly improved and is roughly 3x faster than the implementations in version `2.1.5`. 
- The implementation of sequential designs (with two vignettes) of (D)GP emulators using different criterion is made available.
- Thanks to @tjmckinley, an internal reordering issue in `plot()` is fixed.
- `init_py()` now allows users to reinstall and uninstall the underlying Python environment.
- A bug that occurs when a linked DGP emulator involves a DGP emulator with external inputs is fixed.
- `Intel SVML` will now be installed with the Python environment automatically for Intel users for faster implementations.

# dgpsi 2.1.5

- Initial release of the R interface to the Python package `dgpsi v2.1.5`.
