# Changelog

## dgpsi 2.6.0-9000 (development version)

- Resolved Python environment import failure relating to `libcblas.so.3`
  on Intel CPUs when using MKL BLAS.
- Fixed Python import errors on Linux caused by `libstdc++.so.6` by
  prompting users during
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  to update `R_LD_LIBRARY_PATH` (automatically or manually) to
  prioritize the Python environment’s `lib`.
- Fixed CRAN check error caused by missing vignette images.
- Resolved compatibility issues when importing emulator objects saved
  from older releases.
- Fixed a bug in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md)
  affecting MAP estimation when `prior = "inv_ga"`.  
- Added support for MAP estimation with no prior in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) by
  setting `prior = NULL`.
- Reverted to the `accelerate` BLAS on Apple Silicon (for all macOS
  versions) due to performance issues observed with the `newaccelerate`
  BLAS using the Vecchia approximation.

## dgpsi 2.6.0

CRAN release: 2025-10-15

- Prediction speed with
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  enhanced for small testing data sets by reducing overhead caused by
  the multi-threading implementation.
- The Python environment now installs packages exclusively from
  conda-forge whenever possible. Packages from other channels will only
  be used if they are unavailable on conda-forge.
- A bug in
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md),
  affecting a bundle of emulators that includes GP emulators, has now
  been fixed.
- The column names from the training input and output provided to
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) are
  retained in the relevant slots of the returned objects, as well as in
  any updated objects produced by the downstream functions that operate
  on them.
- The column names from the testing input and output supplied to
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  and
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  are retained in the relevant slots of the returned objects.
- Improved numerical stability and achieved ~30x faster speed for DGP
  emulators using heteroskedastic likelihoods with replicates, with or
  without the Vecchia approximation.
- Enhanced initialization of DGP emulators with heteroskedastic and
  categorical likelihoods for improved performance.
- Removed the `mode` argument from
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  for DGP emulators with categorical likelihoods. Predictions of class
  probabilities can now be obtained using either the `"mean_var"` or
  `"sampling"` method.
- Set the default `method` for
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md),
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md),
  and
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md) to
  `"mean_var"`.
- Redesigned the output of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  for `dgp` objects with `likelihood = "Categorical"`. See
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  documentation for details.
- Added support for the `nugget_est` argument in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) to
  control whether nuggets of GP nodes feeding into the likelihood node
  are estimated when `likelihood` is not `NULL`.
- Updated initial nugget values when `nugget_est = TRUE` in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md). If
  `likelihood = NULL`, all initial GP nuggets default to `1e-6`;
  otherwise, GP nodes feeding into the likelihood node default to `1e-4`
  and all others to `1e-6`.
- Added the `accuracy` metric to the figures produced by
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  for DGP emulators with categorical likelihoods.
- Fixed the confusion matrix visualization (`style = 2` in
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md))
  so that the diagonal is drawn from top-left to bottom-right.
- Updated
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  to handle errors related to TOS acceptance when installing Miniconda,
  and to automate TOS acceptance for required channels.
- Enabled use of the `newaccelerate` BLAS library on Apple Silicon
  machines running macOS \> 13.3.
- Added the `decouple` argument to
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) to
  allow likelihood parameters to be modeled using separate deep Gaussian
  process hierarchies when `depth > 2`.
- Added the `link` argument to
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) to
  support binary classification using either logit or probit link
  function when `likelihood = "Categorical"`.
- Inference for (D)GPs with homogeneous noise and replicates in the
  training data has been significantly enhanced, achieving over 10×
  speed-up.

## dgpsi 2.5.0

CRAN release: 2024-12-14

- Training times for DGP emulators are now approximately 30%-40% faster.
- The computation of (D)GP predictions and Leave-One-Out (LOO)
  evaluations is now 6-7 times faster.
- The `nb_parallel` argument has been removed from relevant functions,
  as multi-threading is now integrated by default.
- A Vecchia approximation, implemented under the SI framework, is now
  available across all functions to support large-scale emulations.
- Two new functions,
  [`get_thread_num()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/get_thread_num.md)
  and
  [`set_thread_num()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_thread_num.md),
  allow users to inspect and adjust the number of threads used for
  multi-threaded computations.
- A new function,
  [`set_vecchia()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_vecchia.md),
  enables users to easily add or remove the Vecchia approximation for
  GP, DGP, or linked (D)GP emulators.
- Documentation now includes lifecycle status badges to highlight
  deprecated and newly introduced functions and arguments.
- The default value of the `nugget` parameter in DGP emulators with
  likelihood layers has been adjusted from `1e-6` to `1e-4`.
- A `Categorical` likelihood option has been added to the
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md)
  function’s `likelihood` argument, enabling DGP-based classification.
- An issue related to the `LD_LIBRARY` environment variable on Linux
  systems has been resolved via the
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  function.
- The [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md)
  function has been enhanced to accept connection information among
  emulators in the form of a data frame, streamlining linked emulation
  setup.
- A new function,
  [`set_id()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_id.md),
  allows users to assign unique IDs to emulators.
- The
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  function has been updated to accommodate predictions from DGP
  classifiers.  
- The
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  function has been updated to generate validation plots for DGP
  classifiers (i.e., DGP emulators with categorical likelihoods) and
  linked emulators created by
  [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md)
  using the new data frame form for `struc`.  
- The
  [`summary()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  function has been redesigned to provide both summary tables and
  visualizations of structure and model specifications for (D)GP and
  linked (D)GP emulators.  
- A `sample_size` argument has been added to the
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  and
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  functions, allowing users to adjust the number of samples used for
  validation when the validation method is set to `sampling`.
- [`combine()`](https://dplyr.tidyverse.org/reference/combine.html) and
  `set_linked_idx()` are deprecated as of this version and will be
  removed in the next release. These two functions are no longer
  maintained. Please refer to the updated package documentation for
  alternative workflows.
- The basic node functions
  [`kernel()`](https://rdrr.io/r/stats/kernel.html), `Hetero()`,
  `Poisson()`, and `NegBin()`, along with the `struc` argument in the
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md)
  functions, have been removed as of this version. Customization of
  (D)GP specifications can be achieved by modifying the other arguments
  in [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
- The
  [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md)
  function has been updated for instances of the `bundle` class to allow
  drawing of design and evaluation plots of all emulators in a single
  figure.  
- The
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  function has been updated for linked emulators generated by
  [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md)
  using the new data frame form for `struc`.  
- The
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  function has been redesigned to allow new specifications of the
  user-supplied `method` function.  
- The `batch_size` argument has been added to
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  to enable locating multiple design points in a single iteration of the
  sequential design. This argument is compatible with all built-in
  `method` functions:
  [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md),
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md),
  and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md).
- The [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md)
  and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md)
  functions have been redesigned to support continuous search for the
  next design point or search from a discrete candidate set passed
  through the `x_cand` argument.
- The [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md),
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md),
  and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md)
  functions have been updated to output the locations of identified
  design points when a discrete candidate set is not supplied.  
- The `pei()` function has been removed from the package for
  re-engineering and will be added back in a future version.  
- The default of the `refit` argument in the
  [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  function has been changed from `FALSE` to `TRUE`.  
- The
  [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
  function now allows `light = TRUE` for both GP emulators and bundles
  of GP emulators.  
- Two new functions,
  [`serialize()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/serialize.md)
  and
  [`deserialize()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/deserialize.md),
  have been added to allow users to export emulators to multi-session
  workers for parallel processing.
- Additional vignettes are available, showcasing large-scale DGP
  emulation, DGP classification, and Bayesian optimization using (D)GP
  emulators.
- Enhanced clarity and consistency across the documentation.
- Improved examples and explanations in vignettes for better user
  guidance.

## dgpsi 2.4.0

CRAN release: 2024-01-14

- One can now use
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  to implement sequential designs using `f` and a fixed candidate set
  passed to `x_cand` with `y_cand = NULL`.
- The sizes of `.pkl` files written by
  [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
  are significantly reduced.
- One can now set different kernel functions to nodes in different
  layers in a DGP emulator by passing a vector of kernel function names
  to `name` argument of
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
- The default number of imputations `B` in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) and
  [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md) is
  changed to `10` for faster validations and predictions.
- The default method for sequential designs in
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  is changed to
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md).
- A new argument `new_wave` is added to
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  to allow users to resume sequential designs with or without a separate
  wave.
- A bug in
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md) is
  fixed when `object` is an instance of the `bundle` class and
  `batch_size` is greater than one.
- Static and dynamic pruning of DGP structures are implemented in
  [`prune()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/prune.md)
  and
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  (via the new arguments `pruning` and `control`) respectively.
- Some redundant codes are removed from
  [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  which makes
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  slightly faster.
- `limits` argument in
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  is now required when `x_cand` is not supplied to avoid under-sampling
  using the limits inferred from the training data.
- [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  now supports `f` that produce `NA` as outputs. This is useful to
  prevent the sequential design from stopping due to errors or `NA`
  outputs from a simulator at some input locations identified by the
  sequential design process.
- A bug is fixed in
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  when `x_cand` is supplied and the input dimension is one.
- [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md),
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md),
  `pei()`, and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md)
  now accept separate candidate sets (even with different number of
  candidate points) via `x_cand` for bundle emulators.
- A slot called `id` is added to instances of `gp`, `dgp`, `lgp`, and
  `bundle` classes to uniquely identify the emulators. `id` can also be
  passed to instances of `gp`, `dgp`,`lgp`, and `bundle` classes by the
  new `id` argument in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md),
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md),
  [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md), and
  [`pack()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/pack.md).
- [`pack()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/pack.md)
  can now accept a list of (D)GP emulators as the input.
- The `check_point` argument is removed from
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  and replaced by `autosave`.
- Automatic saving of emulators during the sequential design is added to
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  through the new argument `autosave`.
- When a customized evaluation function is provided to
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  via `eval`, the design information in previous waves will be retained
  as long as the previous waves of the sequential design also use
  customized evaluation functions. If different customized evaluation
  functions are supplied to
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  in different waves, the trace plot of RMSEs produced by
  [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md)
  will show RMSEs from different evaluation functions in different
  waves.
- One can now link the same emulator multiple times in a chain via
  [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md) by
  setting different linking information for the emulator via
  `set_linked_idx()`.
- Updates of documentations and vignettes.

## dgpsi 2.3.0

CRAN release: 2023-09-03

- A bug from the underlying Python implementations is fixed when
  `name = 'matern2.5'` in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
- Thanks to [@yyimingucl](https://github.com/yyimingucl), a bug from the
  underlying Python implementations for the MICE sequential design
  criterion
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md) is
  fixed.
- An argument `reset` is added to
  [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  and
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  to reset hyperparameters of a (D)GP emulator to their initial values
  (that were specified when the emulator is initialized) after the input
  and output of the emulator are updated and before the emulator is
  refitted. This argument can be useful for sequential designs in cases
  where the hyperparameters of a (D)GP emulator get caught in suboptimal
  estimates. In such circumstances, one can set `reset = TRUE` to
  reinitialize the (D)GP emulator in some steps of the sequential
  designs as a strategy to escape the poor estimates.
- The refitting of an emulator in the final step of a sequential design
  is no longer forced in
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md).
- An argument `type` is added to
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md) to
  allow users to draw OOS validation plots with testing data shown as a
  line instead of individual points when the emulator’s input is
  one-dimensional and `style = 1`.
- Thanks to [@tjmckinley](https://github.com/tjmckinley), an issue
  relating to `libstdc++.so.6` on Linux machines when R is restarting
  after the installation of the package is fixed.
- [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md) and
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md)
  can locate new design points for stochastic simulators with (D)GP or
  bundle emulators that can deal with stochastic outputs.
- [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  can be used to construct (D)GP or bundle emulators adaptively by
  utilizing multiple realizations from a stochastic simulator at the
  same design positions through the new argument `reps` when
  `method = alm` or `method = mice`.
- A new slot called `specs` is added to the objects returned by
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) that
  contains the key information of the kernel functions used in the
  constructions of GP and DGP emulators.
- Due to a bug in the latest version of an underlying Python package,
  the emulators saved by
  [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
  in version `2.1.6` and `2.2.0` may not work properly with
  [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  and
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  when they are loaded back by
  [`read()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/read.md) in
  this version. This bug has been addressed in this version so emulators
  saved in this version would not have the compatibility issue in future
  version.
- A new sequential design criterion, called the Variance of Improvement
  for Global Fit (VIGF), is added to the package with the function
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md).
- The sampling from an existing candidate set `x_cand` in
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  is changed from a random sampling to a conditioned Latin Hypercube
  sampling in `clhs` package.
- The python environment is now automatically installed or invoked when
  a function from the package is executed. One does not need to run
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  to activate the required python environment but
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  is still useful to re-install and uninstall the underlying python
  environment. A `verb` argument is added to
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  to switch on/off the trace information.

## dgpsi 2.2.0

CRAN release: 2023-06-05

- The efficiency and speed of imputations involved in the training and
  predictions of DGP emulators are significantly improved (achieving
  roughly 3x faster training and imputations) by utilizing blocked Gibbs
  sampling that imputes latent variables layer-wise rather than
  node-wise. The blocked Gibbs sampling is now the default method for
  DGP emulator inference and can be changed back to the old node-wise
  approach by setting `blocked_gibbs = FALSE` in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
- One can now optimize GP components that are contained in the same
  layer of a DGP emulator in parallel during the DGP emulator training,
  using multiple cores by setting the new argument `cores` in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
  This option is useful and can accelerate the training speed when the
  input dimension is moderately large (in which case there is a large
  number of GP components to be optimized) and the optimization of GP
  components is computationally expensive, e.g., when `share = FALSE` in
  which case input dimensions to individual GP components have different
  lengthscales.
- Thanks to [@tjmckinley](https://github.com/tjmckinley), a bug in
  [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  when the `object` is an instance of the `dgp` class (that has been
  trimmed by
  [`window()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/window.md))
  is fixed.
- Thanks to [@tjmckinley](https://github.com/tjmckinley), some R memory
  issues due to the underlying Python implementations are rectified.
- [`set_seed()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_seed.md)
  function is added to ensure reproducible results from the package.
- A bug is fixed when candidate sets `x_cand` and `y_cand` are provided
  to
  [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md).
- One can choose different color palettes using the new argument `color`
  in [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  when `style = 2`.
- `set_linked_idx()` allows constructions of different (D)GP emulators
  (in terms of different connections to the feeding layers) from a same
  (D)GP emulator.

## dgpsi 2.1.6

CRAN release: 2023-02-08

- A bug is found in multi-core predictions in
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  when `object` is an instance of `lgp` class and `x` is a list. This
  bug has been fixed in this version.  
- Thanks to [@tjmckinley](https://github.com/tjmckinley), an issue
  (`/usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'GLIBCXX_3.4.30' not found`)
  encountered in Linux machines is fixed automatically during the
  execution of
  [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md).
- [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md)
  allow users to specify the value of scale parameters and whether to
  estimate the parameters.
- [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) and
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md)
  allow users to specify the bounds of lengthscales.
- The jointly robust prior (Gu, 2019) is implemented as the default
  inference approach for GP emulators in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md).
- The default value of `lengthscale` in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) is
  changed from `0.2` to `0.1`, and the default value for `nugget` in
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) is
  changed from `1e-6` to `1e-8` if `nugget_est = FALSE`.
- One can now specify the number of GP nodes in each layer (except for
  the final layer) of a DGP emulator via the `node` argument in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
- Training data are now contained in the S3 classes `gp` and `dgp`.
- The RMSEs (without the min-max normalization) of emulators are now
  contained in the S3 classes `gp`, `dgp`, and `lgp` after the execution
  of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).
- [`window()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/window.md)
  function is added to trim the traces and obtain new point estimates of
  DGP model parameters for predictions.
- The min-max normalization can now be switched off in
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md) by
  setting the value of `min_max`.
- The default number of imputations `B` for
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) is
  changed from `50` to `30` to better balance the uncertainty and the
  speed of DGP emulator predictions. A new function
  [`set_imp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_imp.md)
  is made available to change the number of imputations of a trained DGP
  emulator so one can either achieve faster predictions by further
  reducing the number of imputations, or account for more imputation
  uncertainties by increasing the number of imputations, without
  re-training the emulator.
- The default number of imputations `B` for
  [`continue()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/continue.md)
  is set to `NULL`, in which case the same number of imputations used in
  `object` will be applied.
- `nugget` argument of
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) now
  specifies the nugget values for GP nodes in different layers rather
  than GP nodes in the final layer.
- The speed of predictions from DGP emulators with squared exponential
  kernels is significantly improved and is roughly 3x faster than the
  implementations in version `2.1.5`.
- The implementation of sequential designs (with two vignettes) of (D)GP
  emulators using different criterion is made available.
- Thanks to [@tjmckinley](https://github.com/tjmckinley), an internal
  reordering issue in
  [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md) is
  fixed.
- [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  now allows users to reinstall and uninstall the underlying Python
  environment.
- A bug that occurs when a linked DGP emulator involves a DGP emulator
  with external inputs is fixed.
- `Intel SVML` will now be installed with the Python environment
  automatically for Intel users for faster implementations.

## dgpsi 2.1.5

CRAN release: 2022-09-29

- Initial release of the R interface to the Python package
  `dgpsi v2.1.5`.
