# Sequential design of a (D)GP emulator or a bundle of (D)GP emulators

This function implements sequential design and active learning for a
(D)GP emulator or a bundle of (D)GP emulators, supporting an array of
popular methods as well as user-specified approaches. It can also be
used as a wrapper for Bayesian optimization methods.

## Usage

``` r
design(
  object,
  N,
  x_cand,
  y_cand,
  n_sample,
  limits,
  f,
  reps,
  freq,
  x_test,
  y_test,
  reset,
  target,
  method,
  batch_size,
  eval,
  verb,
  autosave,
  new_wave,
  M_val,
  cores,
  ...
)

# S3 method for class 'gp'
design(
  object,
  N,
  x_cand = NULL,
  y_cand = NULL,
  n_sample = 200,
  limits = NULL,
  f = NULL,
  reps = 1,
  freq = c(1, 1),
  x_test = NULL,
  y_test = NULL,
  reset = FALSE,
  target = NULL,
  method = vigf,
  batch_size = 1,
  eval = NULL,
  verb = TRUE,
  autosave = list(),
  new_wave = TRUE,
  M_val = 50,
  cores = 1,
  ...
)

# S3 method for class 'dgp'
design(
  object,
  N,
  x_cand = NULL,
  y_cand = NULL,
  n_sample = 200,
  limits = NULL,
  f = NULL,
  reps = 1,
  freq = c(1, 1),
  x_test = NULL,
  y_test = NULL,
  reset = FALSE,
  target = NULL,
  method = vigf,
  batch_size = 1,
  eval = NULL,
  verb = TRUE,
  autosave = list(),
  new_wave = TRUE,
  M_val = 50,
  cores = 1,
  train_N = NULL,
  refit_cores = 1,
  pruning = TRUE,
  control = list(),
  ...
)

# S3 method for class 'bundle'
design(
  object,
  N,
  x_cand = NULL,
  y_cand = NULL,
  n_sample = 200,
  limits = NULL,
  f = NULL,
  reps = 1,
  freq = c(1, 1),
  x_test = NULL,
  y_test = NULL,
  reset = FALSE,
  target = NULL,
  method = vigf,
  batch_size = 1,
  eval = NULL,
  verb = TRUE,
  autosave = list(),
  new_wave = TRUE,
  M_val = 50,
  cores = 1,
  train_N = NULL,
  refit_cores = 1,
  ...
)
```

## Arguments

- object:

  can be one of the following:

  - the S3 class `gp`.

  - the S3 class `dgp`.

  - the S3 class `bundle`.

- N:

  the number of iterations for the sequential design.

- x_cand:

  a matrix (with each row being a design point and column being an input
  dimension) that gives a candidate set from which the next design
  points are determined. Defaults to `NULL`.

- y_cand:

  a matrix (with each row being a simulator evaluation and column being
  an output dimension) that gives the realizations from the simulator at
  input positions in `x_cand`. Defaults to `NULL`.

- n_sample:

  an integer that gives the size of a sub-set to be sampled from the
  candidate set `x_cand` at each step of the sequential design to
  determine the next design point, if `x_cand` is not `NULL`.

  Defaults to `200`.

- limits:

  a two-column matrix that gives the ranges of each input dimension, or
  a vector of length two if there is only one input dimension. If a
  vector is provided, it will be converted to a two-column row matrix.
  The rows of the matrix correspond to input dimensions, and its first
  and second columns correspond to the minimum and maximum values of the
  input dimensions. Set `limits = NULL` if `x_cand` is supplied. This
  argument is only used when `x_cand` is not supplied, i.e.,
  `x_cand = NULL`. Defaults to `NULL`. If you provide a custom `method`
  function with an argument called `limits`, the value of `limits` will
  be passed to your function.

- f:

  an R function representing the simulator. `f` must adhere to the
  following rules:

  - **First argument**: a matrix where rows correspond to different
    design points, and columns represent input dimensions.

  - **Function output**:

    - a matrix where rows correspond to different outputs (matching the
      input design points) and columns represent output dimensions. If
      there is only one output dimension, the function should return a
      matrix with a single column.

    - alternatively, a list where:

      - the first element is the output matrix as described above.

      - additional named elements can optionally update values of
        arguments with matching names passed via `...`. This list output
        is useful if additional arguments to `f`, `method`, or `eval`
        need to be updated after each sequential design iteration.

  See the *Note* section below for additional details. This argument is
  required and must be supplied when `y_cand = NULL`. Defaults to
  `NULL`.

- reps:

  an integer that gives the number of repetitions of the located design
  points to be created and used for evaluations of `f`. Set the argument
  to an integer greater than `1` only if `f` is a stochastic function
  that can generate different responses given for the same input and the
  supplied emulator `object` can deal with stochastic responses, e.g., a
  (D)GP emulator with `nugget_est = TRUE` or a DGP emulator with a
  likelihood layer. The argument is only used when `f` is supplied.
  Defaults to `1`.

- freq:

  a vector of two integers with the first element indicating the number
  of iterations taken between re-estimating the emulator
  hyperparameters, and the second element defining the number of
  iterations to take between re-calculation of evaluating metrics on the
  validation set (see `x_test` below) via the `eval` function. Defaults
  to `c(1, 1)`.

- x_test:

  a matrix (with each row being an input testing data point and each
  column being an input dimension) that gives the testing input data to
  evaluate the emulator after each `freq[2]` iterations of the
  sequential design. Set to `NULL` for LOO-based emulator validation.
  Defaults to `NULL`. This argument is only used if `eval = NULL`.

- y_test:

  the testing output data corresponding to `x_test` for emulator
  validation after each `freq[2]` iterations of the sequential design:

  - if `object` is an instance of the `gp` class, `y_test` is a matrix
    with only one column and each row contains a testing output data
    point from the corresponding row of `x_test`.

  - if `object` is an instance of the `dgp` class, `y_test` is a matrix
    with its rows containing testing output data points corresponding to
    the same rows of `x_test` and columns representing the output
    dimensions.

  - if `object` is an instance of the `bundle` class, `y_test` is a
    matrix with each row representing the outputs for the corresponding
    row of `x_test` and each column representing the output of the
    different emulators in the bundle.

  Set to `NULL` for LOO-based emulator validation. Defaults to `NULL`.
  This argument is only used if `eval = NULL`.

- reset:

  A bool or a vector of bools indicating whether to reset the
  hyperparameters of the emulator(s) to their initial values (as set
  during initial construction) before re-fitting. The re-fitting occurs
  based on the frequency specified by `freq[1]`. This option is useful
  when hyperparameters are suspected to have converged to a local
  optimum affecting validation performance.

  - If a single bool is provided, it applies to every iteration of the
    sequential design.

  - If a vector is provided, its length must equal `N` (even if the
    re-fit frequency specified in `freq[1]` is not 1) and it will apply
    to the corresponding iterations of the sequential design.

  Defaults to `FALSE`.

- target:

  a number or vector specifying the target evaluation metric value(s) at
  which the sequential design should terminate. Defaults to `NULL`, in
  which case the sequential design stops after `N` steps. See the *Note*
  section below for further details about `target`.

- method:

  an R function that determines the next design points to be evaluated
  by `f`. The function must adhere to the following rules:

  - **First argument**: an emulator object, which can be one of the
    following:

    - an instance of the `gp` class (produced by
      [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md));

    - an instance of the `dgp` class (produced by
      [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md));

    - an instance of the `bundle` class (produced by
      [`pack()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/pack.md)).

  - **Second argument** (if `x_cand` is not `NULL`): a *candidate
    matrix* representing a set of potential design points from which the
    `method` function selects the next points.

  - **Function output**:

    - If `x_cand` is not `NULL`:

      - for `gp` or `dgp` objects, the output must be a vector of row
        indices corresponding to the selected design points from the
        *candidate matrix* (the second argument).

      - for `bundle` objects, the output must be a matrix containing the
        row indices of the selected design points from the *candidate
        matrix*. Each column corresponds to the indices for an
        individual emulator in the bundle.

    - If `x_cand` is `NULL`:

      - for `gp` or `dgp` objects, the output must be a matrix where
        each row represents a new design point to be added.

      - for `bundle` objects, the output must be a list with a length
        equal to the number of emulators in the bundle. Each element in
        the list is a matrix where rows represent the new design points
        for the corresponding emulator.

  See [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md),
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md),
  and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md)
  for examples of built-in `method` functions. Defaults to
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md).

- batch_size:

  an integer specifying the number of design points to select in a
  single iteration. Defaults to `1`. This argument is used by the
  built-in `method` functions
  [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md),
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md),
  and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md).
  If you provide a custom `method` function with an argument named
  `batch_size`, the value of `batch_size` will be passed to your
  function.

- eval:

  an R function that computes a customized metric for evaluating
  emulator performance. The function must adhere to the following rules:

  - **First argument**: an emulator object, which can be one of the
    following:

    - an instance of the `gp` class (produced by
      [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md));

    - an instance of the `dgp` class (produced by
      [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md));

    - an instance of the `bundle` class (produced by
      [`pack()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/pack.md)).

  - **Function output**:

    - for `gp` objects, the output must be a single metric value.

    - for `dgp` objects, the output can be a single metric value or a
      vector of metric values with a length equal to the number of
      output dimensions.

    - for `bundle` objects, the output can be a single metric value or a
      vector of metric values with a length equal to the number of
      emulators in the bundle.

  If no custom function is provided, a built-in evaluation metric (RMSE
  or log-loss, in the case of DGP emulators with categorical
  likelihoods) will be used. Defaults to `NULL`. See the *Note* section
  below for additional details.

- verb:

  a bool indicating if trace information will be printed during the
  sequential design. Defaults to `TRUE`.

- autosave:

  a list that contains configuration settings for the automatic saving
  of the emulator:

  - `switch`: a bool indicating whether to enable automatic saving of
    the emulator during sequential design. When set to `TRUE`, the
    emulator in the final iteration is always saved. Defaults to
    `FALSE`.

  - `directory`: a string specifying the directory path where the
    emulators will be stored. Emulators will be stored in a
    sub-directory of `directory` named 'emulator-`id`'. Defaults to
    './check_points'.

  - `fname`: a string representing the base name for the saved emulator
    files. Defaults to 'check_point'.

  - `save_freq`: an integer indicating the frequency of automatic saves,
    measured in the number of iterations. Defaults to `5`.

  - `overwrite`: a bool value controlling the file saving behavior. When
    set to `TRUE`, each new automatic save overwrites the previous one,
    keeping only the latest version. If `FALSE`, each automatic save
    creates a new file, preserving all previous versions. Defaults to
    `FALSE`.

- new_wave:

  a bool indicating whether the current call to `design()` will create a
  new wave of sequential designs or add the next sequence of designs to
  the most recent wave. This argument is relevant only if waves already
  exist in the emulator. Creating new waves can improve the
  visualization of sequential design performance across different calls
  to `design()` via
  [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md),
  and allows for specifying a different evaluation frequency in `freq`.
  However, disabling this option can help limit the number of waves
  visualized in
  [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md) to
  avoid issues such as running out of distinct colors for large numbers
  of waves. Defaults to `TRUE`.

- M_val:

  an integer that gives the size of the conditioning set for the Vecchia
  approximation in emulator validations. This argument is only used if
  the emulator `object` was constructed under the Vecchia approximation.
  Defaults to `50`.

- cores:

  an integer that gives the number of processes to be used for emulator
  validation. If set to `NULL`, the number of processes is set to
  `max physical cores available %/% 2`. Defaults to `1`. This argument
  is only used if `eval = NULL`.

- ...:

  Any arguments with names that differ from those used in `design()` but
  are required by `f`, `method`, or `eval` can be passed here.
  `design()` will forward relevant arguments to `f`, `method`, and
  `eval` based on the names of the additional arguments provided.

- train_N:

  the number of training iterations to be used for re-fitting the DGP
  emulator at each step of the sequential design:

  - If `train_N` is an integer, the DGP emulator will be re-fitted at
    each step (based on the re-fit frequency specified in `freq[1]`)
    using `train_N` iterations.

  - If `train_N` is a vector, its length must be `N`, even if the re-fit
    frequency specified in `freq[1]` is not 1.

  - If `train_N` is `NULL`, the DGP emulator will be re-fitted at each
    step (based on the re-fit frequency specified in `freq[1]`) using:

    - `100` iterations if the DGP emulator was constructed without the
      Vecchia approximation, or

    - `50` iterations if the Vecchia approximation was used.

  Defaults to `NULL`.

- refit_cores:

  the number of processes to be used to re-fit GP components (in the
  same layer of a DGP emulator) at each M-step during the re-fitting. If
  set to `NULL`, the number of processes is set to
  `(max physical cores available - 1)` if the DGP emulator was
  constructed without the Vecchia approximation. Otherwise, the number
  of processes is set to `max physical cores available %/% 2`. Only use
  multiple processes when there is a large number of GP components in
  different layers and optimization of GP components is computationally
  expensive. Defaults to `1`.

- pruning:

  a bool indicating if dynamic pruning of DGP structures will be
  implemented during the sequential design after the total number of
  design points exceeds `min_size` in `control`. The argument is only
  applicable to DGP emulators (i.e., `object` is an instance of `dgp`
  class) produced by
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
  Defaults to `TRUE`.

- control:

  a list that can supply any of the following components to control the
  dynamic pruning of the DGP emulator:

  - `min_size`, the minimum number of design points required to trigger
    dynamic pruning. Defaults to 10 times the number of input
    dimensions.

  - `threshold`, the \\R^2\\ value above which a GP node is considered
    redundant. Defaults to `0.97`.

  - `nexceed`, the minimum number of consecutive iterations that the
    \\R^2\\ value of a GP node must exceed `threshold` to trigger the
    removal of that node from the DGP structure. Defaults to `3`.

  The argument is only used when `pruning = TRUE`.

## Value

An updated `object` is returned with a slot called `design` that
contains:

- *S* slots, named `wave1, wave2,..., waveS`, that contain information
  of *S* waves of sequential design that have been applied to the
  emulator. Each slot contains the following elements:

  - `N`, an integer that gives the numbers of iterations implemented in
    the corresponding wave;

  - `rmse`, a matrix providing the evaluation metric values for
    emulators constructed during the corresponding wave, when
    `eval = NULL`. Each row of the matrix represents an iteration.

    - for an `object` of class `gp`, the matrix contains a single column
      of RMSE values.

    - for an `object` of class `dgp` without a categorical likelihood,
      each row contains mean/median squared errors corresponding to
      different output dimensions.

    - for an `object` of class `dgp` with a categorical likelihood, the
      matrix contains a single column of log-loss values.

    - for an `object` of class `bundle`, each row contains either
      mean/median squared errors or log-loss values for the emulators in
      the bundle.

  - `metric`: a matrix providing the values of custom evaluation
    metrics, as computed by the user-supplied `eval` function, for
    emulators constructed during the corresponding wave.

  - `freq`, an integer that gives the frequency that the emulator
    validations are implemented during the corresponding wave.

  - `enrichment`, a vector of size `N` that gives the number of new
    design points added after each step of the sequential design (if
    `object` is an instance of the `gp` or `dgp` class), or a matrix
    that gives the number of new design points added to emulators in a
    bundle after each step of the sequential design (if `object` is an
    instance of the `bundle` class).

  If `target` is not `NULL`, the following additional elements are also
  included:

  - `target`: the target evaluating metric computed by the `eval` or
    built-in function to stop the sequential design.

  - `reached`: indicates whether the `target` was reached at the end of
    the sequential design:

    - a bool if `object` is an instance of the `gp` or `dgp` class.

    - a vector of bools if `object` is an instance of the `bundle`
      class, with its length determined as follows:

      - equal to the number of emulators in the bundle when
        `eval = NULL`.

      - equal to the length of the output from `eval` when a custom
        `eval` function is provided.

- a slot called `type` that gives the type of validation:

  - either LOO ('loo') or OOS ('oos') if `eval = NULL`. See
    [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
    for more information about LOO and OOS.

  - 'customized' if a customized R function is provided to `eval`.

- two slots called `x_test` and `y_test` that contain the data points
  for the OOS validation if the `type` slot is 'oos'.

- If `y_cand = NULL` and `x_cand` is supplied, and there are `NA`s
  returned from the supplied `f` during the sequential design, a slot
  called `exclusion` is included that records the located design
  positions that produced `NA`s via `f`. The sequential design will use
  this information to avoid re-visiting the same locations in later runs
  of `design()`.

See *Note* section below for further information.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- Validation of an emulator is forced after the final step of a
  sequential design even if `N` is not a multiple of the second element
  in `freq`.

- Any `loo` or `oos` slot that already exists in `object` will be
  cleaned, and a new slot called `loo` or `oos` will be created in the
  returned object depending on whether `x_test` and `y_test` are
  provided. The new slot gives the validation information of the
  emulator constructed in the final step of the sequential design. See
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  for more information about the slots `loo` and `oos`.

- If `object` has previously been used by `design()` for sequential
  design, the information of the current wave of the sequential design
  will replace those of old waves and be contained in the returned
  object, unless

  - the validation type (LOO or OOS depending on whether `x_test` and
    `y_test` are supplied or not) of the current wave of the sequential
    design is the same as the validation types (shown in the `type` of
    the `design` slot of `object`) in previous waves, and if the
    validation type is OOS, `x_test` and `y_test` in the current wave
    must also be identical to those in the previous waves;

  - both the current and previous waves of the sequential design supply
    customized evaluation functions to `eval`. Users need to ensure the
    customized evaluation functions are consistent among different
    waves. Otherwise, the trace plot of RMSEs produced by
    [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md)
    will show values of different evaluation metrics in different waves.

  For the above two cases, the information of the current wave of the
  sequential design will be added to the `design` slot of the returned
  object under the name `waveS`.

- If `object` is an instance of the `gp` class and `eval = NULL`, the
  matrix in the `rmse` slot is single-columned. If `object` is an
  instance of the `dgp` or `bundle` class and `eval = NULL`, the matrix
  in the `rmse` slot can have multiple columns that correspond to
  different output dimensions or different emulators in the bundle.

- If `object` is an instance of the `gp` class and `eval = NULL`,
  `target` needs to be a single value giving the RMSE threshold. If
  `object` is an instance of the `dgp` or `bundle` class and
  `eval = NULL`, `target` can be a vector of values that gives the
  thresholds of evaluating metrics for different output dimensions or
  different emulators. If a single value is provided, it will be used as
  the threshold for all output dimensions (if `object` is an instance of
  the `dgp`) or all emulators (if `object` is an instance of the
  `bundle`). If a customized function is supplied to `eval` and `target`
  is given as a vector, the user needs to ensure that the length of
  `target` is equal to that of the output from `eval`.

- When defining `f`, it is important to ensure that:

  - the column order of the first argument of `f` is consistent with the
    training input used for the emulator;

  - the column order of the output matrix of `f` is consistent with the
    order of emulator output dimensions (if `object` is an instance of
    the `dgp` class), or the order of emulators placed in `object` (if
    `object` is an instance of the `bundle` class).

- The output matrix produced by `f` may include `NA`s. This is
  especially beneficial as it allows the sequential design process to
  continue without interruption, even if errors or `NA` outputs are
  encountered from `f` at certain input locations identified by the
  sequential design. Users should ensure that any errors within `f` are
  handled by appropriately returning `NA`s.

- When defining `eval`, the output metric needs to be positive if
  [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md) is
  used with `log = T`. And one needs to ensure that a lower metric value
  indicates a better emulation performance if `target` is set.

## Examples

``` r
if (FALSE) { # \dontrun{

# load packages and the Python env
library(lhs)
library(dgpsi)

# construct a 2D non-stationary function that takes a matrix as the input
f <- function(x) {
  sin(1/((0.7*x[,1,drop=F]+0.3)*(0.7*x[,2,drop=F]+0.3)))
}

# generate the initial design
X <- maximinLHS(5,2)
Y <- f(X)

# generate the validation data
validate_x <- maximinLHS(30,2)
validate_y <- f(validate_x)

# training a 2-layered DGP emulator with the initial design
m <- dgp(X, Y)

# specify the ranges of the input dimensions
lim_1 <- c(0, 1)
lim_2 <- c(0, 1)
lim <- rbind(lim_1, lim_2)

# 1st wave of the sequential design with 10 steps
m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y)

# 2nd wave of the sequential design with 10 steps
m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y)

# 3rd wave of the sequential design with 10 steps
m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y)

# draw the design created by the sequential design
draw(m,'design')

# inspect the trace of RMSEs during the sequential design
draw(m,'rmse')

# reduce the number of imputations for faster OOS
m_faster <- set_imp(m, 5)

# plot the OOS validation with the faster DGP emulator
plot(m_faster, x_test = validate_x, y_test = validate_y)
} # }
```
