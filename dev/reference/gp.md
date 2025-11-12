# Gaussian process emulator construction

**\[updated\]**

This function builds and trains a GP emulator.

## Usage

``` r
gp(
  X,
  Y,
  name = "sexp",
  lengthscale = rep(0.1, ncol(X)),
  bounds = NULL,
  prior = "ref",
  nugget_est = FALSE,
  nugget = ifelse(nugget_est, 0.01, 1e-08),
  scale_est = TRUE,
  scale = 1,
  training = TRUE,
  verb = TRUE,
  check_rep = TRUE,
  vecchia = FALSE,
  M = 25,
  ord = NULL,
  id = NULL
)
```

## Arguments

- X:

  a matrix where each row is an input data point and each column is an
  input dimension.

- Y:

  a matrix with only one column and each row being an output data point.

- name:

  kernel function to be used. Either `"sexp"` for squared exponential
  kernel or `"matern2.5"` for Matérn-2.5 kernel. Defaults to `"sexp"`.

- lengthscale:

  initial values of lengthscales in the kernel function. It can be a
  single numeric value or a vector of length `ncol(X)`:

  - if it is a single numeric value, it is assumed that kernel functions
    across input dimensions share the same lengthscale;

  - if it is a vector, it is assumed that kernel functions across input
    dimensions have different lengthscales.

  Defaults to a vector of `0.1`.

- bounds:

  the lower and upper bounds of lengthscales in the kernel function. It
  is a vector of length two where the first element is the lower bound
  and the second element is the upper bound. The bounds will be applied
  to all lengthscales in the kernel function. Defaults to `NULL` where
  no bounds are specified for the lengthscales.

- prior:

  **\[updated\]** prior to be used for Maximum a Posterior (MAP)
  estimation of the GP lengthscales and nugget: accepts no prior
  (`NULL`), a gamma prior (`"ga"`), an inverse-gamma prior (`"inv_ga"`),
  or a jointly robust reference prior (`"ref"`). Defaults to `"ref"`.
  See the reference below for the jointly robust prior.

- nugget_est:

  a bool indicating if the nugget term is to be estimated:

  1.  `FALSE`: the nugget term is fixed to `nugget`.

  2.  `TRUE`: the nugget term will be estimated.

  Defaults to `FALSE`.

- nugget:

  the initial nugget value. If `nugget_est = FALSE`, the assigned value
  is fixed during the training. Set `nugget` to a small value (e.g.,
  `1e-8`) and the corresponding bool in `nugget_est` to `FALSE` for
  deterministic computer models where the emulator should interpolate
  the training data points. Set `nugget` to a larger value and the
  corresponding bool in `nugget_est` to `TRUE` for stochastic emulation
  where the computer model outputs are assumed to follow a homogeneous
  Gaussian distribution. Defaults to `1e-8` if `nugget_est = FALSE` and
  `0.01` if `nugget_est = TRUE`.

- scale_est:

  a bool indicating if the variance is to be estimated:

  1.  `FALSE`: the variance is fixed to `scale`.

  2.  `TRUE`: the variance term will be estimated.

  Defaults to `TRUE`.

- scale:

  the initial variance value. If `scale_est = FALSE`, the assigned value
  is fixed during the training. Defaults to `1`.

- training:

  a bool indicating if the initialized GP emulator will be trained. When
  set to `FALSE`, `gp()` returns an untrained GP emulator, to which one
  can apply
  [`summary()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  to inspect its specification or apply
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  to check its emulation performance before the training. Defaults to
  `TRUE`.

- verb:

  a bool indicating if the trace information on GP emulator construction
  and training will be printed during function execution. Defaults to
  `TRUE`.

- check_rep:

  **\[new\]** a bool indicating whether to check for repetitions in the
  dataset, i.e., if one input position has multiple outputs. Defaults to
  `TRUE`.

- vecchia:

  a bool indicating whether to use Vecchia approximation for large-scale
  GP emulator construction and prediction. Defaults to `FALSE`. The
  Vecchia approximation implemented for the GP emulation largely follows
  Katzfuss et al. (2022). See reference below.

- M:

  the size of the conditioning set for the Vecchia approximation in the
  GP emulator training. Defaults to `25`.

- ord:

  an R function that returns the ordering of the input to the GP
  emulator for the Vecchia approximation. The function must satisfy the
  following basic rules:

  - the first argument represents the input scaled by the lengthscales.

  - the output of the function is a vector of indices that gives the
    ordering of the input to the GP emulator.

  If `ord = NULL`, the default random ordering is used. Defaults to
  `NULL`.

- id:

  an ID to be assigned to the GP emulator. If an ID is not provided
  (i.e., `id = NULL`), a UUID (Universally Unique Identifier) will be
  automatically generated and assigned to the emulator. Default to
  `NULL`.

## Value

An S3 class named `gp` that contains five slots:

- `id`: A number or character string assigned through the `id` argument.

- `data`: a list that contains two elements: `X` and `Y` which are the
  training input and output data respectively.

- `specs`: a list that contains seven elements:

  1.  `kernel`: the type of the kernel function used. Either `"sexp"`
      for squared exponential kernel or `"matern2.5"` for Matérn-2.5
      kernel.

  2.  `lengthscales`: a vector of lengthscales in the kernel function.

  3.  `scale`: the variance value in the kernel function.

  4.  `nugget`: the nugget value in the kernel function.

  5.  `vecchia`: whether the Vecchia approximation is used for the GP
      emulator training.

  6.  `M`: the size of the conditioning set for the Vecchia
      approximation in the GP emulator training.

- `constructor_obj`: a 'python' object that stores the information of
  the constructed GP emulator.

- `container_obj`: a 'python' object that stores the information for the
  linked emulation.

- `emulator_obj`: a 'python' object that stores the information for the
  predictions from the GP emulator.

The returned `gp` object can be used by

- [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  for GP predictions.

- [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  for LOO and OOS validations.

- [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  for validation plots.

- [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md) for
  linked (D)GP emulator constructions.

- [`summary()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  to summarize the trained GP emulator.

- [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
  to save the GP emulator to a `.pkl` file.

- [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  for sequential designs.

- [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  to update the GP emulator with new inputs and outputs.

- [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md),
  [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md),
  and
  [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md) to
  locate next design points.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

Any R vector detected in `X` and `Y` will be treated as a column vector
and automatically converted into a single-column R matrix. Thus, if `X`
is a single data point with multiple dimensions, it must be given as a
matrix.

## References

- Gu, M. (2019). Jointly robust prior for Gaussian stochastic process in
  emulation, calibration and variable selection. *Bayesian Analysis*,
  **14(3)**, 857-885.

- Katzfuss, M., Guinness, J., & Lawrence, E. (2022). Scaled Vecchia
  approximation for fast computer-model emulation. *SIAM/ASA Journal on
  Uncertainty Quantification*, **10(2)**, 537-554.

## Examples

``` r
if (FALSE) { # \dontrun{
# load the package and the Python env
library(dgpsi)

# construct a step function
f <- function(x) {
   if (x < 0.5) return(-1)
   if (x >= 0.5) return(1)
  }

# generate training data
X <- seq(0, 1, length = 10)
Y <- sapply(X, f)

# training
m <- gp(X, Y)

# summarizing
summary(m)

# LOO cross validation
m <- validate(m)
plot(m)

# prediction
test_x <- seq(0, 1, length = 200)
m <- predict(m, x = test_x)

# OOS validation
validate_x <- sample(test_x, 10)
validate_y <- sapply(validate_x, f)
plot(m, validate_x, validate_y)

# write and read the constructed emulator
write(m, 'step_gp')
m <- read('step_gp')
} # }
```
