# Deep Gaussian process emulator construction

**\[updated\]**

This function builds and trains a DGP emulator.

## Usage

``` r
dgp(
  X,
  Y,
  depth = 2,
  node = ncol(X),
  name = "sexp",
  lengthscale = 1,
  bounds = NULL,
  prior = "ga",
  share = TRUE,
  nugget_est = FALSE,
  nugget = NULL,
  scale_est = TRUE,
  scale = 1,
  connect = NULL,
  likelihood = NULL,
  training = TRUE,
  verb = TRUE,
  check_rep = TRUE,
  vecchia = FALSE,
  M = 25,
  ord = NULL,
  N = ifelse(vecchia, 200, 500),
  cores = 1,
  blocked_gibbs = TRUE,
  ess_burn = 10,
  burnin = NULL,
  B = 10,
  id = NULL,
  decouple = FALSE,
  link = NULL
)
```

## Arguments

- X:

  a matrix where each row is an input training data point and each
  column represents an input dimension.

- Y:

  a matrix containing observed training output data. The matrix has its
  rows being output data points and columns representing output
  dimensions. When `likelihood` (see below) is not `NULL`, `Y` must be a
  matrix with a single column.

- depth:

  number of layers (including the likelihood layer) for a DGP structure.
  `depth` must be at least `2`. Defaults to `2`.

- node:

  number of GP nodes in each layer (except for the final layer or the
  layer feeding the likelihood node) of the DGP. Defaults to `ncol(X)`.

- name:

  a character or a vector of characters that indicates the kernel
  functions (either `"sexp"` for squared exponential kernel or
  `"matern2.5"` for Matérn-2.5 kernel) used in the DGP emulator:

  1.  if a single character is supplied, the corresponding kernel
      function will be used for all GP nodes in the DGP hierarchy.

  2.  if a vector of characters is supplied, each character of the
      vector specifies the kernel function that will be applied to all
      GP nodes in the corresponding layer.

  Defaults to `"sexp"`.

- lengthscale:

  initial lengthscales for GP nodes in the DGP emulator. It can be a
  single numeric value or a vector:

  1.  if it is a single numeric value, the value will be applied as the
      initial lengthscales for all GP nodes in the DGP hierarchy.

  2.  if it is a vector, each element of the vector specifies the
      initial lengthscales that will be applied to all GP nodes in the
      corresponding layer. The vector should have a length of `depth` if
      `likelihood = NULL` or a length of `depth - 1` if `likelihood` is
      not `NULL`.

  Defaults to a numeric value of `1.0`.

- bounds:

  the lower and upper bounds of lengthscales in GP nodes. It can be a
  vector or a matrix:

  1.  if it is a vector, the lower bound (the first element of the
      vector) and upper bound (the second element of the vector) will be
      applied to lengthscales for all GP nodes in the DGP hierarchy.

  2.  if it is a matrix, each row of the matrix specifies the lower and
      upper bounds of lengthscales for all GP nodes in the corresponding
      layer. The matrix should have its row number equal to `depth` if
      `likelihood = NULL` or to `depth - 1` if `likelihood` is not
      `NULL`.

  Defaults to `NULL` where no bounds are specified for the lengthscales.

- prior:

  **\[updated\]** prior to be used for MAP estimation of lengthscales
  and nuggets of all GP nodes in the DGP hierarchy:

  - no prior (`NULL`),

  - gamma prior (`"ga"`),

  - inverse gamma prior (`"inv_ga"`), or

  - jointly robust prior (`"ref"`).

  Defaults to `"ga"`.

- share:

  a bool indicating if all input dimensions of a GP node share a common
  lengthscale. Defaults to `TRUE`.

- nugget_est:

  **\[updated\]** a bool or a bool vector indicating whether the nuggets
  of GP nodes in the final layer (or the layer feeding the likelihood
  node) should be estimated. If a bool is provided, it is applied to all
  GP nodes in that layer. If a bool vector is provided, its length must
  match the number of GP nodes:

  - `ncol(Y)` if `likelihood = NULL`

  - `3` if `likelihood` is `"ZINB"`

  - `2` if `likelihood` is `"Hetero"` or `"NegBin"` or `"ZIP"`

  - `1` if `likelihood` is `"Poisson"` or `"Categorical"` with two
    classes

  - the number of classes if `likelihood` is `"Categorical"` with more
    than two classes.

  Each element of the vector is applied to the corresponding GP node in
  the final layer (or the layer feeding the likelihood node). The value
  of a bool has following effects:

  - `FALSE`: the nugget of the corresponding GP is fixed to the
    corresponding value defined in `nugget` (see below).

  - `TRUE`: the nugget of the corresponding GP will be estimated with
    the initial value given by the correspondence in `nugget` (see
    below).

  Defaults to `FALSE`.

- nugget:

  the initial nugget value(s) of GP nodes (if any) in each layer:

  1.  if it is a single numeric value, the value will be applied as the
      initial nugget for all GP nodes in the DGP hierarchy.

  2.  if it is a vector, each element of the vector specifies the
      initial nugget that will be applied to all GP nodes in the
      corresponding layer. The vector should have a length of `depth` if
      `likelihood = NULL` or a length of `depth - 1` if `likelihood` is
      not `NULL`.

  Set `nugget` to a small value and the bools in `nugget_est` to `FALSE`
  for deterministic emulation, where the emulator interpolates the
  training data points. Set `nugget` to a larger value and the bools in
  `nugget_est` to `TRUE` for stochastic emulation where the computer
  model outputs are assumed to follow a homogeneous Gaussian
  distribution. Defaults to `1e-6` if `likelihood` is `NULL`. If
  `likelihood` is not `NULL`, the nuggets of GPs that feed into the
  likelihood layer default to `1e-4`, while those of all other GPs
  default to `1e-6`.

- scale_est:

  **\[updated\]** a bool or a bool vector indicating whether the
  variances of GP nodes in the final layer (or the layer feeding the
  likelihood node) should be estimated. If a bool is provided, it is
  applied to all GP nodes in that layer. If a bool vector is provided,
  its length must match the number of GP nodes:

  - `ncol(Y)` if `likelihood = NULL`

  - `3` if `likelihood` is `"ZINB"`

  - `2` if `likelihood` is `"Hetero"` or `"NegBin"` or `"ZIP"`

  - `1` if `likelihood` is `"Poisson"` or `"Categorical"` with two
    classes

  - the number of classes if `likelihood` is `"Categorical"` with more
    than two classes.

  The value of a bool has following effects:

  - `FALSE`: the variance of the corresponding GP is fixed to the
    corresponding value defined in `scale` (see below).

  - `TRUE`: the variance of the corresponding GP will be estimated with
    the initial value given by the correspondence in `scale` (see
    below).

  Defaults to `TRUE`.

- scale:

  **\[updated\]** the initial variance value(s) of GP nodes in the final
  layer (or the layer feeding the likelihood node). If it is a single
  numeric value, it will be applied to all GP nodes in the final layer
  (or the layer feeding the likelihood node). If it is a vector, its
  length must match the number of GP nodes:

  - `ncol(Y)` if `likelihood = NULL`

  - `3` if `likelihood` is `"ZINB"`

  - `2` if `likelihood` is `"Hetero"` or `"NegBin"` or `"ZIP"`

  - `1` if `likelihood` is `"Poisson"` or `"Categorical"` with two
    classes

  - the number of classes if `likelihood` is `"Categorical"` with more
    than two classes.

  Each numeric in the vector will be applied to the corresponding GP
  node.

  Defaults to `1`.

- connect:

  a bool indicating whether to apply global input connections in the DGP
  structure. Setting this to `FALSE` may yield a better emulator in some
  cases. When set to `NULL`, the value defaults to `FALSE` if
  `likelihood = "Categorical"` and to `TRUE` otherwise. Defaults to
  `NULL`.

- likelihood:

  **\[updated\]** the likelihood type of a DGP emulator:

  1.  `NULL`: no likelihood layer is included in the emulator.

  2.  `"Hetero"`: a heteroskedastic Gaussian likelihood layer is added
      for stochastic emulation where the computer model outputs are
      assumed to follow a heteroskedastic Gaussian distribution (i.e.,
      the computer model outputs have input-dependent noise).

  3.  `"Poisson"`: a Poisson likelihood layer is added for emulation
      where the computer model outputs are counts and a Poisson
      distribution is used to model them.

  4.  `"NegBin"`: a negative Binomial likelihood layer is added for
      emulation where the computer model outputs are counts and a
      negative Binomial distribution is used to capture dispersion
      variability in input space.

  5.  `"ZIP"`: a zero-inflated Poisson likelihood layer is added for
      emulation where the computer model outputs are counts with excess
      zeros relative to a standard Poisson distribution, modeled via a
      mixture of structural zeros and a Poisson component.

  6.  `"ZINB"`: a zero-inflated negative Binomial likelihood layer is
      added for emulation where the computer model outputs are counts
      exhibiting both over-dispersion and excess zeros, combining a
      structural-zero component with a negative Binomial component.

  7.  `"Categorical"`: a categorical likelihood layer is added for
      emulation (classification), where the computer model output is
      categorical.

  Defaults to `NULL`.

- training:

  a bool indicating if the initialized DGP emulator will be trained.
  When set to `FALSE`, `dgp()` returns an untrained DGP emulator, to
  which one can apply
  [`summary()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  to inspect its specifications or apply
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  to check its emulation performance before training. Defaults to
  `TRUE`.

- verb:

  a bool indicating if the trace information on DGP emulator
  construction and training will be printed during the function
  execution. Defaults to `TRUE`.

- check_rep:

  a bool indicating whether to check for repetitions in the dataset,
  i.e., if one input position has multiple outputs. Defaults to `TRUE`.

- vecchia:

  a bool indicating whether to use Vecchia approximation for large-scale
  DGP emulator construction and prediction. Defaults to `FALSE`.

- M:

  the size of the conditioning set for the Vecchia approximation in the
  DGP emulator training. Defaults to `25`.

- ord:

  an R function that returns the ordering of the input to each GP node
  contained in the DGP emulator for the Vecchia approximation. The
  function must satisfy the following basic rules:

  - the first argument represents the input to a GP node scaled by its
    lengthscales.

  - the output of the function is a vector of indices that gives the
    ordering of the input to the GP node.

  If `ord = NULL`, the default random ordering is used. Defaults to
  `NULL`.

- N:

  number of iterations for the training. Defaults to `500` if
  `vecchia = FALSE` and `200` if `vecchia = TRUE`. This argument is only
  used when `training = TRUE`.

- cores:

  the number of processes to be used to optimize GP components (in the
  same layer) at each M-step of the training. If set to `NULL`, the
  number of processes is set to `(max physical cores available - 1)` if
  `vecchia = FALSE` and `max physical cores available %/% 2` if
  `vecchia = TRUE`. Only use multiple processes when there is a large
  number of GP components in different layers and optimization of GP
  components is computationally expensive. Defaults to `1`.

- blocked_gibbs:

  a bool indicating if the latent variables are imputed layer-wise using
  ESS-within-Blocked-Gibbs. ESS-within-Blocked-Gibbs would be faster and
  more efficient than ESS-within-Gibbs that imputes latent variables
  node-wise because it reduces the number of components to be sampled
  during Gibbs steps, especially when there is a large number of GP
  nodes in layers due to higher input dimensions. Default to `TRUE`.

- ess_burn:

  number of burnin steps for the ESS-within-Gibbs at each I-step of the
  training. Defaults to `10`. This argument is only used when
  `training = TRUE`.

- burnin:

  the number of training iterations to be discarded for point estimates
  of model parameters. Must be smaller than the training iterations `N`.
  If this is not specified, only the last 25% of iterations are used.
  Defaults to `NULL`. This argument is only used when `training = TRUE`.

- B:

  the number of imputations used to produce predictions. Increase the
  value to refine the representation of imputation uncertainty. Defaults
  to `10`.

- id:

  an ID to be assigned to the DGP emulator. If an ID is not provided
  (i.e., `id = NULL`), a UUID (Universally Unique Identifier) will be
  automatically generated and assigned to the emulator. Default to
  `NULL`.

- decouple:

  **\[new\]** A boolean indicating whether the model parameters for the
  heteroskedastic Gaussian likelihood, negative Binomial likelihood, and
  categorical likelihood (when the number of categories is greater
  than 2) should be modeled using separate deep Gaussian process
  hierarchies when `depth` is greater than 2. Defaults to `FALSE`.

- link:

  **\[new\]** The link function used for classification when
  `likelihood = "Categorical"`. Supported options are `"logit"` and
  `"probit"` for binary classification, and `"softmax"` or `"robustmax"`
  for multi-class classification. If set to `NULL`, the default is
  `"logit"` for binary classification and `"softmax"` for multi-class
  classification. Defaults to `NULL`.

## Value

An S3 class named `dgp` that contains five slots:

- `id`: A number or character string assigned through the `id` argument.

- `data`: a list that contains two elements: `X` and `Y` which are the
  training input and output data respectively.

- `specs`: a list that contains

  1.  *L* (i.e., the number of layers in the DGP hierarchy) sub-lists
      named `layer1, layer2,..., layerL`. Each sub-list contains *D*
      (i.e., the number of GP/likelihood nodes in the corresponding
      layer) sub-lists named `node1, node2,..., nodeD`. If a sub-list
      corresponds to a likelihood node, it contains one element called
      `type` that gives the name (`Hetero`, `Poisson`, `NegBin`, `ZIP`,
      `ZINB`, or `Categorical`) of the likelihood node. If a sub-list
      corresponds to a GP node, it contains four elements:

      - `kernel`: the type of the kernel function used for the GP node.

      - `lengthscales`: a vector of lengthscales in the kernel function.

      - `scale`: the variance value in the kernel function.

      - `nugget`: the nugget value in the kernel function.

  2.  `seed`: the random seed generated to produce imputations. This
      information is stored for reproducibility when the DGP emulator
      (that was saved by
      [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
      with the light option `light = TRUE`) is loaded back to R by
      [`read()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/read.md).

  3.  `B`: the number of imputations used to generate the emulator.

  4.  `vecchia`: whether the Vecchia approximation is used for the GP
      emulator training.

  5.  `M`: the size of the conditioning set for the Vecchia
      approximation in the DGP emulator training. `M` is generated only
      when `vecchia = TRUE`.

- `constructor_obj`: a 'python' object that stores the information of
  the constructed DGP emulator.

- `container_obj`: a 'python' object that stores the information for the
  linked emulation.

- `emulator_obj`: a 'python' object that stores the information for the
  predictions from the DGP emulator.

The returned `dgp` object can be used by

- [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  for DGP predictions.

- [`continue()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/continue.md)
  for additional DGP training iterations.

- [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  for LOO and OOS validations.

- [`plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  for validation plots.

- [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md) for
  linked (D)GP emulator constructions.

- [`window()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/window.md)
  for model parameter trimming.

- [`summary()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  to summarize the trained DGP emulator.

- [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
  to save the DGP emulator to a `.pkl` file.

- [`set_imp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_imp.md)
  to change the number of imputations.

- [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  for sequential design.

- [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  to update the DGP emulator with new inputs and outputs.

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

# set a random seed
set_seed(999)

# training a DGP emulator
m <- dgp(X, Y)

# continue for further training iterations
m <- continue(m)

# summarizing
summary(m)

# trace plot
trace_plot(m)

# trim the traces of model parameters
m <- window(m, 800)

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
write(m, 'step_dgp')
m <- read('step_dgp')
} # }
```
