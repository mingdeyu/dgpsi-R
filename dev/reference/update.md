# Update a GP or DGP emulator

This function updates the training input and output of a GP or DGP
emulator with an option to refit the emulator.

## Usage

``` r
update(object, X, Y, refit, reset, verb, ...)

# S3 method for class 'dgp'
update(
  object,
  X,
  Y,
  refit = TRUE,
  reset = FALSE,
  verb = TRUE,
  N = NULL,
  cores = 1,
  ess_burn = 10,
  B = NULL,
  ...
)

# S3 method for class 'gp'
update(object, X, Y, refit = TRUE, reset = FALSE, verb = TRUE, ...)
```

## Arguments

- object:

  can be one of the following:

  - the S3 class `gp`.

  - the S3 class `dgp`.

- X:

  the new input data which is a matrix where each row is an input
  training data point and each column represents an input dimension.

- Y:

  the new output data:

  - If `object` is an instance of the `gp` class, `Y` is a matrix with
    only one column and each row being an output data point.

  - If `object` is an instance of the `dgp` class, `Y` is a matrix with
    its rows being output data points and columns being output
    dimensions. When `likelihood` (see below) is not `NULL`, `Y` must be
    a matrix with only one column.

- refit:

  a bool indicating whether to re-fit the emulator `object` after the
  training input and output are updated. Defaults to `TRUE`.

- reset:

  a bool indicating whether to reset hyperparameters of the emulator
  `object` to the initial values first obtained when the emulator was
  constructed. Use if it is suspected that a local mode for the
  hyperparameters has been reached through successive updates. Defaults
  to `FALSE`.

- verb:

  a bool indicating if trace information will be printed during the
  function execution. Defaults to `TRUE`.

- ...:

  N/A.

- N:

  number of training iterations used to re-fit the emulator `object` if
  it is an instance of the `dgp` class. If set to `NULL`, the number of
  iterations is set to `100` if the DGP emulator was constructed without
  the Vecchia approximation, and is set to `50` if Vecchia approximation
  was used. Defaults to `NULL`.

- cores:

  the number of processes to be used to re-fit GP components (in the
  same layer) at each M-step during the re-fitting. If set to `NULL`,
  the number of processes is set to `(max physical cores available - 1)`
  if `vecchia = FALSE` and `max physical cores available %/% 2` if
  `vecchia = TRUE`. Only use multiple processes when there is a large
  number of GP components in different layers and optimization of GP
  components is computationally expensive. Defaults to `1`.

- ess_burn:

  number of burnin steps for the ESS-within-Gibbs sampler at each I-step
  of the training of the emulator `object` if it is an instance of the
  `dgp` class. Defaults to `10`.

- B:

  the number of imputations for predictions from the updated emulator
  `object` if it is an instance of the `dgp` class. This overrides the
  number of imputations set in `object`. Set to `NULL` to use the same
  number of imputations set in `object`. Defaults to `NULL`.

## Value

An updated `object`.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- The following slots:

  - `loo` and `oos` created by
    [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md);

  - `results` created by
    [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md);
    and

  - `design` created by
    [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)

  in `object` will be removed and not contained in the returned object.

## Examples

``` r
if (FALSE) { # \dontrun{

# See alm(), mice(), or vigf() for an example.
} # }
```
