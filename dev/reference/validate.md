# Validate a constructed GP, DGP, or linked (D)GP emulator

This function calculates Leave-One-Out (LOO) cross validation or
Out-Of-Sample (OOS) validation statistics for a constructed GP, DGP, or
linked (D)GP emulator.

## Usage

``` r
validate(
  object,
  x_test,
  y_test,
  method,
  sample_size,
  verb,
  M,
  force,
  cores,
  ...
)

# S3 method for class 'gp'
validate(
  object,
  x_test = NULL,
  y_test = NULL,
  method = "mean_var",
  sample_size = 50,
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)

# S3 method for class 'dgp'
validate(
  object,
  x_test = NULL,
  y_test = NULL,
  method = "mean_var",
  sample_size = 50,
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)

# S3 method for class 'lgp'
validate(
  object,
  x_test = NULL,
  y_test = NULL,
  method = "mean_var",
  sample_size = 50,
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)
```

## Arguments

- object:

  can be one of the following:

  - the S3 class `gp`.

  - the S3 class `dgp`.

  - the S3 class `lgp`.

- x_test:

  OOS testing input data:

  - if `object` is an instance of the `gp` or `dgp` class, `x_test` is a
    matrix where each row is a new input location to be used for
    validating the emulator and each column is an input dimension.

  - if `object` is an instance of the `lgp` class, `x_test` must be a
    matrix representing the global input, where each row corresponds to
    a test data point and each column represents a global input
    dimension. The column indices in `x_test` must align with the
    indices specified in the `From_Output` column of the `struc` data
    frame (used in
    [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md)),
    corresponding to rows where the `From_Emulator` column is
    `"Global"`.

  `x_test` must be provided if `object` is an instance of the `lgp`.
  `x_test` must also be provided if `y_test` is provided. Defaults to
  `NULL`, in which case LOO validation is performed.

- y_test:

  the OOS output data corresponding to `x_test`:

  - if `object` is an instance of the `gp` class, `y_test` is a matrix
    with only one column where each row represents the output
    corresponding to the matching row of `x_test`.

  - if `object` is an instance of the `dgp` class, `y_test` is a matrix
    where each row represents the output corresponding to the matching
    row of `x_test` and with columns representing output dimensions.

  - if `object` is an instance of the `lgp` class, `y_test` can be a
    single matrix or a list of matrices:

    - if `y_test` is a single matrix, then there should be only one
      emulator in the final layer of the linked emulator system and
      `y_test` represents the emulator's output with rows being testing
      positions and columns being output dimensions.

    - if `y_test` is a list, then `y_test` should have *L* matrices,
      where *L* is the number of emulators in the final layer of the
      system. Each matrix has its rows corresponding to testing
      positions and columns corresponding to output dimensions of the
      associated emulator in the final layer.

  `y_test` must be provided if `object` is an instance of the `lgp`.
  `y_test` must also be provided if `x_test` is provided. Defaults to
  `NULL`, in which case LOO validation is performed.

- method:

  **\[updated\]** the prediction approach to use for validation: either
  the mean-variance approach (`"mean_var"`) or the sampling approach
  (`"sampling"`). For details see
  [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md).
  Defaults to `"mean_var"`.

- sample_size:

  the number of samples to draw for each given imputation if
  `method = "sampling"`. Defaults to `50`.

- verb:

  a bool indicating if trace information for validation should be
  printed during function execution. Defaults to `TRUE`.

- M:

  the size of the conditioning set for the Vecchia approximation in
  emulator validation. This argument is only used if the emulator
  `object` was constructed under the Vecchia approximation. Defaults to
  `50`.

- force:

  a bool indicating whether to force LOO or OOS re-evaluation when the
  `loo` or `oos` slot already exists in `object`. When `force = FALSE`,
  `validate()` will only re-evaluate the emulators if the `x_test` and
  `y_test` are not identical to the values in the `oos` slot. If the
  existing `loo` or `oos` validation used a different `M` in a Vecchia
  approximation or a different `method` to the one prescribed in this
  call, the emulator will be re-evaluated. Set `force` to `TRUE` when
  LOO or OOS re-evaluation is required. Defaults to `FALSE`.

- cores:

  the number of processes to be used for validation. If set to `NULL`,
  the number of processes is set to
  `max physical cores available %/% 2`. Defaults to `1`.

- ...:

  N/A.

## Value

- If `object` is an instance of the `gp` class, an updated `object` is
  returned with an additional slot called `loo` (for LOO cross
  validation) or `oos` (for OOS validation) that contains:

  - two slots called `x_train` (or `x_test`) and `y_train` (or `y_test`)
    that contain the validation data points for LOO (or OOS).

  - a column matrix called `mean`, if `method = "mean_var"`, or
    `median`, if `method = "sampling"`, that contains the predictive
    means or medians of the GP emulator at validation positions.

  - three column matrices called `std`, `lower`, and `upper` that
    contain the predictive standard deviations and credible intervals of
    the GP emulator at validation positions. If `method = "mean_var"`,
    the upper and lower bounds of a credible interval are two standard
    deviations above and below the predictive mean. If
    `method = "sampling"`, the upper and lower bounds of a credible
    interval are 2.5th and 97.5th percentiles.

  - a numeric value called `rmse` that contains the root mean/median
    squared error of the GP emulator.

  - a numeric value called `nrmse` that contains the (max-min)
    normalized root mean/median squared error of the GP emulator. The
    max-min normalization uses the maximum and minimum values of the
    validation outputs contained in `y_train` (or `y_test`).

  - an integer called `M` that contains the size of the conditioning set
    used for the Vecchia approximation, if used, for emulator
    validation.

  - an integer called `sample_size` that contains the number of samples
    used for validation if `method = "sampling"`.

  The rows of matrices (`mean`, `median`, `std`, `lower`, and `upper`)
  correspond to the validation positions.

- If `object` is an instance of the `dgp` class, an updated `object` is
  returned with an additional slot called `loo` (for LOO cross
  validation) or `oos` (for OOS validation) that contains:

  - two slots called `x_train` (or `x_test`) and `y_train` (or `y_test`)
    that contain the validation data points for LOO (or OOS).

  - a matrix called `mean`, if `method = "mean_var"`, or `median`, if
    `method = "sampling"`, that contains the predictive means or medians
    of the DGP emulator at validation positions.

  - three matrices called `std`, `lower`, and `upper` that contain the
    predictive standard deviations and credible intervals of the DGP
    emulator at validation positions. If `method = "mean_var"`, the
    upper and lower bounds of a credible interval are two standard
    deviations above and below the predictive mean. If
    `method = "sampling"`, the upper and lower bounds of a credible
    interval are 2.5th and 97.5th percentiles.

  - a vector called `rmse` that contains the root mean/median squared
    errors of the DGP emulator across different output dimensions.

  - a vector called `nrmse` that contains the (max-min) normalized root
    mean/median squared errors of the DGP emulator across different
    output dimensions. The max-min normalization uses the maximum and
    minimum values of the validation outputs contained in `y_train` (or
    `y_test`).

  - an integer called `M` that contains size of the conditioning set
    used for the Vecchia approximation, if used, for emulator
    validation.

  - an integer called `sample_size` that contains the number of samples
    used for validation if `method = "sampling"`.

  The rows and columns of matrices (`mean`, `median`, `std`, `lower`,
  and `upper`) correspond to the validation positions and DGP emulator
  output dimensions, respectively.

- **\[updated\]** If `object` is an instance of the `dgp` class with a
  categorical likelihood, an updated `object` is returned with an
  additional slot called `loo` (for LOO cross validation) or `oos` (for
  OOS validation) that contains:

  - two slots called `x_train` (or `x_test`) and `y_train` (or `y_test`)
    that contain the validation data points for LOO (or OOS).

  - a vector called `label` that contains predictive labels from the DGP
    emulator at validation positions.

  - a matrix called `probability` that contains mean predictive
    probabilities for each class from the DGP emulator at validation
    positions. The matrix has its rows corresponding to validation
    positions and columns corresponding to different classes.

  - a scalar called `log_loss` that represents the log loss of the
    trained DGP classifier. Log loss measures the accuracy of
    probabilistic predictions, with lower values indicating better
    classification performance. `log_loss` ranges from `0` to positive
    infinity, where a value closer to `0` suggests more confident and
    accurate predictions.

  - a scalar called `accuracy` that represents the accuracy of the
    trained DGP classifier. Accuracy measures the proportion of
    correctly classified instances among all predictions, with higher
    values indicating better classification performance. accuracy ranges
    from `0` to `1`, where a value closer to `1` suggests more reliable
    and precise predictions.

  - a slot named `method` indicating whether the matrix in the
    `probability` slot were obtained using the `"mean-var"` method or
    the `"sampling"` method.

  - an integer called `M` that contains size of the conditioning set
    used for the Vecchia approximation, if used, in emulator validation.

  - an integer called `sample_size` that contains the number of samples
    used for validation.

- If `object` is an instance of the `lgp` class, an updated `object` is
  returned with an additional slot called `oos` (for OOS validation)
  that contains:

  - two slots called `x_test` and `y_test` that contain the validation
    data points for OOS.

  - a list called `mean`, if `method = "mean_var"`, or `median`, if
    `method = "sampling"`, that contains the predictive means or medians
    of the linked (D)GP emulator at validation positions.

  - three lists called `std`, `lower`, and `upper` that contain the
    predictive standard deviations and credible intervals of the linked
    (D)GP emulator at validation positions. If `method = "mean_var"`,
    the upper and lower bounds of a credible interval are two standard
    deviations above and below the predictive mean. If
    `method = "sampling"`, the upper and lower bounds of a credible
    interval are 2.5th and 97.5th percentiles.

  - a list called `rmse` that contains the root mean/median squared
    errors of the linked (D)GP emulator.

  - a list called `nrmse` that contains the (max-min) normalized root
    mean/median squared errors of the linked (D)GP emulator. The max-min
    normalization uses the maximum and minimum values of the validation
    outputs contained in `y_test`.

  - an integer called `M` that contains size of the conditioning set
    used for the Vecchia approximation, if used, in emulator validation.

  - an integer called `sample_size` that contains the number of samples
    used for validation if `method = "sampling"`.

  Each element in `mean`, `median`, `std`, `lower`, `upper`, `rmse`, and
  `nrmse` corresponds to a (D)GP emulator in the final layer of the
  linked (D)GP emulator.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- When both `x_test` and `y_test` are `NULL`, LOO cross validation will
  be implemented. Otherwise, OOS validation will be implemented. LOO
  validation is only applicable to a GP or DGP emulator (i.e., `object`
  is an instance of the `gp` or `dgp` class). If a linked (D)GP emulator
  (i.e., `object` is an instance of the `lgp` class) is provided,
  `x_test` and `y_test` must also be provided for OOS validation.

## Examples

``` r
if (FALSE) { # \dontrun{

# See gp(), dgp(), or lgp() for an example.
} # }
```
