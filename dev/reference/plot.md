# Validation plots of a constructed GP, DGP, or linked (D)GP emulator

This function draws validation plots of a GP, DGP, or linked (D)GP
emulator.

## Usage

``` r
# S3 method for class 'dgp'
plot(
  x,
  x_test = NULL,
  y_test = NULL,
  dim = NULL,
  method = "mean_var",
  sample_size = 50,
  style = 1,
  min_max = TRUE,
  normalize = TRUE,
  color = "turbo",
  type = "points",
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)

# S3 method for class 'lgp'
plot(
  x,
  x_test = NULL,
  y_test = NULL,
  dim = NULL,
  method = "mean_var",
  sample_size = 50,
  style = 1,
  min_max = TRUE,
  color = "turbo",
  type = "points",
  M = 50,
  verb = TRUE,
  force = FALSE,
  cores = 1,
  ...
)

# S3 method for class 'gp'
plot(
  x,
  x_test = NULL,
  y_test = NULL,
  dim = NULL,
  method = "mean_var",
  sample_size = 50,
  style = 1,
  min_max = TRUE,
  color = "turbo",
  type = "points",
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)
```

## Arguments

- x:

  can be one of the following emulator classes:

  - the S3 class `gp`.

  - the S3 class `dgp`.

  - the S3 class `lgp`.

- x_test:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- y_test:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- dim:

  if `dim = NULL`, the index of an emulator's input within the design
  will be shown on the x-axis in validation plots. Otherwise, `dim`
  indicates which dimension of an emulator's input will be shown on the
  x-axis in validation plots:

  - If `x` is an instance of the `gp` of `dgp` class, `dim` is an
    integer.

  - If `x` is an instance of the `lgp` class, `dim` is an integer
    referring to the dimension of the global input to the linked
    emulator system.

  This argument is only used when `style = 1`. Defaults to `NULL`.

- method:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- sample_size:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- style:

  either `1` or `2`, indicating two different plotting styles for
  validation.

- min_max:

  a bool indicating if min-max normalization will be used to scale the
  testing output, RMSE, predictive mean and std from the emulator.
  Defaults to `TRUE`. This argument is not applicable to DGP emulators
  with categorical likelihoods.

- normalize:

  **\[new\]** a bool indicating if normalization will be used to scale
  the counts in validation plots of DGP emulators with categorical
  likelihoods when `style = 2`. Defaults to `TRUE`.

- color:

  a character string indicating the color map to use when `style = 2`:

  - `'magma'` (or `'A'`)

  - `'inferno'` (or `'B'`)

  - `'plasma'` (or '`C`')

  - `'viridis'` (or `'D'`)

  - `'cividis'` (or `'E'`)

  - `'rocket'` (or `'F'`)

  - `'mako'` (or `'G'`)

  - `'turbo'` (or `'H'`)

  Defaults to `'turbo'` (or `'H'`).

- type:

  either `'line'` or `'points`, indicating whether to draw testing data
  in the OOS validation plot as a line or individual points when the
  input of the emulator is one-dimensional and `style = 1`. This
  argument is not applicable to DGP emulators with categorical
  likelihoods. Defaults to `'points'`

- verb:

  a bool indicating if trace information on plotting will be printed
  during execution. Defaults to `TRUE`.

- M:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- force:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- cores:

  same as that of
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md).

- ...:

  N/A.

## Value

A `patchwork` object.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- `plot()` calls
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  internally to obtain validation results for plotting. However,
  `plot()` will not export the emulator object with validation results.
  Instead, it only returns the plotting object. For small-scale
  validations (i.e., small training or testing data points), direct
  execution of `plot()` works well. However, for moderate- to
  large-scale validation, it is recommended to first run
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  to obtain and store validation results in the emulator object, and
  then supply the object to `plot()`. `plot()` checks the object's `loo`
  and `oos` slots prior to calling
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  and will not perform further calculation if the required information
  is already stored.

- `plot()` will only use stored OOS validation if `x_test` and `y_test`
  are identical to those used by
  [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  to produce the data contained in the object's `oos` slot, otherwise
  `plot()` will re-evaluate OOS validation before plotting.

- The returned
  [patchwork::patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
  object contains the
  [ggplot2::ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
  objects. One can modify the included individual ggplots by accessing
  them with double-bracket indexing. See
  <https://patchwork.data-imaginist.com/> for further information.

## Examples

``` r
if (FALSE) { # \dontrun{

# See gp(), dgp(), or lgp() for an example.
} # }
```
