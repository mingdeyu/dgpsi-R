# Validation and diagnostic plots for a sequential design

This function draws diagnostic and validation plots for a sequential
design of a (D)GP emulator or a bundle of (D)GP emulators.

## Usage

``` r
draw(object, ...)

# S3 method for class 'gp'
draw(object, type = "rmse", log = FALSE, ...)

# S3 method for class 'dgp'
draw(object, type = "rmse", log = FALSE, ...)

# S3 method for class 'bundle'
draw(object, type = "rmse", log = FALSE, emulator = NULL, ...)
```

## Arguments

- object:

  can be one of the following emulator classes:

  - the S3 class `gp`.

  - the S3 class `dgp`.

  - the S3 class `bundle`.

- ...:

  N/A.

- type:

  specifies the type of plot or visualization to generate:

  - `"rmse"`: generates a trace plot of RMSEs, log-losses for DGP
    emulators with categorical likelihoods, or custom evaluation metrics
    specified via the `"eval"` argument in the `[design()]` function.

  - `"design"`: shows visualizations of input designs created by the
    sequential design procedure.

  Defaults to `"rmse"`.

- log:

  a bool indicating whether to plot RMSEs, log-losses (for DGP emulators
  with categorical likelihoods), or custom evaluation metrics on a log
  scale when `type = "rmse"`. Defaults to `FALSE`.

- emulator:

  an index or vector of indices of emulators packed in `object`. This
  argument is only used if `object` is an instance of the `bundle`
  class. When set to `NULL`, all emulators in the bundle are drawn.
  Defaults to `NULL`.

## Value

A `patchwork` object.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Examples

``` r
if (FALSE) { # \dontrun{

# See design() for an example.
} # }
```
