# Prediction from GP, DGP, or linked (D)GP emulators

This function implements prediction from GP, DGP, or linked (D)GP
emulators.

## Usage

``` r
# S3 method for class 'dgp'
predict(
  object,
  x,
  method = "mean_var",
  full_layer = FALSE,
  sample_size = 50,
  M = 50,
  cores = 1,
  chunks = NULL,
  ...
)

# S3 method for class 'lgp'
predict(
  object,
  x,
  method = "mean_var",
  full_layer = FALSE,
  sample_size = 50,
  M = 50,
  cores = 1,
  chunks = NULL,
  ...
)

# S3 method for class 'gp'
predict(
  object,
  x,
  method = "mean_var",
  sample_size = 50,
  M = 50,
  cores = 1,
  chunks = NULL,
  ...
)
```

## Arguments

- object:

  an instance of the `gp`, `dgp`, or `lgp` class.

- x:

  the testing input data:

  - if `object` is an instance of the `gp` or `dgp` class, `x` is a
    matrix where each row is an input testing data point and each column
    is an input dimension.

  - if `object` is an instance of the `lgp` class, `x` must be a matrix
    representing the global input, where each row corresponds to a test
    data point and each column represents a global input dimension. The
    column indices in `x` must align with the indices specified in the
    `From_Output` column of the `struc` data frame (used in
    [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md)),
    corresponding to rows where the `From_Emulator` column is
    `"Global"`.

- method:

  **\[updated\]** the prediction approach to use: either the
  mean-variance approach (`"mean_var"`) or the sampling approach
  (`"sampling"`). The mean-variance approach returns the means and
  variances of the predictive distributions, while the sampling approach
  generates samples from predictive distributions using the derived
  means and variances. Defaults to `"mean_var"`.

- full_layer:

  a bool indicating whether to output the predictions of all layers.
  Defaults to `FALSE`. Only used when `object` is a DGP or a linked
  (D)GP emulator.

- sample_size:

  the number of samples to draw for each given imputation if
  `method = "sampling"`. Defaults to `50`.

- M:

  the size of the conditioning set for the Vecchia approximation in the
  emulator prediction. Defaults to `50`. This argument is only used if
  the emulator `object` was constructed under the Vecchia approximation.

- cores:

  the number of processes to be used for prediction. If set to `NULL`,
  the number of processes is set to
  `max physical cores available %/% 2`. Defaults to `1`.

- chunks:

  the number of chunks that the testing input matrix `x` will be divided
  into for multi-cores to work on. Only used when `cores` is not `1`. If
  not specified (i.e., `chunks = NULL`), the number of chunks is set to
  the value of `cores`. Defaults to `NULL`.

- ...:

  N/A.

## Value

- If `object` is an instance of the `gp` class:

  1.  if `method = "mean_var"`: an updated `object` is returned with an
      additional slot called `results` that contains two matrices named
      `mean` for the predictive means and `var` for the predictive
      variances. Each matrix has only one column with its rows
      corresponding to testing positions (i.e., rows of `x`).

  2.  if `method = "sampling"`: an updated `object` is returned with an
      additional slot called `results` that contains a matrix whose rows
      correspond to testing positions and columns correspond to
      `sample_size` number of samples drawn from the predictive
      distribution of GP.

- **\[updated\]** If `object` is an instance of the `dgp` class:

  1.  if `method = "mean_var"` and `full_layer = FALSE`: an updated
      `object` is returned with an additional slot called `results` that
      contains two matrices named `mean` for the predictive means and
      `var` for the predictive variances respectively. Each matrix has
      its rows corresponding to testing positions and columns
      corresponding to DGP global output dimensions (i.e., the number of
      GP/likelihood nodes in the final layer). If the likelihood node is
      categorical, the matrices contain the predictive means and
      variances of the class probabilities, with columns corresponding
      to different classes.

  2.  if `method = "mean_var"` and `full_layer = TRUE`: an updated
      `object` is returned with an additional slot called `results` that
      contains two sub-lists named `mean` for the predictive means and
      `var` for the predictive variances respectively. Each sub-list
      contains *L* (i.e., the number of layers) matrices named
      `layer1, layer2,..., layerL`. Each matrix has its rows
      corresponding to testing positions and columns corresponding to
      output dimensions (i.e., the number of GP/likelihood nodes from
      the associated layer). If the likelihood node is categorical, the
      matrices named `layerL` in both `mean` and `var` contain the
      predictive means and variances of the class probabilities,
      respectively, with columns corresponding to different classes.

  3.  if `method = "sampling"` and `full_layer = FALSE`: an updated
      `object` is returned with an additional slot called `results` that
      contains *D* (i.e., the number of GP/likelihood nodes in the final
      layer) matrices named `output1, output2,..., outputD`. If the
      likelihood node in the final layer is categorical, `results`
      contains *D* matrices (where *D* is the number of classes) of
      sampled class probabilities, each named according to its
      corresponding class label. Each matrix in `results` has its rows
      corresponding to testing positions and columns corresponding to
      samples of size: `B * sample_size`, where `B` is the number of
      imputations specified in
      [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).

  4.  if `method = "sampling"` and `full_layer = TRUE`: an updated
      `object` is returned with an additional slot called `results` that
      contains *L* (i.e., the number of layers) sub-lists named
      `layer1, layer2,..., layerL`. Each sub-list represents samples
      drawn from the GP/likelihood nodes in the corresponding layer, and
      contains *D* (i.e., the number of GP/likelihood nodes in the
      corresponding layer) matrices named
      `output1, output2,..., outputD`. If the likelihood node in the
      final layer is categorical, `layerL` contains *D* matrices (where
      *D* is the number of classes) of sampled class probabilities, each
      named according to its corresponding class label. Each matrix has
      its rows corresponding to testing positions and columns
      corresponding to samples of size: `B * sample_size`, where `B` is
      the number of imputations specified in
      [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).

- If `object` is an instance of the `lgp` class:

  1.  if `method = "mean_var"` and `full_layer = FALSE`: an updated
      `object` is returned with an additional slot called `results` that
      contains two sub-lists named `mean` for the predictive means and
      `var` for the predictive variances respectively. Each sub-list
      contains *K* (same number of emulators in the final layer of the
      system) matrices named using the `ID`s of the corresponding
      emulators in the final layer. Each matrix has rows corresponding
      to global testing positions and columns corresponding to output
      dimensions of the associated emulator in the final layer.

  2.  if `method = "mean_var"` and `full_layer = TRUE`: an updated
      `object` is returned with an additional slot called `results` that
      contains two sub-lists named `mean` for the predictive means and
      `var` for the predictive variances respectively. Each sub-list
      contains *L* (i.e., the number of layers in the emulated system)
      components named `layer1, layer2,..., layerL`. Each component
      represents a layer and contains *K* (same number of emulators in
      the corresponding layer of the system) matrices named using the
      `ID`s of the corresponding emulators in that layer. Each matrix
      has its rows corresponding to global testing positions and columns
      corresponding to output dimensions of the associated GP/DGP
      emulator in the corresponding layer.

  3.  if `method = "sampling"` and `full_layer = FALSE`: an updated
      `object` is returned with an additional slot called `results` that
      contains *K* (same number of emulators in the final layer of the
      system) sub-lists named using the `ID`s of the corresponding
      emulators in the final layer. Each sub-list contains *D* matrices,
      named `output1, output2,..., outputD`, that correspond to the
      output dimensions of the GP/DGP emulator. Each matrix has rows
      corresponding to testing positions and columns corresponding to
      samples of size: `B * sample_size`, where `B` is the number of
      imputations specified in
      [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md).

  4.  if `method = "sampling"` and `full_layer = TRUE`: an updated
      `object` is returned with an additional slot called `results` that
      contains *L* (i.e., the number of layers of the emulated system)
      sub-lists named `layer1, layer2,..., layerL`. Each sub-list
      represents a layer and contains *K* (same number of emulators in
      the corresponding layer of the system) components named using the
      `ID`s of the corresponding emulators in that layer. Each component
      contains *D* matrices, named `output1, output2,..., outputD`, that
      correspond to the output dimensions of the GP/DGP emulator. Each
      matrix has its rows corresponding to testing positions and columns
      corresponding to samples of size: `B * sample_size`, where `B` is
      the number of imputations specified in
      [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md).

The `results` slot will also include:

- the value of `M`, which represents the size of the conditioning set
  for the Vecchia approximation, if used, in the emulator prediction.

- the value of `sample_size` if `method = "sampling"`.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Examples

``` r
if (FALSE) { # \dontrun{

# See gp(), dgp(), or lgp() for an example.
} # }
```
