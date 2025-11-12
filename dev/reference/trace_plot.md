# Trace plot for DGP hyperparameters

This function draws trace plots for the hyperparameters of a chosen GP
node in a DGP emulator.

## Usage

``` r
trace_plot(object, layer = NULL, node = 1)
```

## Arguments

- object:

  an instance of the `dgp` class.

- layer:

  the index of a layer. Defaults to `NULL` for the final layer.

- node:

  the index of a GP node in the layer specified by `layer`. Defaults to
  `1` for the first GP node in the corresponding layer.

## Value

A `ggplot` object.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Examples

``` r
if (FALSE) { # \dontrun{

# See dgp() for an example.
} # }
```
