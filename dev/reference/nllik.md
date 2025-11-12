# Calculate the predictive negative log-likelihood

This function computes the predictive negative log-likelihood from a DGP
emulator with a likelihood layer.

## Usage

``` r
nllik(object, x, y)
```

## Arguments

- object:

  an instance of the `dgp` class and it should be produced by
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) with
  `likelihood` not being `NULL`;

- x:

  a matrix where each row is an input testing data point and each column
  is an input dimension.

- y:

  a matrix with only one column where each row is a scalar-valued
  testing output data point.

## Value

An updated `object` is returned with an additional slot named `NLL` that
contains two elements. The first one, named `meanNLL`, is a scalar that
gives the average negative predicted log-likelihood across all testing
data points. The second one, named `allNLL`, is a vector that gives the
negative predicted log-likelihood for each testing data point.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.
