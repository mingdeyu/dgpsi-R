# Reset number of imputations for a DGP emulator

This function resets the number of imputations for prediction from a DGP
emulator.

## Usage

``` r
set_imp(object, B = 5)
```

## Arguments

- object:

  an instance of the S3 class `dgp`.

- B:

  the number of imputations to produce predictions from `object`.
  Increase the value to improve imputation uncertainty quantification.
  Decrease the value to improve speed of prediction. Defaults to `5`.

## Value

An updated `object` with the information of `B` incorporated.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- This function is useful when a DGP emulator has been trained and one
  wants to make faster predictions by decreasing the number of
  imputations without rebuilding the emulator.

- The following slots:

  - `loo` and `oos` created by
    [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md);
    and

  - `results` created by
    [`predict()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
    in `object` will be removed and not contained in the returned
    object.

## Examples

``` r
if (FALSE) { # \dontrun{

# See design() for an example.
} # }
```
