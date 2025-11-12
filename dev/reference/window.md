# Trim the sequence of hyperparameter estimates within a DGP emulator

This function trims the sequence of hyperparameter estimates within a
DGP emulator generated during training.

## Usage

``` r
window(object, start, end = NULL, thin = 1)
```

## Arguments

- object:

  an instance of the S3 class `dgp`.

- start:

  the first iteration before which all iterations are trimmed from the
  sequence.

- end:

  the last iteration after which all iterations are trimmed from the
  sequence. Set to `NULL` to keep all iterations after (including)
  `start`. Defaults to `NULL`.

- thin:

  the interval between the `start` and `end` iterations to thin out the
  sequence. Defaults to 1.

## Value

An updated `object` with a trimmed sequence of hyperparameters.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- This function is useful when a DGP emulator has been trained and one
  wants to trim the sequence of hyperparameters estimated and to use the
  trimmed sequence to generate point estimates of the DGP model
  parameters for prediction.

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

# See dgp() for an example.
} # }
```
