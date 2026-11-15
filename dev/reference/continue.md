# Continue training a DGP emulator

This function implements additional training iterations for a DGP
emulator.

## Usage

``` r
continue(object, N = NULL, ess_burn = 10, verb = TRUE, burnin = NULL, B = NULL)
```

## Arguments

- object:

  an instance of the `dgp` class.

- N:

  additional number of iterations to train the DGP emulator. If set to
  `NULL`, the number of iterations is set to `500` if the DGP emulator
  was constructed without the Vecchia approximation, and is set to `200`
  if Vecchia approximation was used. Defaults to `NULL`.

- ess_burn:

  number of burnin steps for ESS-within-Gibbs at each I-step of the
  training. Defaults to `10`.

- verb:

  a bool indicating if a progress bar will be printed during training.
  Defaults to `TRUE`.

- burnin:

  the number of training iterations to be discarded for point estimates
  calculation. Must be smaller than the overall training iterations
  so-far implemented. If this is not specified, only the last 25% of
  iterations are used. This overrides the value of `burnin` set in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).
  Defaults to `NULL`.

- B:

  the number of imputations to produce predictions. Increase the value
  to account for more imputation uncertainty. This overrides the value
  of `B` set in
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) if
  `B` is not `NULL`. Defaults to `NULL`.

## Value

An updated `object`.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

- One can also use this function to fit an untrained DGP emulator
  constructed by
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) with
  `training = FALSE`.

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
