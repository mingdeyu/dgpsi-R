# Add or remove the Vecchia approximation

This function adds or removes the Vecchia approximation from a GP, DGP
or linked (D)GP emulator constructed by
[`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md),
[`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md) or
[`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md).

## Usage

``` r
set_vecchia(object, vecchia = TRUE, M = 25, ord = NULL)
```

## Arguments

- object:

  an instance of the S3 class `gp`, `dgp`, or `lgp`.

- vecchia:

  a bool to indicate the addition or removal of the Vecchia
  approximation:

  - if `object` is an instance of the `gp` or `dgp` class, `vecchia`
    indicates either addition (`vecchia = TRUE`) or removal
    (`vecchia = FALSE`) of the Vecchia approximation from `object`.

  - if `object` is an instance of the `lgp` class, `vecchia` indicates
    either addition (`vecchia = TRUE`) or removal (`vecchia = FALSE`) of
    the Vecchia approximation from all individual (D)GP emulators
    contained in `object`.

  Defaults to `TRUE`.

- M:

  the size of the conditioning set for the Vecchia approximation in the
  (D)GP emulator training. Defaults to `25`.

- ord:

  an R function that returns the ordering of the input to the (D)GP
  emulator for the Vecchia approximation. The function must satisfy the
  following basic rules:

  - the first argument represents the lengthscale-scaled input to the GP
    emulator or the lengthscale-scaled input to a GP node of the DGP
    emulator.

  - the output of the function is a vector of indices that gives the
    ordering of the input to the GP emulator or the input to the GP
    nodes of the DGP emulator.

  If `ord = NULL`, the default random ordering is used. Defaults to
  `NULL`.

## Value

An updated `object` with the Vecchia approximation either added or
removed.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

This function is useful for quickly switching between Vecchia and
non-Vecchia approximations for an existing emulator without the need to
reconstruct the emulator. If the emulator was built without the Vecchia
approximation, the function can add it, and if the emulator was built
with the Vecchia approximation, the function can remove it. If the
current state already matches the requested state, the emulator remains
unchanged.
