# Save the constructed emulator

This function saves the constructed emulator to a `.pkl` file.

## Usage

``` r
write(object, pkl_file, light = TRUE)
```

## Arguments

- object:

  an instance of the S3 class `gp`, `dgp`, `lgp`, or `bundle`.

- pkl_file:

  the path to and the name of the `.pkl` file to which the emulator
  `object` is saved.

- light:

  a bool indicating if a light version of the constructed emulator (that
  requires less disk space to store) will be saved. Defaults to `TRUE`.

## Value

No return value. `object` will be saved to a local `.pkl` file specified
by `pkl_file`.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

Since emulators built from the package are 'python' objects,
[`save()`](https://rdrr.io/r/base/save.html) from R will not work as it
would for R objects. If `object` was processed by
[`set_vecchia()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_vecchia.md)
to add or remove the Vecchia approximation, `light` should be set to
`FALSE` to ensure reproducibility after the saved emulator is reloaded
by [`read()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/read.md).

## Examples

``` r
if (FALSE) { # \dontrun{

# See gp(), dgp(), lgp(), or pack() for an example.
} # }
```
