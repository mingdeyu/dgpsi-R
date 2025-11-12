# Unpack a bundle of (D)GP emulators

This function unpacks a bundle of (D)GP emulators safely so that any
further manipulations of unpacked individual emulators will not impact
those in the bundle.

## Usage

``` r
unpack(object)
```

## Arguments

- object:

  an instance of the class `bundle`.

## Value

A named list that contains individual emulators (named
`emulator1,...,emulatorS`) packed in `object`, where `S` is the number
of emulators in `object`.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Examples

``` r
if (FALSE) { # \dontrun{

# See pack() for an example.
} # }
```
