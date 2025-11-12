# Summary of a constructed GP, DGP, or linked (D)GP emulator

This function provides a summary of key information for a GP, DGP, or
linked (D)GP emulator by generating either a table or an interactive
plot of the emulator’s structure.

## Usage

``` r
# S3 method for class 'gp'
summary(object, type = "plot", ...)

# S3 method for class 'dgp'
summary(object, type = "plot", ...)

# S3 method for class 'lgp'
summary(object, type = "plot", group_size = 1, ...)
```

## Arguments

- object:

  can be one of the following:

  - the S3 class `gp`.

  - the S3 class `dgp`.

  - the S3 class `lgp`.

- type:

  a character string, either `"table"` or `"plot"`, indicating the
  format of the output. If set to `"table"`, the function returns a
  summary in table. If set to `"plot"`, the function returns an
  interactive visualization. Defaults to `"plot"`.

- ...:

  Any arguments that can be passed to
  [`kableExtra::kbl()`](https://rdrr.io/pkg/kableExtra/man/kbl.html)
  when `type = "table"`.

- group_size:

  an integer specifying the number of consecutive layers to be grouped
  together in the interactive visualization of linked emulators when
  `type = "plot"`. This argument is only applicable if `object` is an
  instance of the `lgp` class. Defaults to `1`.

## Value

Either a summary table (returned as `kableExtra` object) or an
interactive visualization (returned as a `visNetwork` object) of the
emulator. The visualization is compatible with R Markdown documents and
the RStudio Viewer. The summary table can be further customized by
[kableExtra::kableExtra](https://rdrr.io/pkg/kableExtra/man/kableExtra-package.html)
package. The resulting `visNetwork` object can be saved as an HTML file
using
[`visNetwork::visSave()`](https://rdrr.io/pkg/visNetwork/man/visSave.html)
from the
[visNetwork::visNetwork](https://rdrr.io/pkg/visNetwork/man/visNetwork.html)
package.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Examples

``` r
if (FALSE) { # \dontrun{

# See gp(), dgp(), or lgp() for an example.
} # }
```
