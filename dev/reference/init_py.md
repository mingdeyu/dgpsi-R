# 'python' environment initialization

**\[updated\]**

This function initializes the 'python' environment for the package.

## Usage

``` r
init_py(
  py_ver = NULL,
  dgpsi_ver = NULL,
  conda = NULL,
  reinstall = FALSE,
  uninstall = FALSE,
  verb = TRUE,
  show_config = FALSE
)
```

## Arguments

- py_ver:

  a string that gives the 'python' version to be installed. Supported
  versions are 3.10, 3.11, and 3.12. If `py_ver = NULL`, the default
  'python' version '3.10' will be installed.

- dgpsi_ver:

  a string that gives the 'python' version of 'dgpsi' to be used. If
  `dgpsi_ver = NULL`,

  - the latest 'python' version of 'dgpsi' will be used, if the package
    is installed from CRAN;

  - the development 'python' version of 'dgpsi' will be used, if the
    package is installed from GitHub.

- conda:

  optional path to a conda binary. If `NULL`, `init_py()` attempts to
  locate conda automatically; see
  [`reticulate::conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  for the search order. If provided, the path is used directly and must
  point to an existing conda executable.

- reinstall:

  a bool that indicates whether to reinstall the 'python' version of
  'dgpsi' specified in `dgpsi_ver` if it has already been installed.
  This argument is useful when the development version of the R package
  is installed and one may want to regularly update the development
  'python' version of 'dgpsi'. Defaults to `FALSE`.

- uninstall:

  a bool that indicates whether to uninstall the 'python' version of
  'dgpsi' specified in `dgpsi_ver` if it has already been installed.
  This argument is useful when the 'python' environment is corrupted and
  one wants to completely uninstall and reinstall it. Defaults to
  `FALSE`.

- verb:

  a bool indicating if trace information will be printed during function
  execution. Defaults to `TRUE`.

- show_config:

  a bool indicating whether to print the Python backend configuration
  when `verb = TRUE`. Defaults to `FALSE`.

## Value

No return value, called to install required 'python' environment.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

On Linux, `init_py()` may ask for permission during installation. To
auto-accept these prompts in non-interactive settings, set the
environment variable `DGPSI_AUTO_YES` to `"TRUE"` before calling
`init_py()` for the first time. Accepted values include `"TRUE"`,
`"true"`, `"yes"`, `"y"`, and `"1"`. If unset, `init_py()` will continue
to prompt interactively.

## Examples

``` r
if (FALSE) { # \dontrun{

# See gp(), dgp(), or lgp() for an example.
} # }
```
