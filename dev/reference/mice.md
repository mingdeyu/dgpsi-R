# Locate the next design point for a (D)GP emulator or a bundle of (D)GP emulators using MICE

This function searches from a candidate set to locate the next design
point(s) to be added to a (D)GP emulator or a bundle of (D)GP emulators
using the Mutual Information for Computer Experiments (MICE), see the
reference below.

## Usage

``` r
mice(object, ...)

# S3 method for class 'gp'
mice(
  object,
  x_cand = NULL,
  n_cand = 200,
  batch_size = 1,
  M = 50,
  nugget_s = 1e-06,
  workers = 1,
  limits = NULL,
  int = FALSE,
  ...
)

# S3 method for class 'dgp'
mice(
  object,
  x_cand = NULL,
  n_cand = 200,
  batch_size = 1,
  M = 50,
  nugget_s = 1e-06,
  workers = 1,
  limits = NULL,
  int = FALSE,
  aggregate = NULL,
  ...
)

# S3 method for class 'bundle'
mice(
  object,
  x_cand = NULL,
  n_cand = 200,
  batch_size = 1,
  M = 50,
  nugget_s = 1e-06,
  workers = 1,
  limits = NULL,
  int = FALSE,
  aggregate = NULL,
  ...
)
```

## Arguments

- object:

  can be one of the following:

  - the S3 class `gp`.

  - the S3 class `dgp`.

  - the S3 class `bundle`.

- ...:

  any arguments (with names different from those of arguments used in
  `mice()`) that are used by `aggregate` can be passed here.

- x_cand:

  a matrix (with each row being a design point and column being an input
  dimension) that gives a candidate set from which the next design
  point(s) are determined. If `object` is an instance of the `bundle`
  class and `aggregate` is not supplied, `x_cand` can also be a list.
  The list must have a length equal to the number of emulators in
  `object`, with each element being a matrix representing the candidate
  set for a corresponding emulator in the bundle. Defaults to `NULL`.

- n_cand:

  an integer specifying the size of the candidate set to be generated
  for selecting the next design point(s). This argument is used only
  when `x_cand` is `NULL`. Defaults to `200`.

- batch_size:

  an integer that gives the number of design points to be chosen.
  Defaults to `1`.

- M:

  the size of the conditioning set for the Vecchia approximation in the
  criterion calculation. This argument is only used if the emulator
  `object` was constructed under the Vecchia approximation. Defaults to
  `50`.

- nugget_s:

  the value of the smoothing nugget term used by MICE. Defaults to
  `1e-6`.

- workers:

  the number of processes to be used for the criterion calculation. If
  set to `NULL`, the number of processes is set to
  `max physical cores available %/% 2`. Defaults to `1`.

- limits:

  a two-column matrix that gives the ranges of each input dimension, or
  a vector of length two if there is only one input dimension. If a
  vector is provided, it will be converted to a two-column row matrix.
  The rows of the matrix correspond to input dimensions, and its first
  and second columns correspond to the minimum and maximum values of the
  input dimensions. This argument is only used when `x_cand = NULL`.
  Defaults to `NULL`.

- int:

  a bool or a vector of bools that indicates if an input dimension is an
  integer type. If a single bool is given, it will be applied to all
  input dimensions. If a vector is provided, it should have a length
  equal to the input dimensions and will be applied to individual input
  dimensions. This argument is only used when `x_cand = NULL`. Defaults
  to `FALSE`.

- aggregate:

  an R function that aggregates scores of the MICE across different
  output dimensions (if `object` is an instance of the `dgp` class) or
  across different emulators (if `object` is an instance of the `bundle`
  class). The function should be specified in the following basic form:

  - the first argument is a matrix representing scores. The rows of the
    matrix correspond to different design points. The number of columns
    of the matrix equals to:

    - the emulator output dimension if `object` is an instance of the
      `dgp` class; or

    - the number of emulators contained in `object` if `object` is an
      instance of the `bundle` class.

  - the output should be a vector that gives aggregate scores at
    different design points.

  Set to `NULL` to disable aggregation. Defaults to `NULL`.

## Value

1.  If `x_cand` is not `NULL`:

    - When `object` is an instance of the `gp` class, a vector of length
      `batch_size` is returned, containing the positions (row numbers)
      of the next design points from `x_cand`.

    - When `object` is an instance of the `dgp` class, a vector of
      length `batch_size * D` is returned, containing the positions (row
      numbers) of the next design points from `x_cand` to be added to
      the DGP emulator.

      - `D` is the number of output dimensions of the DGP emulator if no
        likelihood layer is included.

      - For a DGP emulator with a `Hetero` or `NegBin` likelihood layer,
        `D = 2`.

      - For a DGP emulator with a `Categorical` likelihood layer,
        `D = 1` for binary output or `D = K` for multi-class output with
        `K` classes.

    - When `object` is an instance of the `bundle` class, a matrix is
      returned with `batch_size` rows and a column for each emulator in
      the bundle, containing the positions (row numbers) of the next
      design points from `x_cand` for individual emulators.

2.  If `x_cand` is `NULL`:

    - When `object` is an instance of the `gp` class, a matrix with
      `batch_size` rows is returned, giving the next design points to be
      evaluated.

    - When `object` is an instance of the `dgp` class, a matrix with
      `batch_size * D` rows is returned, where:

      - `D` is the number of output dimensions of the DGP emulator if no
        likelihood layer is included.

      - For a DGP emulator with a `Hetero` or `NegBin` likelihood layer,
        `D = 2`.

      - For a DGP emulator with a `Categorical` likelihood layer,
        `D = 1` for binary output or `D = K` for multi-class output with
        `K` classes.

    - When `object` is an instance of the `bundle` class, a list is
      returned with a length equal to the number of emulators in the
      bundle. Each element of the list is a matrix with `batch_size`
      rows, where each row represents a design point to be added to the
      corresponding emulator.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

The first column of the matrix supplied to the first argument of
`aggregate` must correspond to the first output dimension of the DGP
emulator if `object` is an instance of the `dgp` class, and so on for
subsequent columns and dimensions. If `object` is an instance of the
`bundle` class, the first column must correspond to the first emulator
in the bundle, and so on for subsequent columns and emulators.

## References

Beck, J., & Guillas, S. (2016). Sequential design with mutual
information for computer experiments (MICE): emulation of a tsunami
model. *SIAM/ASA Journal on Uncertainty Quantification*, **4(1)**,
739-766.

## Examples

``` r
if (FALSE) { # \dontrun{

# load packages and the Python env
library(lhs)
library(dgpsi)

# construct a 1D non-stationary function
f <- function(x) {
 sin(30*((2*x-1)/2-0.4)^5)*cos(20*((2*x-1)/2-0.4))
}

# generate the initial design
X <- maximinLHS(10,1)
Y <- f(X)

# training a 2-layered DGP emulator with the global connection off
m <- dgp(X, Y, connect = F)

# generate a candidate set
x_cand <- maximinLHS(200,1)

# locate the next design point using MICE
next_point <- mice(m, x_cand = x_cand)
X_new <- x_cand[next_point,,drop = F]

# obtain the corresponding output at the located design point
Y_new <- f(X_new)

# combine the new input-output pair to the existing data
X <- rbind(X, X_new)
Y <- rbind(Y, Y_new)

# update the DGP emulator with the new input and output data and refit
m <- update(m, X, Y, refit = TRUE)

# plot the LOO validation
plot(m)
} # }
```
