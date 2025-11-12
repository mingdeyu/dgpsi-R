# Pack GP and DGP emulators into a bundle

This function packs GP emulators and DGP emulators into a `bundle` class
for sequential designs if each emulator emulates one output dimension of
the underlying simulator.

## Usage

``` r
pack(..., id = NULL)
```

## Arguments

- ...:

  a sequence or a list of emulators produced by
  [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md) or
  [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md).

- id:

  an ID to be assigned to the bundle emulator. If an ID is not provided
  (i.e., `id = NULL`), a UUID (Universally Unique Identifier) will be
  automatically generated and assigned to the emulator. Default to
  `NULL`.

## Value

An S3 class named `bundle` to be used by
[`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
for sequential designs. It has:

- a slot called `id` that is assigned through the `id` argument.

- *N* slots named `emulator1,...,emulatorN`, each of which contains a GP
  or DGP emulator, where *N* is the number of emulators that are
  provided to the function.

- a slot called `data` which contains two elements `X` and `Y`. `X`
  contains *N* matrices named `emulator1,...,emulatorN` that are
  training input data for different emulators. `Y` contains *N*
  single-column matrices named `emulator1,...,emulatorN` that are
  training output data for different emulators.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Examples

``` r
if (FALSE) { # \dontrun{

# load packages
library(lhs)
library(dgpsi)

# construct a function with a two-dimensional output
f <- function(x) {
 y1 = sin(30*((2*x-1)/2-0.4)^5)*cos(20*((2*x-1)/2-0.4))
 y2 = 1/3*sin(2*(2*x - 1))+2/3*exp(-30*(2*(2*x-1))^2)+1/3
 return(cbind(y1,y2))
}

# generate the initial design
X <- maximinLHS(10,1)
Y <- f(X)

# generate the validation data
validate_x <- maximinLHS(30,1)
validate_y <- f(validate_x)

# training a 2-layered DGP emulator with respect to each output with the global connection off
m1 <- dgp(X, Y[,1], connect = F)
m2 <- dgp(X, Y[,2], connect = F)

# specify the range of the input dimension
lim <- c(0, 1)

# pack emulators to form an emulator bundle
m <- pack(m1, m2)

# 1st wave of the sequential design with 10 iterations and the target RMSE of 0.01
m <- design(m, N = 10, limits = lim, f = f, x_test = validate_x, y_test = validate_y, target = 0.01)

# 2rd wave of the sequential design with additional 10 iterations and the same target
m <- design(m, N = 10, limits = lim, f = f, x_test = validate_x, y_test = validate_y, target = 0.01)

# draw sequential designs of the two packed emulators
draw(m, type = 'design')

# inspect the traces of RMSEs of the two packed emulators
draw(m, type = 'rmse')

# write and read the constructed emulator bundle
write(m, 'bundle_dgp')
m <- read('bundle_dgp')

# unpack the bundle into individual emulators
m_unpacked <- unpack(m)

# plot OOS validations of individual emulators
plot(m_unpacked[[1]], x_test = validate_x, y_test = validate_y[,1])
plot(m_unpacked[[2]], x_test = validate_x, y_test = validate_y[,2])
} # }
```
