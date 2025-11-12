# Serialize the constructed emulator

This function serializes the constructed emulator.

## Usage

``` r
serialize(object, light = TRUE)
```

## Arguments

- object:

  an instance of the S3 class `gp`, `dgp`, `lgp`, or `bundle`.

- light:

  a bool indicating if a light version of the constructed emulator (that
  requires a small storage) will be serialized. Defaults to `TRUE`.

## Value

A serialized version of `object`.

## Details

See further examples and tutorials at
<https://mingdeyu.github.io/dgpsi-R/dev/>.

## Note

Since the constructed emulators are 'python' objects, they cannot be
directly exported to other R processes for parallel processing. This
function provides a solution by converting the emulators into serialized
objects, which can be restored using
[`deserialize()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/deserialize.md)
for multi-process parallel implementation.

## Examples

``` r
if (FALSE) { # \dontrun{

library(parallel)
library(dgpsi)

# model
f <- function(x) {
 (sin(7.5*x)+1)/2
}

# training data
X <- seq(0, 1, length = 10)
Y <- sapply(X, f)

# train a DGP emulator
m <- dgp(X, Y, name = "matern2.5")

# testing input data
X_dgp <- seq(0, 1, length = 100)

# serialize the DGP emulator
m_serialized <- serialize(m)

# create a cluster with 8 workers for parallel predictions
cl <- makeCluster(8)

# export objects to the cluster
clusterExport(cl, varlist = c("m_serialized", "X_dgp"))

# initialize deserialized object on each worker
res <- clusterEvalQ(cl, {
  library(dgpsi)
  assign("m_deserialized", deserialize(m_serialized), envir = .GlobalEnv)
})

# perform parallel predictions
results <- parLapply(cl, 1:length(X_dgp), function(i) {
  mean_i <- predict(m_deserialized, X_dgp[i])$results$mean
})

# reset the cluster
stopCluster(cl)

# combine mean predictions
pred_mean <- do.call(rbind, results)
} # }
```
