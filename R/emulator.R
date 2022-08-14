#' @title Initialize the emulator object
#'
#' @description This function constructs an object for deep Gaussian process predictions.
#'
#' @param all_layer the trained DGP model produced by [estimate()].
#' @param N the number of imputation to produce the predictions. Increase the value to account for
#'     more imputation uncertainties. Defaults to `50`.
#' @param nb_parallel whether to use the multi-threading to accelerate the predictions. Defaults to `FALSE`.
#'
#' @return An emulator object to be used by [predict()] or [ppredict()] for DGP predictions.
#' @details See examples in tutorials at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
emulator <- function(all_layer, N = 50, nb_parallel = FALSE) {
  N <- as.integer(N)
  res <- pkg.env$dgpsi$emulator(all_layer, N, nb_parallel)
  return(res)
}
