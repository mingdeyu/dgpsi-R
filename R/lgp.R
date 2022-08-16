#' @title Initialize the linked (D)GP object
#'
#' @description This function constructs an object for linked (D)GP predictions.
#'
#' @param all_layer a list contains *L* (the number of layers of a systems of computer models) sub-lists,
#'     each of which represents a layer and contains the GP/DGP emulators of computer models represented by
#'     containers constructed by [container()]. The sub-lists are placed in the list in the same order of
#'     the specified computer model system.
#' @param N the number of imputation to produce the predictions. Increase the value to account for more
#'     imputation uncertainties. If the system consists only GP emulators, `N` is set to `1` automatically.
#'     Defaults to `50`.
#' @param nb_parallel whether to use the multi-threading to accelerate the predictions. Defaults to `FALSE`.
#'
#' @return A linked (D)GP object to be used by [predict()] or [ppredict()] for linked (D)GP predictions.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
lgp <- function(all_layer, N = 50, nb_parallel = FALSE) {
  N <- as.integer(N)
  res <- pkg.env$dgpsi$lgp(all_layer, N, nb_parallel)
  return(res)
}
