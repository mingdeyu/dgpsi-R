#' @title Linked (D)GP emulator construction
#'
#' @description This function constructs a linked (D)GP emulator.
#'
#' @param struc a list contains *L* (the number of layers in a systems of computer models) sub-lists,
#'     each of which represents a layer and contains (D)GP emulators (represented by
#'     instances of S3 class `gp` or `dgp`) of computer models. The sub-lists are placed in the list
#'     in the same order of the specified computer model system's hierarchy.
#' @param B the number of imputations to produce the predictions. Increase the value to account for more
#'     imputation uncertainties. Decrease the value for lower imputation uncertainties but faster predictions.
#'     If the system consists only GP emulators, `B` is set to `1` automatically. Defaults to `50`.
#'
#' @return An S3 class named `lgp` to be used by [predict()] for linked (D)GP predictions.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
lgp <- function(struc, B = 50) {
  B <- as.integer(B)
  L <- length(struc)
  extracted_struc <- list()
  for ( l in 1:L ) {
    layer <- list()
    K <- length(struc[[l]])
    for (k in 1:K) {
      cont <- struc[[l]][[k]]$container_obj
      if ( is.null(cont$local_input_idx) ){
        stop(sprintf("Emulator %i in Layer %i has no 'linked_idx' specified. Use set_linked_idx() to specify this attribute.", k, l), call. = FALSE)
      }
      layer[[k]] <- cont
    }
    extracted_struc[[l]] <- layer
  }
  obj <- pkg.env$dgpsi$lgp(all_layer = extracted_struc, N = B)
  res <- list(emulator_obj = obj)
  class(res) <- "lgp"
  return(res)
}
