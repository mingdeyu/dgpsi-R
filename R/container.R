#' @title Initialize the container object
#'
#' @description This function constructs an object that contains the trained GP or DGP emulator of a computer model
#'     for linked (D)GP emulation.
#'
#' @param structure a list that contains the trained structure of GP or DGP of a computer model. For GP,
#'     this is the list exported from [export()]. For DGP, this is the list exported from [estimate()] .
#' @param local_input_idx a vector that specifies the indices (i.e., the columns) of outputs (a matrix)
#'     produced by all models in the feeding layer that are input to the emulator
#'     represented by the `structure` argument. The indices should be ordered in such a way that the extracted
#'     output from the feeding layer is sorted in the same order as the training input used for the GP/DGP
#'     emulator that the `structure` argument represents. When the emulator is in the first
#'     layer, `local_input_idx` gives the indices of its input in the global testing input set, see [predict()]
#'     for descriptions of the global testing input set. Defaults to `NULL`. When the
#'     argument is `NULL`, one needs to set its value by applying [set_local_input()] to the returned container object.
#'
#' @return A container object to be used for linked GP emulation.
#' @md
#' @export
container <- function(structure, local_input_idx = NULL) {
  if(!is.null(local_input_idx)){
    local_input_idx <- reticulate::np_array(as.integer(local_input_idx - 1))
  }
  res <- pkg.env$dgpsi$container(structure, local_input_idx)
  return(res)
}


#' @title Set local input indices
#'
#' @description This function sets the value of `local_input_idx` to a container object if the value is not set
#'     when the container object is constructed by [container()].
#'
#' @param obj the container object produced by [container()].
#' @param idx same as `local_input_idx` in [container()].
#'
#' @return A container object to be used for linked GP emulation.
#'
#' @details This function is useful when different models are emulated by different teams. Each team can create the container
#'     of their model even without knowing how different models are connected together. When this information is available and
#'     containers of different emulators are collected, the connections between emulators can then be set by assigning
#'     values to `local_input_idx` of each container with this function.
#' @md
#' @export
set_local_input <- function(obj, idx) {
  idx <- reticulate::np_array(as.integer(idx - 1))
  obj$set_local_input(idx)
  return(obj)
}
