#' @title Combine layers
#'
#' @description This function combines layers into one list as a DGP or linked (D)GP structure.
#'
#' @param ... a sequence of lists, each of which contains the GP nodes (produced by [kernel()]),
#'     likelihood nodes (e.g., produced by [Poisson()]), or containers (produced by [container()])
#'     in that layer.
#'
#' @return A list of layers defining the DGP or linked (D)GP structure.
#'
#' @md
#' @export
combine <- function(...) {
  arguments = list(...)
  res <- do.call(pkg.env$dgpsi$combine, arguments)
  return(res)
}


#' @title Save the constructed emulator object
#'
#' @description This function saves the constructed emulator to a `.pkl` file.
#'
#' @param obj an emulator object. For GP, it is the one produced by [gp()] (either before
#'     or after applying [train()]). For DGP, it can be the one produced by [dgp()] (either before
#'     or after applying [train()]), or the one produced by [emulator()]. For linked (D)GP,
#'     it is the one produced by [lgp()].
#' @param pkl_file the path to and the name of the `.pkl` file to which
#'     the emulator object `obj` is saved.
#'
#' @md
#' @export
write <- function(obj, pkl_file) {
  pkg.env$dgpsi$write(obj, pkl_file)
}


#' @title Load the stored emulator object
#'
#' @description This function loads the `.pkl` file that stores the emulator.
#'
#' @param pkl_file the path to and the name of the `.pkl` file where the emulator is stored.
#'
#' @return an emulator object. See [write()].
#'
#' @md
#' @export
read <- function(pkl_file) {
  res <- pkg.env$dgpsi$read(pkl_file)
  return(res)
}


#' @title Summary of constructed GP, DGP, and linked (D)GP structures
#'
#' @description This function summarizes key information of GP, DGP, and linked (D)GP structures.
#'
#' @param obj can be one of the following:
#'     * the object produced by [kernel()].
#'     * the object produced by [gp()] (either before or after applying [train()]).
#'     * the object produced by [dgp()] before applying [train()].
#'     * the object produced by [emulator()].
#'     * the object produced by [lgp()].
#'
#' @return a table summarizing key information contained in `obj`.
#'
#' @md
#' @export
summary <- function(obj) {
  pkg.env$dgpsi$summary(obj)
}
