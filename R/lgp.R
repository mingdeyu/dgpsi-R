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
#' @return An S3 class named `lgp` that contains a slot called `emulator_obj`, which is a 'python' object that
#'     stores the information for predictions from the linked emulator. The returned `lgp` object can be used by
#' * [predict()] for linked (D)GP predictions.
#' * [validate()] for the OOS validation.
#' * [plot()] for the validation plots.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # load the package and the Python env
#' library(dgpsi)
#' init_py()
#'
#' # model 1
#' f1 <- function(x) {
#'   (sin(7.5*x)+1)/2
#' }
#' # model 2
#' f2 <- function(x) {
#'   2/3*sin(2*(2*x - 1))+4/3*exp(-30*(2*(2*x-1))^2)-1/3
#' }
#' # linked model
#' f12 <- function(x) {
#'   f2(f1(x))
#' }
#'
#' # training data for Model 1
#' X1 <- seq(0, 1, length = 9)
#' Y1 <- sapply(X1, f1)
#' # training data for Model 2
#' X2 <- seq(0, 1, length = 13)
#' Y2 <- sapply(X2, f2)
#'
#' # emulation of model 1
#' m1 <- gp(X1, Y1, name = "matern2.5", linked_idx = c(1))
#' # emulation of model 2
#' m2 <- dgp(X2, Y2, depth = 2, name = "matern2.5")
#' # assign linking information after the emulation construction
#' m2 <- set_linked_idx(m2, c(1))
#'
#' # emulation of the linked model
#' struc <- combine(list(m1), list(m2))
#' m_link <- lgp(struc)
#'
#' # summarizing
#' summary(m_link)
#'
#' # prediction
#' test_x <- seq(0, 1, length = 300)
#' m_link <- predict(m_link, x = test_x)
#'
#' # OOS validation
#' validate_x <- sample(test_x, 20)
#' validate_y <- sapply(validate_x, f12)
#' plot(m_link, validate_x, validate_y, style = 2)
#'
#' # write and read the constructed linked emulator
#' write(m_link, 'linked_emulator')
#' m_link <- read('linked_emulator')
#' }
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
