#' @title Update a GP or DGP emulator
#'
#' @description This function updates the training input and output of a GP or DGP emulator with an option to refit the emulator.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' @param X the new input data which is a matrix where each row is an input training data point and each column is an input dimension.
#' @param Y the new output data:
#' * If `object` is an instance of the `gp` class, `Y` is a matrix with only one column and each row being an output data point.
#' * If `object` is an instance of the `dgp` class, `Y` is a matrix with its rows being output data points and columns being
#'     output dimensions. When `likelihood` (see below) is not `NULL`, `Y` must be a matrix with only one column.
#' @param refit a bool indicating whether to re-fit the emulator `object` after the training input and output are updated. Defaults to `FALSE`.
#' @param verb a bool indicating if the trace information will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param N number of training iterations used to re-fit the emulator `object` if it is an instance of the `dgp` class. Defaults to `100`.
#' @param ess_burn number of burnin steps for the ESS-within-Gibbs at each I-step in training the emulator `object` if it is an
#'     instance of the `dgp` class. Defaults to `10`.
#' @param B the number of imputations for predictions from the updated emulator `object` if it is an instance of the `dgp` class.
#'     This overrides the number of imputations set in `object`. Set to `NULL` to use the same number of imputations set
#'     in `object`. Defaults to `NULL`.
#' @param ... N/A.
#'
#' @return An updated `object`.
#'
#' @note
#' * The following slots:
#'   - `loo` and `oos` created by [validate()];
#'   - `results` created by [predict()]; and
#'   - `design` created by [design()]
#'   in `object` will be removed and not contained in the returned object.
#' * Any R vector detected in `X` and `Y` will be treated as a column vector and automatically converted into a single-column
#'   R matrix. Thus, if `X` is a single data point with multiple dimensions, it must be given as a matrix.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See locate() for an example.
#' }
#' @md
#' @name update
#' @export
update <- function(object, X, Y, refit, verb, ...){
  UseMethod("update")
}

#' @rdname update
#' @method update dgp
#' @export
update.dgp <- function(object, X, Y, refit = FALSE, verb = TRUE, N = 100, ess_burn = 10, B = NULL, ...) {
  #check class
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  N <- as.integer(N)
  ess_burn <- as.integer(ess_burn)
  if ( is.null(B) ){
    B <- as.integer(length(object$emulator_obj$all_layer_set))
  } else {
    B <- as.integer(B)
  }

  if ( !is.matrix(X)&!is.vector(X) ) stop("'X' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(Y)&!is.vector(Y) ) stop("'Y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(X) ) X <- as.matrix(X)
  if ( is.vector(Y) ) Y <- as.matrix(Y)

  if ( nrow(X)!=nrow(Y) ) stop("'X' and 'Y' have different number of data points.", call. = FALSE)

  linked_idx <- object$container_obj$local_input_idx

  if ( verb ) message("Updating ...", appendLF = FALSE)
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  constructor_obj_cp$update_xy(X,Y)
  if ( verb ) {
    Sys.sleep(0.2)
    message(" done")
  }

  if ( refit ){
    if ( verb ) {
      disable <- FALSE
      message("Re-fitting:")
    } else {
      disable <- TRUE
    }
    N0 <- constructor_obj_cp$N
    constructor_obj_cp$train(N, ess_burn, disable)
    burnin <- N0 + as.integer(0.75*N)
  } else {
    burnin <- constructor_obj_cp$burnin
  }

  if ( isTRUE(verb) ) message("Imputing ...", appendLF = FALSE)
  est_obj <- constructor_obj_cp$estimate(burnin)

  new_object <- list()
  new_object[['data']][['X']] <- unname(X)
  new_object[['data']][['Y']] <- unname(Y)
  new_object[['constructor_obj']] <- constructor_obj_cp
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx)
  class(new_object) <- "dgp"
  if ( isTRUE(verb) ) message(" done")

  return(new_object)
}

#' @rdname update
#' @method update gp
#' @export
update.gp <- function(object, X, Y, refit = FALSE, verb = TRUE, ...) {
  #check class
  if ( !inherits(object,"gp") ){
    stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  }
  if ( !is.matrix(X)&!is.vector(X) ) stop("'X' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(Y)&!is.vector(Y) ) stop("'Y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(X) ) X <- as.matrix(X)
  if ( is.vector(Y) ) Y <- as.matrix(Y)

  if ( nrow(X)!=nrow(Y) ) stop("'X' and 'Y' have different number of data points.", call. = FALSE)

  if ( ncol(Y) != 1 ) {
    stop("'Y' must be a vector or a matrix with only one column for a GP emulator.", call. = FALSE)
  }

  linked_idx <- object$container_obj$local_input_idx

  if ( verb ) message("Updating ...", appendLF = FALSE)
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  constructor_obj_cp$update_xy(X,Y)
  if ( verb ) {
    Sys.sleep(0.5)
    message(" done")
  }

  if ( refit ){
    if ( verb ) message("Re-fitting ...", appendLF = FALSE)
    constructor_obj_cp$train()
    if ( verb ) message(" done")
  }

  new_object <- list()
  new_object[['data']][['X']] <- unname(X)
  new_object[['data']][['Y']] <- unname(Y)
  new_object[['constructor_obj']] <- constructor_obj_cp
  new_object[['container_obj']] <- pkg.env$dgpsi$container(constructor_obj_cp$export(), linked_idx)
  new_object[['emulator_obj']] <- constructor_obj_cp
  class(new_object) <- "gp"

  return(new_object)
}



