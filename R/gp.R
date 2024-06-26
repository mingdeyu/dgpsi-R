#' @title Gaussian process emulator construction
#'
#' @description This function builds and trains a GP emulator.
#'
#' @param X a matrix where each row is an input data point and each column is an input dimension.
#' @param Y a matrix with only one column and each row being an output data point.
#' @param struc an object produced by [kernel()] that gives a user-defined GP specifications. When `struc = NULL`,
#'     the GP specifications are automatically generated using information provided in `name`, `lengthscale`,
#'     `nugget_est`, `nugget`, `scale_est`, `scale`,and `internal_input_idx`. Defaults to `NULL`.
#' @param name kernel function to be used. Either `"sexp"` for squared exponential kernel or
#'     `"matern2.5"` for Matérn-2.5 kernel. Defaults to `"sexp"`. This argument is only used when `struc = NULL`.
#' @param lengthscale initial values of lengthscales in the kernel function. It can be a single numeric value or a vector:
#' * if it is a single numeric value, it is assumed that kernel functions across input dimensions share the same lengthscale;
#' * if it is a vector (which must have a length of `ncol(X)`), it is assumed that kernel functions across input dimensions have different lengthscales.
#'
#' Defaults to a vector of `0.1`. This argument is only used when `struc = NULL`.
#' @param bounds the lower and upper bounds of lengthscales in the kernel function. It is a vector of length two where the first element is
#'    the lower bound and the second element is the upper bound. The bounds will be applied to all lengthscales in the kernel function. Defaults
#'    to `NULL` where no bounds are specified for the lengthscales. This argument is only used when `struc = NULL`.
#' @param prior prior to be used for Maximum a Posterior for lengthscales and nugget of the GP: gamma prior (`"ga"`), inverse gamma prior (`"inv_ga"`),
#'     or jointly robust prior (`"ref"`). Defaults to `"ref"`. This argument is only used when `struc = NULL`. See the reference below for the jointly
#'     robust prior.
#' @param nugget_est a bool indicating if the nugget term is to be estimated:
#' 1. `FALSE`: the nugget term is fixed to `nugget`.
#' 2. `TRUE`: the nugget term will be estimated.
#'
#' Defaults to `FALSE`. This argument is only used when `struc = NULL`.
#' @param nugget the initial nugget value. If `nugget_est = FALSE`, the assigned value is fixed during the training.
#'     Set `nugget` to a small value (e.g., `1e-8`) and the corresponding bool in `nugget_est` to `FASLE` for deterministic emulations where the emulator
#'     interpolates the training data points. Set `nugget` to a reasonable larger value and the corresponding bool in `nugget_est` to `TRUE` for stochastic
#'     emulations where the computer model outputs are assumed to follow a homogeneous Gaussian distribution. Defaults to `1e-8` if `nugget_est = FALSE` and
#'     `0.01` if `nugget_est = TRUE`. This argument is only used when `struc = NULL`.
#' @param scale_est a bool indicating if the variance is to be estimated:
#' 1. `FALSE`: the variance is fixed to `scale`.
#' 2. `TRUE`: the variance term will be estimated.
#'
#' Defaults to `TRUE`. This argument is only used when `struc = NULL`.
#' @param scale the initial variance value. If `scale_est = FALSE`, the assigned value is fixed during the training.
#'     Defaults to `1`. This argument is only used when `struc = NULL`.
#' @param training a bool indicating if the initialized GP emulator will be trained.
#'     When set to `FALSE`, [gp()] returns an untrained GP emulator, to which one can apply [summary()] to inspect its specifications
#'     (especially when a customized `struc` is provided) or apply [predict()] to check its emulation performance before the training. Defaults to `TRUE`.
#' @param verb a bool indicating if the trace information on GP emulator construction and training will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param vecchia a bool indicating whether to use Vecchia approximation for large-scale GP emulator construction and prediction. Defaults to `FALSE`.
#'     The Vecchia approximation implemented for the GP emulation largely follows Katzfuss et al. (2022). See reference below.
#' @param M the size of the conditioning set for the Vecchia approximation in the GP emulator training. Defaults to `25`.
#' @param ord an R function that returns the ordering of the input to the GP emulator for the Vecchia approximation. The function must satisfy the following basic rules:
#' * the first argument represents the input scaled by the lengthscales.
#' * the output of the function is a vector of indices that gives the ordering of the input to the GP emulator.
#'
#' If `ord = NULL`, the default random ordering is used. Defaults to `NULL`.
#' @param internal_input_idx the column indices of `X` that are generated by the linked emulators in the preceding layers.
#'     Set `internal_input_idx = NULL` if the GP emulator is in the first layer of a system or all columns in `X` are
#'     generated by the linked emulators in the preceding layers. Defaults to `NULL`. This argument is only used when `struc = NULL`.
#' @param linked_idx either a vector or a list of vectors:
#' * If `linked_idx` is a vector, it gives indices of columns in the pooled output matrix (formed by column-combined outputs of all
#'   emulators in the feeding layer) that feed into the GP emulator. The length of the vector shall equal to the length of `internal_input_idx` when `internal_input_idx`
#'   is not `NULL`. If the GP emulator is in the first layer of a linked emulator system, the vector gives the column indices of the global
#'   input (formed by column-combining all input matrices of emulators in the first layer) that the GP emulator will use. If the GP emulator is to be used in both the first
#'   and subsequent layers, one should initially set `linked_idx` to the appropriate values for the situation where the emulator is not in the first layer. Then, use the
#'   function [set_linked_idx()] to reset the linking information when the emulator is in the first layer.
#' * When the GP emulator is not in the first layer of a linked emulator system, `linked_idx` can be a list that gives the information on connections
#'   between the GP emulator and emulators in all preceding layers. The length of the list should equal to the number of layers before
#'   the GP emulator. Each element of the list is a vector that gives indices of columns in the pooled output matrix (formed by column-combined outputs
#'   of all emulators) in the corresponding layer that feed into the GP emulator. If the GP emulator has no connections to any emulator in a certain layer,
#'   set `NULL` in the corresponding position of the list. The order of input dimensions in `X[,internal_input_idx]` should be consistent with `linked_idx`.
#'   For example, a GP emulator in the second layer that is fed by the output dimension 1 and 3 of emulators in layer 1 should have `linked_idx = list( c(1,3) )`.
#'   In addition, the first and second columns of `X[,internal_input_idx]` should correspond to the output dimensions 1 and 3 from layer 1.
#'
#' Set `linked_idx = NULL` if the GP emulator will not be used for linked emulations. However, if this is no longer the case, one can use [set_linked_idx()]
#' to add linking information to the GP emulator. Defaults to `NULL`.
#' @param id an ID to be assigned to the GP emulator. If an ID is not provided (i.e., `id = NULL`), a UUID (Universally Unique Identifier) will be automatically generated
#'    and assigned to the emulator. Default to `NULL`.
#'
#' @return An S3 class named `gp` that contains five slots:
#' * `id`: A number or character string assigned through the `id` argument.
#' * `data`: a list that contains two elements: `X` and `Y` which are the training input and output data respectively.
#' * `specs`: a list that contains seven elements:
#'    1. `kernel`: the type of the kernel function used. Either `"sexp"` for squared exponential kernel or `"matern2.5"` for Matérn-2.5 kernel.
#'    2. `lengthscales`: a vector of lengthscales in the kernel function.
#'    3. `scale`: the variance value in the kernel function.
#'    4. `nugget`: the nugget value in the kernel function.
#'    5. `internal_dims`: the column indices of `X` that correspond to the linked emulators in the preceding layers of a linked system.
#'    6. `external_dims`: the column indices of `X` that correspond to global inputs to the linked system of emulators. It is shown as `FALSE` if `internal_input_idx = NULL`.
#'    7. `linked_idx`: the value passed to argument `linked_idx`. It is shown as `FALSE` if the argument `linked_idx` is `NULL`.
#'    8. `vecchia`: whether the Vecchia approximation is used for the GP emulator training.
#'    9. `M`: the size of the conditioning set for the Vecchia approximation in the GP emulator training.
#'
#'    `internal_dims` and `external_dims` are generated only when `struc = NULL`.
#' * `constructor_obj`: a 'python' object that stores the information of the constructed GP emulator.
#' * `container_obj`: a 'python' object that stores the information for the linked emulation.
#' * `emulator_obj`: a 'python' object that stores the information for the predictions from the GP emulator.
#'
#' The returned `gp` object can be used by
#' * [predict()] for GP predictions.
#' * [validate()] for LOO and OOS validations.
#' * [plot()] for validation plots.
#' * [lgp()] for linked (D)GP emulator constructions.
#' * [summary()] to summarize the trained GP emulator.
#' * [write()] to save the GP emulator to a `.pkl` file.
#' * [set_linked_idx()] to add the linking information to the GP emulator for linked emulations.
#' * [design()] for sequential designs.
#' * [update()] to update the GP emulator with new inputs and outputs.
#' * [alm()], [mice()], [pei()], and [vigf()] to locate next design points.
#'
#' @references
#' - Gu, M. (2019). Jointly robust prior for Gaussian stochastic process in emulation, calibration and variable selection. *Bayesian Analysis*, **14(3)**, 857-885.
#' - Katzfuss, M., Guinness, J., & Lawrence, E. (2022). Scaled Vecchia approximation for fast computer-model emulation. *SIAM/ASA Journal on Uncertainty Quantification*, **10(2)**, 537-554.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note Any R vector detected in `X` and `Y` will be treated as a column vector and automatically converted into a single-column
#'     R matrix. Thus, if `X` is a single data point with multiple dimensions, it must be given as a matrix.
#' @examples
#' \dontrun{
#' # load the package and the Python env
#' library(dgpsi)
#'
#' # construct a step function
#' f <- function(x) {
#'    if (x < 0.5) return(-1)
#'    if (x >= 0.5) return(1)
#'   }
#'
#' # generate training data
#' X <- seq(0, 1, length = 10)
#' Y <- sapply(X, f)
#'
#' # training
#' m <- gp(X, Y)
#'
#' # summarizing
#' summary(m)
#'
#' # LOO cross validation
#' m <- validate(m)
#' plot(m)
#'
#' # prediction
#' test_x <- seq(0, 1, length = 200)
#' m <- predict(m, x = test_x)
#'
#' # OOS validation
#' validate_x <- sample(test_x, 10)
#' validate_y <- sapply(validate_x, f)
#' plot(m, validate_x, validate_y)
#'
#' # write and read the constructed emulator
#' write(m, 'step_gp')
#' m <- read('step_gp')
#' }
#'
#' @md
#' @export
gp <- function(X, Y, struc = NULL, name = 'sexp', lengthscale = rep(0.1, ncol(X)), bounds = NULL, prior = 'ref', nugget_est = FALSE, nugget = ifelse(nugget_est, 0.01, 1e-8), scale_est = TRUE, scale = 1., training = TRUE, verb = TRUE, vecchia = FALSE, M = 25, ord = NULL, internal_input_idx = NULL, linked_idx = NULL, id = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !is.matrix(X)&!is.vector(X) ) stop("'X' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(Y)&!is.vector(Y) ) stop("'Y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(X) ) X <- as.matrix(X)
  if ( is.vector(Y) ) Y <- as.matrix(Y)

  if ( nrow(X)!=nrow(Y) ) stop("'X' and 'Y' have different number of data points.", call. = FALSE)
  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)
  if ( n_dim_Y != 1 ) {
    stop("'Y' must be a vector or a matrix with only one column for a GP emulator.", call. = FALSE)
  }

  M <- as.integer(M)

  if ( !is.null(ord) ) {
    ord_wrapper <- function(x) {
      return( as.integer(ord(x) - 1) )
    }
    ord_wrapper <- reticulate::py_func(ord_wrapper)
  } else {
    ord_wrapper <- NULL
  }

  if ( is.null(struc) ) {
    is.null.struc <- TRUE
  } else {
    is.null.struc <- FALSE
  }

  if ( name!='sexp' & name!='matern2.5' ) stop("'name' can only be either 'sexp' or 'matern2.5'.", call. = FALSE)
  if ( prior!='ga' & prior!='inv_ga' & prior!='ref') stop("'prior' can only be 'ga', 'inv_ga', or 'ref'.", call. = FALSE)

  linked_idx_py <- linked_idx_r_to_py(linked_idx)

  if ( is.null.struc ) {
    if ( verb ) message("Auto-generating a GP structure ...", appendLF = FALSE)

    if ( length(lengthscale) != 1 & length(lengthscale) != n_dim_X) {
      stop("length(lengthscale) must be 1 or ncol(X).", call. = FALSE)
    }

    if ( !is.null(bounds) ){
      if ( !is.vector(bounds) ) {
        bounds <- as.vector(bounds)
      }
      if ( length(bounds)!=2 ) {
        stop(sprintf("length(bounds) must equal to %i.", 2), call. = FALSE)
      }
      if ( bounds[1]>bounds[2] ) stop("The second element of 'bounds' must be greater than the first one.", call. = FALSE)
      bounds <- reticulate::np_array(bounds)
    }

    if( !is.null(internal_input_idx) ) {
      external_input_idx <- setdiff(1:n_dim_X, internal_input_idx)
      if ( length(external_input_idx) == 0) {
        internal_input_idx = NULL
        external_input_idx = NULL
      } else {
        internal_input_idx <- reticulate::np_array(as.integer(internal_input_idx - 1))
        external_input_idx <- reticulate::np_array(as.integer(external_input_idx - 1))
      }
    } else {
      external_input_idx = NULL
    }

    struc <- pkg.env$dgpsi$kernel(length = reticulate::np_array(lengthscale), name = name, prior_name = prior, bds = bounds, scale = scale, scale_est = scale_est, nugget = nugget, nugget_est = nugget_est,
                                  input_dim = internal_input_idx, connect = external_input_idx)

    if ( verb ) {
      message(" done")
      Sys.sleep(0.5)
    }
  }

  if ( verb ) message("Initializing the GP emulator ...", appendLF = FALSE)

  obj <- pkg.env$dgpsi$gp(X, Y, struc, vecchia, M, ord_wrapper)

  if ( verb ) {
    Sys.sleep(0.5)
    message(" done")
  }

  if ( training ) {
    if ( verb ){
      Sys.sleep(0.5)
      message("Training the GP emulator ...", appendLF = FALSE)
    }
    obj$train()
    if ( verb ) message(" done")
  }

  res <- list()
  res[['id']] <- if (is.null(id)) uuid::UUIDgenerate() else id
  res[['data']][['X']] <- unname(X)
  res[['data']][['Y']] <- unname(Y)
  res[['specs']] <- extract_specs(obj, "gp")
  if ( is.null.struc ) {
    res[['specs']][['internal_dims']] <- if( is.null(internal_input_idx) ) 1:n_dim_X else as.integer(reticulate::py_to_r(internal_input_idx)+1)
    res[['specs']][['external_dims']] <- if( is.null(internal_input_idx) ) FALSE else as.integer(reticulate::py_to_r(external_input_idx)+1)
  }
  res[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx
  res[['specs']][['vecchia']] <- vecchia
  res[['specs']][['M']] <- M
  res[['constructor_obj']] <- obj
  res[['container_obj']] <- pkg.env$dgpsi$container(obj$export(), linked_idx_py)
  res[['emulator_obj']] <- obj

  class(res) <- "gp"

  return(res)
}
