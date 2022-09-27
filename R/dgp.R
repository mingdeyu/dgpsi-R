#' @title Deep Gaussian process emulator construction
#'
#' @description This function builds and trains a DGP emulator.
#'
#' @param X a matrix where each row is an input training data point and each column is an input dimension.
#' @param Y a matrix containing observed training output data. The matrix has its rows being output data points and columns being
#'     output dimensions. When `likelihood` (see below) is not `NULL`, `Y` must be a matrix with only one column.
#' @param struc a list that specifies a user-defined DGP structure. It should contain *L* (the number of DGP layers) sub-lists,
#'     each of which represents a layer and contains a number of GP nodes (defined by [kernel()]) in the corresponding layer.
#'     The final layer of the DGP structure (i.e., the final sub-list in `struc`) can be a likelihood
#'     layer that contains a likelihood function (e.g., [Poisson()]). When `struc = NULL`,
#'     the DGP structure is automatically generated and can be checked by applying [summary()] to the output from [dgp()] with `training = FALSE`.
#'     If this argument is used (i.e., user provides a customized DGP structure), arguments `depth`, `name`, `lengthscale`, `nugget_est`, `nugget`,
#'     `connect`, `likelihood`, and `internal_input_idx` will NOT be used. Defaults to `NULL`.
#' @param depth number of layers (including the likelihood layer) for a DGP structure. `depth` must be at least `2`.
#'     Defaults to `2`. This argument is only used when `struc = NULL`.
#' @param name kernel function to be used. Either `"sexp"` for squared exponential kernel or
#'     `"matern2.5"` for Matérn-2.5 kernel. Defaults to `"sexp"`. This argument is only used when `struc = NULL`.
#' @param lengthscale initial lengthscales for GP nodes in the DGP emulator. It can be a single numeric value or a vector:
#' 1. if it is a single numeric value, the value will be applied as the initial lengthscales for all GP nodes in the DGP hierarchy.
#' 2. if it is a vector, each element of the vector specifies the initial lengthscales that will be applied to all GP nodes in the corresponding layer.
#'    The vector should have a length of `depth` if `likelihood = NULL` or a length of `depth - 1` if `likelihood` is not `NULL`.
#'
#' Defaults to a numeric value of `1.0`. This argument is only used when `struc = NULL`.
#' @param share a bool indicating if all input dimensions of a GP node share a common lengthscale. Defaults to `TRUE`. This argument is only used when `struc = NULL`.
#' @param nugget_est a bool or a bool vector that indicates if the nuggets of GP nodes (if any) in the final layer are to be estimated. If a single bool is
#'     provided, it will be applied to all GP nodes (if any) in the final layer. If a bool vector (which must have a length of `ncol(Y)`) is provided, each
#'     bool element in the vector will be applied to the corresponding GP node (if any) in the final layer. The value of a bool has following effects:
#' * `FALSE`: the nugget of the corresponding GP in the final layer is fixed to the corresponding value defined in `nugget` (see below).
#' * `TRUE`: the nugget of the corresponding GP in the final layer will be estimated with the initial value given by the correspondence in `nugget` (see below).
#'
#' Defaults to `FALSE`. This argument is only used when `struc = NULL`.
#' @param nugget the initial nugget value(s) of GP nodes (if any) in the final layer. If it is a single numeric value, it will be applied to all GP nodes (if any)
#'    in the final layer. If it is a vector (which must have a length of `ncol(Y)`), each numeric in the vector will be applied to the corresponding GP node
#'    (if any) in the final layer. Set `nugget` to a small value and the corresponding bool in `nugget_est` to `FASLE` for deterministic emulations where the emulator
#'    interpolates the training data points. Set `nugget` to a reasonable larger value and the corresponding bool in `nugget_est` to `TRUE` for stochastic emulations where
#'    the computer model outputs are assumed to follow a homogeneous Gaussian distribution. Defaults to `1e-6`. This argument is only used when `struc = NULL`.
#' @param connect a bool indicating whether to implement global input connection to the DGP structure. Defaults to `TRUE`.
#'     This argument is only used when `struc = NULL`.
#' @param likelihood the likelihood type of a DGP emulator:
#' 1. `NULL`: no likelihood layer is included in the emulator.
#' 2. `"Hetero"`: a heteroskedastic Gaussian likelihood layer is added for stochastic emulation where the computer model outputs are assumed to follow a heteroskedastic Gaussian distribution
#'    (i.e., the computer model outputs have varying noises).
#' 3. `"Poisson"`: a Poisson likelihood layer is added for stochastic emulation where the computer model outputs are assumed to a Poisson distribution.
#' 4. `"NegBin"`: a negative Binomial likelihood layer is added for stochastic emulation where the computer model outputs are assumed to follow a negative Binomial distribution.
#'
#' When `likelihood` is not `NULL`, the values of `nugget_est` and `nugget` are overridden by `FALSE` and `1e-6` respectively. Defaults to `NULL`. This argument is only used when `struc = NULL`.
#' @param training a bool indicating if the initialized DGP emulator will be trained.
#'     When set to `FALSE`, [dgp()] returns an untrained DGP emulator, to which one can apply [summary()] to inspect its specifications
#'     (especially when a customized `struc` is provided) or apply [predict()] to check its emulation performance before the training. Defaults to `TRUE`.
#' @param verb a bool indicating if the trace information on DGP emulator construction and training will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param check_rep a bool indicating whether to check the repetitions in the dataset, i.e., if one input
#'     position has multiple outputs. Defaults to `TRUE`.
#' @param rff a bool indicating whether to use random Fourier features to approximate the correlation matrices in training. Turning on this option could help accelerate
#'     the training when the training data is relatively large but may reduce the quality of the resulting emulator. Defaults to `FALSE`.
#' @param M the number of features to be used by random Fourier approximation. It is only used
#'     when `rff` is set to `TRUE`. Defaults to `NULL`. If it is `NULL`, `M` is automatically set to
#'     `max(100, ceiling(sqrt(nrow(X))*log(nrow(X))))`.
#' @param N number of iterations for the training. Defaults to `500`. This argument is only used when `training = TRUE`.
#' @param ess_burn number of burnin steps for the ESS-within-Gibbs
#'     at each I-step of the training. Defaults to `10`. This argument is only used when `training = TRUE`.
#' @param burnin the number of training iterations to be discarded for
#'     point estimates of model parameters. Must be smaller than the training iterations `N`. If this is not specified, only the last 25% of iterations
#'     are used. Defaults to `NULL`. This argument is only used when `training = TRUE`.
#' @param B the number of imputations to produce the later predictions. Increase the value to account for
#'     more imputation uncertainties. Decrease the value for lower imputation uncertainties but faster predictions.
#'     Defaults to `50`.
#' @param internal_input_idx column indices of `X` that are generated by the linked emulators in the preceding layers.
#'     Set `internal_input_idx = NULL` if the DGP emulator is in the first layer of a system or all columns in `X` are
#'     generated by the linked emulators in the preceding layers. Defaults to `NULL`. This argument is only used when `struc = NULL`.
#' @param linked_idx either a vector or a list of vectors:
#' * If `linked_idx` is a vector, it gives indices of columns in the pooled output matrix (formed by column-combined outputs of all
#'   emulators in the feeding layer) that feed into the DGP emulator. If the DGP emulator is in the first layer of a linked emulator system,
#'   the vector gives the column indices of the global input (formed by column-combining all input matrices of emulators in the first layer)
#'   that the DGP emulator will use. The length of the vector shall equal to the length of `internal_input_idx` when `internal_input_idx` is not `NULL`.
#' * When the DGP emulator is not in the first layer of a linked emulator system, `linked_idx` can be a list that gives the information on connections
#'   between the DGP emulator and emulators in all preceding layers. The length of the list should equal to the number of layers before
#'   the DGP emulator. Each element of the list is a vector that gives indices of columns in the pooled output matrix (formed by column-combined outputs
#'   of all emulators) in the corresponding layer that feed into the DGP emulator. If the DGP emulator has no connections to any emulator in a certain layer,
#'   set `NULL` in the corresponding position of the list. The order of input dimensions in `X[,internal_input_idx]` should be consistent with `linked_idx`.
#'   For example, a DGP emulator in the 4th-layer that is fed by the output dimension 2 and 4 of emulators in layer 2 and all output dimension 1 to 3 of
#'   emulators in layer 3 should have `linked_idx = list( NULL, c(2,4), c(1,2,3) )`. In addition, the first and second columns of `X[,internal_input_idx]`
#'   should correspond to the output dimensions 2 and 4 from layer 2, and the third to fifth columns of `X[,internal_input_idx]` should
#'   correspond to the output dimensions 1 to 3 from layer 3.
#'
#' Set `linked_idx = NULL` if the DGP emulator will not be used for linked emulations. However, if this is no longer the case, one can use [set_linked_idx()]
#' to add linking information to the DGP emulator. Defaults to `NULL`.
#'
#' @return An S3 class named `dgp` that contains three slots:
#' * `constructor_obj`: a 'python' object that stores the information of the constructed DGP emulator.
#' * `container_obj`: a 'python' object that stores the information for the linked emulation.
#' * `emulator_obj`: a 'python' object that stores the information for the predictions from the DGP emulator.
#'
#' The returned `dgp` object can be used by
#' * [predict()] for DGP predictions.
#' * [continue()] for additional DGP training iterations.
#' * [validate()] for LOO and OOS validations.
#' * [plot()] for validation plots.
#' * [lgp()] for linked (D)GP emulator constructions.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/> and learn how to customize a DGP structure.
#' @note Any R vector detected in `X` and `Y` will be treated as a column vector and automatically converted into a single-column
#'     R matrix.
#' @examples
#' \dontrun{
#'
#' # load the package and the Python env
#' library(dgpsi)
#' init_py()
#'
#' # construct a step function
#' f <- function(x) {
#'   if (x < 0.5) return(-1)
#'   if (x >= 0.5) return(1)
#'   }
#'
#' # generate training data
#' X <- seq(0, 1, length = 10)
#' Y <- sapply(X, f)
#'
#' # training a 3-layered DGP emulator
#' m <- dgp(X, Y, depth = 3)
#'
#' # continue for further training iterations
#' m <- continue(m)
#'
#' # summarizing
#' summary(m)
#'
#' # trace plot
#' trace_plot(m)
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
#' write(m, 'step_dgp')
#' m <- read('step_dgp')
#' }
#' @md
#' @export
dgp <- function(X, Y, struc = NULL, depth = 2, name = 'sexp', lengthscale = 1.0, share = TRUE,
                nugget_est = FALSE, nugget = 1e-6, connect = TRUE, likelihood = NULL, training =TRUE, verb = TRUE, check_rep = TRUE, rff = FALSE, M = NULL, N = 500, ess_burn = 10,
                burnin = NULL, B = 50, internal_input_idx = NULL, linked_idx = NULL) {

  if ( !is.matrix(X)&!is.vector(X) ) stop("'X' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(Y)&!is.vector(Y) ) stop("'Y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(X) ) X <- as.matrix(X)
  if ( is.vector(Y) ) Y <- as.matrix(Y)

  if ( nrow(X)!=nrow(Y) ) stop("'X' and 'Y' have different number of data points.", call. = FALSE)

  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)

  N <- as.integer(N)
  B <- as.integer(B)
  ess_burn <- as.integer(ess_burn)

  if ( !is.null(M) ) {
    M <- as.integer(M)
  }

  if ( !is.null(burnin) ) {
    burnin <- as.integer(burnin)
  }

  if( !is.null(linked_idx) ) {
    if ( !is.list(linked_idx) ) {
      linked_idx <- reticulate::np_array(as.integer(linked_idx - 1))
    } else {
      for ( i in 1:length(linked_idx) ){
        if ( !is.null(linked_idx[[i]]) ) linked_idx[[i]] <- reticulate::np_array(as.integer(linked_idx[[i]] - 1))
      }
    }
  }

  #If struc is NULL
  if ( is.null(struc) ) {
    depth <- as.integer(depth)
    if ( depth < 2 ) stop("'depth' must >= 2. Use gp() if you want a single-layered DGP.", call. = FALSE)

    if ( is.null(likelihood) ) {
      if ( length(lengthscale)==1 ) {
        lengthscale <- rep(lengthscale, depth)
      } else {
        if ( length(lengthscale)!=depth ) {
          stop(sprintf("length(lengthscale) must equal to %i.", depth), call. = FALSE)
        }
      }

      if ( length(nugget_est)==1 ) {
        nugget_est <- rep(nugget_est, n_dim_Y)
      } else {
        if ( length(nugget_est)!=n_dim_Y ) {
          stop(sprintf("length(nugget_est) should equal to %i.", n_dim_Y), call. = FALSE)
        }
      }

      if ( length(nugget)==1 ) {
        nugget <- rep(nugget, n_dim_Y)
      } else {
        if ( length(nugget)!=n_dim_Y ) {
          stop(sprintf("length(nugget) should equal to %i.", n_dim_Y), call. = FALSE)
        }
      }

      no_gp_layer = depth
    } else {
      if ( length(lengthscale)==1 ) {
        lengthscale <- rep(lengthscale, depth-1)
      } else {
        if ( length(lengthscale)!=depth-1 ) {
          stop(sprintf("length(lengthscale) must equal to %i.", depth-1), call. = FALSE)
        }
      }

      if ( likelihood == 'Hetero'|likelihood == 'NegBin' ) {
        nugget <- rep(1e-6, 2)
        nugget_est <- rep(FALSE, 2)
      } else if ( likelihood == 'Poisson' ) {
        nugget <- 1e-6
        nugget_est <- FALSE
      }

      if ( n_dim_Y != 1 ) {
        stop("'Y' must be a vector or a matrix with only one column when 'likelihood' is not NULL.", call. = FALSE)
      }

      no_gp_layer = depth - 1
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

    if ( isTRUE(verb) ) message(sprintf("Auto-generating a %i-layered DGP structure ...", depth ), appendLF = FALSE)

    struc <- list()
    for ( l in 1:no_gp_layer ) {
      layer_l <- list()

      if( l == no_gp_layer ) {
        no_kerenl <- length(nugget)
        } else {
          no_kerenl <- n_dim_X
       }

      for ( k in 1:no_kerenl ) {
        if ( l == no_gp_layer ) {
          if ( no_gp_layer == 1 ) {
            if ( isTRUE(share) ){
              length_scale <- lengthscale[l]
            } else {
              length_scale <- rep(lengthscale[l], n_dim_X)
            }
            layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), name = name, scale_est = TRUE, nugget = nugget[k], nugget_est = nugget_est[k],
                                                 input_dim = internal_input_idx, connect = external_input_idx)
          } else {
            if ( connect ) {
              if ( isTRUE(share) ){
                length_scale <- lengthscale[l]
              } else {
                length_scale <- rep(lengthscale[l], 2*n_dim_X)
              }
              layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), name = name, scale_est = TRUE, nugget = nugget[k], nugget_est = nugget_est[k],
                                                   connect = reticulate::np_array(as.integer(1:n_dim_X - 1)))
            } else {
              if ( isTRUE(share) ){
                length_scale <- lengthscale[l]
              } else {
                length_scale <- rep(lengthscale[l], n_dim_X)
              }
              layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), name = name, scale_est = TRUE, nugget = nugget[k], nugget_est = nugget_est[k])
            }
          }
        } else {
            if ( l == 1 ) {
              if ( isTRUE(share) ){
                length_scale <- lengthscale[l]
              } else {
                length_scale <- rep(lengthscale[l], n_dim_X)
              }
              layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), name = name,
                                                   input_dim = internal_input_idx, connect = external_input_idx)
            } else {
              if ( connect ) {
                if ( isTRUE(share) ){
                  length_scale <- lengthscale[l]
                } else {
                  length_scale <- rep(lengthscale[l], 2*n_dim_X)
                }
                layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), name = name,
                                                     connect = reticulate::np_array(as.integer(1:n_dim_X - 1)))
              } else {
                if ( isTRUE(share) ){
                  length_scale <- lengthscale[l]
                } else {
                  length_scale <- rep(lengthscale[l], n_dim_X)
                }
                layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), name = name)
              }
            }
          }
        }
      struc[[l]] <- layer_l
    }
    if ( !is.null(likelihood) ){
      if ( likelihood == 'Poisson' ) {
        struc[[depth]] <- c(pkg.env$dgpsi$Poisson())
      } else if ( likelihood == 'Hetero' ) {
        struc[[depth]] <- c(pkg.env$dgpsi$Hetero())
      } else if ( likelihood == 'NegBin' ) {
        struc[[depth]] <- c(pkg.env$dgpsi$NegBin())
      }
    }
    if ( isTRUE(verb) ) {
      message(" done")
      Sys.sleep(0.5)
    }
  }

  if ( isTRUE(verb) ) message("Initializing the DGP emulator ...", appendLF = FALSE)

  obj <- pkg.env$dgpsi$dgp(X, Y, struc, check_rep, rff, M)

  if ( isTRUE(verb) ) {
    message(" done")
    Sys.sleep(0.5)
  }

  if ( training ) {
    if ( isTRUE(verb) ){
      disable <- FALSE
      message("Training the DGP emulator: ")
    } else {
      disable <- TRUE
    }
    obj$train(N, ess_burn, disable)
    est_obj <- obj$estimate(burnin)
  } else {
    est_obj <- obj$estimate(NULL)
  }

  emu_obj <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B)

  res <- list()
  res[['constructor_obj']] <- obj
  res[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx)
  res[['emulator_obj']] <- emu_obj

  class(res) <- "dgp"
  return(res)
}


#' @title Continue the training of a DGP emulator
#'
#' @description This function implements additional training iterations for a DGP emulator.
#'
#' @param object an instance of the `dgp` class.
#' @param N additional number of iterations for the DGP emulator training. Defaults to `500`.
#' @param ess_burn number of burnin steps for the ESS-within-Gibbs
#'     at each I-step of the training. Defaults to `10`.
#' @param verb a bool indicating if the progress bar will be printed during the training:
#' 1. `FALSE`: the training progress bar will not be displayed.
#' 2. `TRUE`: the training progress bar will be displayed.
#'
#' Defaults to `TRUE`.
#' @param burnin the number of training iterations to be discarded for
#'     point estimates calculation. Must be smaller than the overall training iterations
#'     so-far implemented. If this is not specified, only the last 25% of iterations
#'     are used. This overrides the value of `burnin` set in [dgp()]. Defaults to `NULL`.
#' @param B the number of imputations to produce the predictions. Increase the value to account for
#'     more imputation uncertainties. This overrides the value of `B` set in [dgp()]. Defaults to `50`.
#'
#' @return An updated `object`.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note One can also use this function to fit an untrained DGP emulator constructed by [dgp()] with `training = FALSE`.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export

continue <- function(object, N = 500, ess_burn = 10, verb = TRUE, burnin = NULL, B = 50) {
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  N <- as.integer(N)
  B <- as.integer(B)
  ess_burn <- as.integer(ess_burn)

  if( !is.null(burnin) ) {
    burnin <- as.integer(burnin)
  }

  if ( verb ) {
    disable <- FALSE
    message("Continue the training:")
  } else {
    disable <- TRUE
  }

  linked_idx <- object$container_obj$local_input_idx
  object$constructor_obj$train(N, ess_burn, disable)
  est_obj <- object$constructor_obj$estimate(burnin)

  new_object <- list()
  new_object[['constructor_obj']] <- object$constructor_obj
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx)
  class(new_object) <- "dgp"
  return(new_object)
}
