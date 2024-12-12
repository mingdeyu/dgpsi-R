#' @title Deep Gaussian process emulator construction
#'
#' @description This function builds and trains a DGP emulator.
#'
#' @param X a matrix where each row is an input training data point and each column represents an input dimension.
#' @param Y a matrix containing observed training output data. The matrix has its rows being output data points and columns representing
#'     output dimensions. When `likelihood` (see below) is not `NULL`, `Y` must be a matrix with a single column.
#' @param depth number of layers (including the likelihood layer) for a DGP structure. `depth` must be at least `2`.
#'     Defaults to `2`.
#' @param node number of GP nodes in each layer (except for the final layer or the layer feeding the likelihood node) of the DGP. Defaults to
#'    `ncol(X)`.
#' @param name a character or a vector of characters that indicates the kernel functions (either `"sexp"` for squared exponential kernel or
#'     `"matern2.5"` for Matérn-2.5 kernel) used in the DGP emulator:
#' 1. if a single character is supplied, the corresponding kernel function will be used for all GP nodes in the DGP hierarchy.
#' 2. if a vector of characters is supplied, each character of the vector specifies the kernel function that will be applied to all GP nodes in the corresponding layer.
#'
#' Defaults to `"sexp"`.
#' @param lengthscale initial lengthscales for GP nodes in the DGP emulator. It can be a single numeric value or a vector:
#' 1. if it is a single numeric value, the value will be applied as the initial lengthscales for all GP nodes in the DGP hierarchy.
#' 2. if it is a vector, each element of the vector specifies the initial lengthscales that will be applied to all GP nodes in the corresponding layer.
#'    The vector should have a length of `depth` if `likelihood = NULL` or a length of `depth - 1` if `likelihood` is not `NULL`.
#'
#' Defaults to a numeric value of `1.0`.
#' @param bounds the lower and upper bounds of lengthscales in GP nodes. It can be a vector or a matrix:
#' 1. if it is a vector, the lower bound (the first element of the vector) and upper bound (the second element of the vector) will be applied to
#'    lengthscales for all GP nodes in the DGP hierarchy.
#' 2. if it is a matrix, each row of the matrix specifies the lower and upper bounds of lengthscales for all GP nodes in the corresponding layer.
#'    The matrix should have its row number equal to `depth` if `likelihood = NULL` or to `depth - 1` if `likelihood` is not `NULL`.
#'
#' Defaults to `NULL` where no bounds are specified for the lengthscales.
#' @param prior prior to be used for MAP estimation of lengthscales and nuggets of all GP nodes in the DGP hierarchy:
#' * gamma prior (`"ga"`),
#' * inverse gamma prior (`"inv_ga"`), or
#' * jointly robust prior (`"ref"`).
#'
#' Defaults to `"ga"`.
#' @param share a bool indicating if all input dimensions of a GP node share a common lengthscale. Defaults to `TRUE`.
#' @param nugget_est a bool or a bool vector that indicates if the nuggets of GP nodes (if any) in the final layer are to be estimated. If a single bool is
#'     provided, it will be applied to all GP nodes (if any) in the final layer. If a bool vector (which must have a length of `ncol(Y)`) is provided, each
#'     bool element in the vector will be applied to the corresponding GP node (if any) in the final layer. The value of a bool has following effects:
#' * `FALSE`: the nugget of the corresponding GP in the final layer is fixed to the corresponding value defined in `nugget` (see below).
#' * `TRUE`: the nugget of the corresponding GP in the final layer will be estimated with the initial value given by the correspondence in `nugget` (see below).
#'
#' Defaults to `FALSE`.
#' @param nugget the initial nugget value(s) of GP nodes (if any) in each layer:
#' 1. if it is a single numeric value, the value will be applied as the initial nugget for all GP nodes in the DGP hierarchy.
#' 2. if it is a vector, each element of the vector specifies the initial nugget that will be applied to all GP nodes in the corresponding layer.
#'    The vector should have a length of `depth` if `likelihood = NULL` or a length of `depth - 1` if `likelihood` is not `NULL`.
#'
#' Set `nugget` to a small value and the bools in `nugget_est` to `FALSE` for deterministic emulation, where the emulator
#'    interpolates the training data points. Set `nugget` to a larger value and the bools in `nugget_est` to `TRUE` for stochastic emulation where
#'    the computer model outputs are assumed to follow a homogeneous Gaussian distribution. Defaults to `1e-6` if `nugget_est = FALSE` and
#'    `0.01` if `nugget_est = TRUE`. If `likelihood` is not `NULL` and `nugget_est = FALSE`, the nuggets of GPs that feed into the likelihood layer default to
#'    `1e-4`.
#' @param scale_est a bool or a bool vector that indicates if the variance of GP nodes (if any) in the final layer are to be estimated. If a single bool is
#'     provided, it will be applied to all GP nodes (if any) in the final layer. If a bool vector (which must have a length of `ncol(Y)`) is provided, each
#'     bool element in the vector will be applied to the corresponding GP node (if any) in the final layer. The value of a bool has following effects:
#' * `FALSE`: the variance of the corresponding GP in the final layer is fixed to the corresponding value defined in `scale` (see below).
#' * `TRUE`: the variance of the corresponding GP in the final layer will be estimated with the initial value given by the correspondence in `scale` (see below).
#'
#' Defaults to `TRUE`.
#' @param scale the initial variance value(s) of GP nodes (if any) in the final layer. If it is a single numeric value, it will be applied to all GP nodes (if any)
#'    in the final layer. If it is a vector (which must have a length of `ncol(Y)`), each numeric in the vector will be applied to the corresponding GP node
#'    (if any) in the final layer. Defaults to `1`.
#' @param connect a bool indicating whether to implement global input connection to the DGP structure. Setting it to `FALSE` may produce a better emulator in some cases at
#'    the cost of slower training. Defaults to `TRUE`.
#' @param likelihood the likelihood type of a DGP emulator:
#' 1. `NULL`: no likelihood layer is included in the emulator.
#' 2. `"Hetero"`: a heteroskedastic Gaussian likelihood layer is added for stochastic emulation where the computer model outputs are assumed to follow a heteroskedastic Gaussian distribution
#'    (i.e., the computer model outputs have input-dependent noise).
#' 3. `"Poisson"`: a Poisson likelihood layer is added for emulation where the computer model outputs are counts and a Poisson distribution is used to model them.
#' 4. `"NegBin"`: a negative Binomial likelihood layer is added for emulation where the computer model outputs are counts and a negative Binomial distribution is used to capture dispersion variability in input space.
#' 5. `r new_badge("new")` `"Categorical"`: a categorical likelihood layer is added for emulation (classification), where the computer model output is categorical.
#'
#' When `likelihood` is not `NULL`, the value of `nugget_est` is overridden by `FALSE`. Defaults to `NULL`.
#' @param training a bool indicating if the initialized DGP emulator will be trained.
#'     When set to `FALSE`, [dgp()] returns an untrained DGP emulator, to which one can apply [summary()] to inspect its specifications
#'     or apply [predict()] to check its emulation performance before training. Defaults to `TRUE`.
#' @param verb a bool indicating if the trace information on DGP emulator construction and training will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param check_rep a bool indicating whether to check for repetitions in the dataset, i.e., if one input
#'     position has multiple outputs. Defaults to `TRUE`.
#' @param vecchia `r new_badge("new")` a bool indicating whether to use Vecchia approximation for large-scale DGP emulator construction and prediction. Defaults to `FALSE`.
#' @param M `r new_badge("new")` the size of the conditioning set for the Vecchia approximation in the DGP emulator training. Defaults to `25`.
#' @param ord `r new_badge("new")` an R function that returns the ordering of the input to each GP node contained in the DGP emulator for the Vecchia approximation. The
#'    function must satisfy the following basic rules:
#' * the first argument represents the input to a GP node scaled by its lengthscales.
#' * the output of the function is a vector of indices that gives the ordering of the input to the GP node.
#'
#' If `ord = NULL`, the default random ordering is used. Defaults to `NULL`.
#' @param N number of iterations for the training. Defaults to `500` if `vecchia = FALSE` and `200` if `vecchia = TRUE`. This argument is only used when `training = TRUE`.
#' @param cores the number of processes to be used to optimize GP components (in the same layer) at each M-step of the training. If set to `NULL`,
#'     the number of processes is set to `(max physical cores available - 1)` if `vecchia = FALSE` and `max physical cores available %/% 2` if `vecchia = TRUE`.
#'     Only use multiple processes when there is a large number of GP components in different layers and optimization of GP components is computationally expensive. Defaults to `1`.
#' @param blocked_gibbs a bool indicating if the latent variables are imputed layer-wise using ESS-within-Blocked-Gibbs. ESS-within-Blocked-Gibbs would be faster and
#'     more efficient than ESS-within-Gibbs that imputes latent variables node-wise because it reduces the number of components to be sampled during Gibbs steps,
#'     especially when there is a large number of GP nodes in layers due to higher input dimensions. Default to `TRUE`.
#' @param ess_burn number of burnin steps for the ESS-within-Gibbs
#'     at each I-step of the training. Defaults to `10`. This argument is only used when `training = TRUE`.
#' @param burnin the number of training iterations to be discarded for
#'     point estimates of model parameters. Must be smaller than the training iterations `N`. If this is not specified, only the last 25% of iterations
#'     are used. Defaults to `NULL`. This argument is only used when `training = TRUE`.
#' @param B the number of imputations used to produce predictions. Increase the value to refine the representation of imputation uncertainty.
#'     Defaults to `10`.
#' @param internal_input_idx `r lifecycle::badge("deprecated")` The argument will be removed in the next release. To set up connections of emulators for linked emulations,
#'     please use the updated [lgp()] function instead.
#'
#' Column indices of `X` that are generated by the linked emulators in the preceding layers.
#'     Set `internal_input_idx = NULL` if the DGP emulator is in the first layer of a system or all columns in `X` are
#'     generated by the linked emulators in the preceding layers. Defaults to `NULL`.
#' @param linked_idx `r lifecycle::badge("deprecated")` The argument will be removed in the next release. To set up connections of emulators for linked emulation,
#'     please use the updated [lgp()] function instead.
#'
#' Either a vector or a list of vectors:
#' * If `linked_idx` is a vector, it gives indices of columns in the pooled output matrix (formed by column-combined outputs of all
#'   emulators in the feeding layer) that feed into the DGP emulator. The length of the vector shall equal to the length of `internal_input_idx`
#'   when `internal_input_idx` is not `NULL`. If the DGP emulator is in the first layer of a linked emulator system, the vector gives the column indices of the global
#'   input (formed by column-combining all input matrices of emulators in the first layer) that the DGP emulator will use. If the DGP emulator is to be used in both the first
#'   and subsequent layers, one should initially set `linked_idx` to the appropriate values for the situation where the emulator is not in the first layer. Then, use the
#'   function [set_linked_idx()] to reset the linking information when the emulator is in the first layer.
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
#' @param id an ID to be assigned to the DGP emulator. If an ID is not provided (i.e., `id = NULL`), a UUID (Universally Unique Identifier) will be automatically generated
#'    and assigned to the emulator. Default to `NULL`.
#'
#' @return An S3 class named `dgp` that contains five slots:
#' * `id`: A number or character string assigned through the `id` argument.
#' * `data`: a list that contains two elements: `X` and `Y` which are the training input and output data respectively.
#' * `specs`: a list that contains
#'   1. *L* (i.e., the number of layers in the DGP hierarchy) sub-lists named `layer1, layer2,..., layerL`. Each sub-list contains *D*
#'      (i.e., the number of GP/likelihood nodes in the corresponding layer) sub-lists named `node1, node2,..., nodeD`. If a sub-list
#'      corresponds to a likelihood node, it contains one element called `type` that gives the name (`Hetero`, `Poisson`, `NegBin`, or `Categorical`) of the likelihood node.
#'      If a sub-list corresponds to a GP node, it contains four elements:
#'      - `kernel`: the type of the kernel function used for the GP node.
#'      - `lengthscales`: a vector of lengthscales in the kernel function.
#'      - `scale`: the variance value in the kernel function.
#'      - `nugget`: the nugget value in the kernel function.
#'   2. `r lifecycle::badge("deprecated")` `internal_dims`: the column indices of `X` that correspond to the linked emulators in the preceding layers of a linked system.
#'      **The slot will be removed in the next release**.
#'   3. `r lifecycle::badge("deprecated")` `external_dims`: the column indices of `X` that correspond to global inputs to the linked system of emulators. It is shown
#'      as `FALSE` if `internal_input_idx = NULL`. **The slot will be removed in the next release**.
#'   4. `r lifecycle::badge("deprecated")` `linked_idx`: the value passed to argument `linked_idx`. It is shown as `FALSE` if the argument `linked_idx` is `NULL`.
#'      **The slot will be removed in the next release**.
#'   5. `seed`: the random seed generated to produce imputations. This information is stored for reproducibility when the DGP emulator (that was saved by [write()]
#'      with the light option `light = TRUE`) is loaded back to R by [read()].
#'   6. `B`: the number of imputations used to generate the emulator.
#'   7. `r new_badge("new")` `vecchia`: whether the Vecchia approximation is used for the GP emulator training.
#'   8. `r new_badge("new")` `M`: the size of the conditioning set for the Vecchia approximation in the DGP emulator training. `M` is generated only when `vecchia = TRUE`.
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
#' * [window()] for model parameter trimming.
#' * [summary()] to summarize the trained DGP emulator.
#' * [write()] to save the DGP emulator to a `.pkl` file.
#' * [set_imp()] to change the number of imputations.
#' * [design()] for sequential design.
#' * [update()] to update the DGP emulator with new inputs and outputs.
#' * [alm()], [mice()], and [vigf()] to locate next design points.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @note Any R vector detected in `X` and `Y` will be treated as a column vector and automatically converted into a single-column
#'     R matrix. Thus, if `X` is a single data point with multiple dimensions, it must be given as a matrix.
#' @examples
#' \dontrun{
#'
#' # load the package and the Python env
#' library(dgpsi)
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
#' # set a random seed
#' set_seed(999)
#'
#' # training a DGP emulator
#' m <- dgp(X, Y)
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
#' # trim the traces of model parameters
#' m <- window(m, 800)
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
dgp <- function(X, Y, depth = 2, node = ncol(X), name = 'sexp', lengthscale = 1.0, bounds = NULL, prior = 'ga', share = TRUE,
                nugget_est = FALSE, nugget = NULL, scale_est = TRUE, scale = 1., connect = TRUE,
                likelihood = NULL, training =TRUE, verb = TRUE, check_rep = TRUE, vecchia = FALSE, M = 25, ord = NULL, N = ifelse(vecchia, 200, 500), cores = 1, blocked_gibbs = TRUE,
                ess_burn = 10, burnin = NULL, B = 10, internal_input_idx = NULL, linked_idx = NULL, id = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  if (!is.null(internal_input_idx)) {
    lifecycle::deprecate_warn(
      when = "2.5.0",
      what = "dgp(internal_input_idx)",
      details = c(i = "The argument will be dropped in the next release.",
                  i = "To set up connections of DGPs for linked emulation, please use the updated `lgp()` function instead."
      )
    )
  }

  if (!is.null(linked_idx)) {
    lifecycle::deprecate_warn(
      when = "2.5.0",
      what = "dgp(linked_idx)",
      details = c(i = "The argument will be dropped in the next release.",
                  i = "To set up connections of DGPs for linked emulation, please use the updated `lgp()` function instead."
      )
    )
  }

  if (  is.null(nugget) ) {
    if ( all(nugget_est) ) {
      nugget = 0.01
    } else {
        if ( is.null(likelihood) ) {
          nugget = 1e-6
        } else {
          nugget = c(rep(1e-6, depth-2), 1e-4)
        }
      }
  }

  if ( !is.matrix(X)&!is.vector(X) ) stop("'X' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(Y)&!is.vector(Y) ) stop("'Y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(X) ) X <- as.matrix(X)
  if ( is.vector(Y) ) Y <- as.matrix(Y)

  if ( nrow(X)!=nrow(Y) ) stop("'X' and 'Y' have different number of data points.", call. = FALSE)

  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)

  rank_num <- pkg.env$np$linalg$matrix_rank(X)
  if (rank_num < n_dim_X) stop("The input matrix is not full rank. This indicates perfect multicollinearity and redundant information. We recommend identifying and removing redundant columns.")

  if ( !is.null(likelihood) ){
    if (likelihood!='Hetero' &  likelihood!='Poisson' & likelihood!='NegBin' & likelihood!='Categorical' ) stop("The provided 'likelihood' is not supported.", call. = FALSE)
  }

  N <- as.integer(N)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }
  B <- as.integer(B)
  ess_burn <- as.integer(ess_burn)

  M <- as.integer(M)
  if ( !is.null(ord) ) {
    ord_wrapper <- function(x) {
      return( as.integer(ord(x) - 1) )
    }
    ord_wrapper <- reticulate::py_func(ord_wrapper)
  } else {
    ord_wrapper <- NULL
  }

  if ( !is.null(burnin) ) {
    burnin <- as.integer(burnin)
  }

  linked_idx_py <- linked_idx_r_to_py(linked_idx)

  #If struc is NULL
    depth <- as.integer(depth)
    if ( depth < 2 ) stop("'depth' must >= 2. Use gp() if you want a single-layered DGP.", call. = FALSE)

    if ( is.null(likelihood) ) {
      if ( length(name)==1 ) {
        if ( name!='sexp' & name!='matern2.5' ) stop("'name' can only be either 'sexp' or 'matern2.5'.", call. = FALSE)
        name <- rep(name, depth)
      } else {
        if ( length(name)!=depth ) {
          stop(sprintf("length(name) must equal to %i.", depth), call. = FALSE)
        }
        if ( !all(name %in% c('sexp', 'matern2.5')) ) stop("'name' can only contain either 'sexp' or 'matern2.5'.", call. = FALSE)
      }

      if ( length(lengthscale)==1 ) {
        lengthscale <- rep(lengthscale, depth)
      } else {
        if ( length(lengthscale)!=depth ) {
          stop(sprintf("length(lengthscale) must equal to %i.", depth), call. = FALSE)
        }
      }

      if ( !is.null(bounds) ){
        if ( is.vector(bounds) ) {
          bounds <- matrix(rep(bounds, depth), ncol = 2,  byrow = TRUE)
        } else {
          if ( nrow(bounds)!=depth ) {
            stop(sprintf("nrow(bounds) must equal to %i.", depth), call. = FALSE)
          }
        }
      }

      if ( length(nugget_est)==1 ) {
        nugget_est <- rep(nugget_est, n_dim_Y)
      } else {
        if ( length(nugget_est)!=n_dim_Y ) {
          stop(sprintf("length(nugget_est) should equal %i.", n_dim_Y), call. = FALSE)
        }
      }

      #if ( length(nugget)==1 ) {
      #  nugget <- rep(nugget, n_dim_Y)
      #} else {
      #  if ( length(nugget)!=n_dim_Y ) {
      #    stop(sprintf("length(nugget) should equal to %i.", n_dim_Y), call. = FALSE)
      #  }
      #}

      if ( length(nugget)==1 ) {
        nugget <- rep(nugget, depth)
      } else {
        if ( length(nugget)!=depth ) {
          stop(sprintf("length(nugget) must equal %i.", depth), call. = FALSE)
        }
      }

      if ( length(scale_est)==1 ) {
        scale_est <- rep(scale_est, n_dim_Y)
      } else {
        if ( length(scale_est)!=n_dim_Y ) {
          stop(sprintf("length(scale_est) should equal %i.", n_dim_Y), call. = FALSE)
        }
      }

      if ( length(scale)==1 ) {
        scale <- rep(scale, n_dim_Y)
      } else {
        if ( length(scale)!=n_dim_Y ) {
          stop(sprintf("length(scale) should equal %i.", n_dim_Y), call. = FALSE)
        }
      }

      no_gp_layer = depth
    } else {
      if ( length(name)==1 ) {
        if ( name!='sexp' & name!='matern2.5' ) stop("'name' can only be either 'sexp' or 'matern2.5'.", call. = FALSE)
        name <- rep(name, depth-1)
      } else {
        if ( length(name)!=(depth-1) ) {
          stop(sprintf("length(name) must equal %i.", depth-1), call. = FALSE)
        }
        if ( !all(name %in% c('sexp', 'matern2.5')) ) stop("'name' can only contain either 'sexp' or 'matern2.5'.", call. = FALSE)
      }

      if ( length(lengthscale)==1 ) {
        lengthscale <- rep(lengthscale, depth-1)
      } else {
        if ( length(lengthscale)!=(depth-1) ) {
          stop(sprintf("length(lengthscale) must equal %i.", depth-1), call. = FALSE)
        }
      }

      if ( length(nugget)==1 ) {
        nugget <- rep(nugget, depth-1)
      } else {
        if ( length(nugget)!=(depth-1) ) {
          stop(sprintf("length(nugget) must equal %i.", depth-1), call. = FALSE)
        }
      }

      if ( !is.null(bounds) ){
        if ( is.vector(bounds) ) {
          bounds <- matrix(rep(bounds, depth-1), ncol = 2,  byrow = TRUE)
        } else {
          if ( nrow(bounds)!=(depth-1) ) {
            stop(sprintf("nrow(bounds) must equal %i.", depth-1), call. = FALSE)
          }
        }
      }

      if ( n_dim_Y != 1 ) {
        stop("'Y' must be a vector or a matrix with only one column when 'likelihood' is not NULL.", call. = FALSE)
      }

      if ( likelihood == 'Hetero'|likelihood == 'NegBin' ) {
        #nugget <- rep(1e-6, 2)
        nugget_est <- rep(FALSE, 2)
        if ( length(scale_est)==1 ) {
          scale_est <- rep(scale_est, 2)
        } else {
          if ( length(scale_est)!=2 ) {
            stop(sprintf("length(scale_est) should equal %i.", 2), call. = FALSE)
          }
        }

        if ( length(scale)==1 ) {
          scale <- rep(scale, 2)
        } else {
          if ( length(scale)!=2 ) {
            stop(sprintf("length(scale) should equal %i.", 2), call. = FALSE)
          }
        }
      } else if ( likelihood == 'Poisson' ) {
        #nugget <- 1e-6
        nugget_est <- FALSE
        if ( length(scale_est)!=1 ) stop(sprintf("length(scale_est) should equal %i.", 1), call. = FALSE)
        if ( length(scale)!=1 ) stop(sprintf("length(scale) should equal %i.", 1), call. = FALSE)
      } else if ( likelihood == 'Categorical' ) {
        num_class <- length(unique(Y[,1]))
        if (num_class==2) {
          nugget_est <- FALSE
          if ( length(scale_est)!=1 ) stop(sprintf("length(scale_est) should equal %i.", 1), call. = FALSE)
          if ( length(scale)!=1 ) stop(sprintf("length(scale) should equal %i.", 1), call. = FALSE)
        } else {
          nugget_est <- rep(FALSE, num_class)
          if ( length(scale_est)==1 ) {
            scale_est <- rep(scale_est, num_class)
          } else {
            if ( length(scale_est)!=num_class ) {
              stop(sprintf("length(scale_est) should equal %i.", num_class), call. = FALSE)
            }
          }

          if ( length(scale)==1 ) {
            scale <- rep(scale, num_class)
          } else {
            if ( length(scale)!=num_class ) {
              stop(sprintf("length(scale) should equal %i.", num_class), call. = FALSE)
            }
          }
        }
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
        no_kerenl <- length(nugget_est)
        } else {
          no_kerenl <- node
        }

      if ( is.null(bounds) ){
        bds <- NULL
      } else {
        bds <- reticulate::np_array(bounds[l,])
      }

      for ( k in 1:no_kerenl ) {
        if ( l == no_gp_layer ) {
          if ( no_gp_layer == 1 ) {
            if ( isTRUE(share) ){
              length_scale <- lengthscale[l]
            } else {
              length_scale <- rep(lengthscale[l], n_dim_X)
            }
            layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), bds = bds, name = name[l], prior_name = prior, scale = scale[k], scale_est = scale_est[k], nugget = nugget[l], nugget_est = nugget_est[k],
                                                 input_dim = internal_input_idx, connect = external_input_idx)
          } else {
            if ( connect ) {
              if ( isTRUE(share) ){
                length_scale <- lengthscale[l]
              } else {
                length_scale <- rep(lengthscale[l], n_dim_X+node)
              }
              connect_info <- as.integer(1:n_dim_X - 1)
              if (!is.null(internal_input_idx)){
                idx <- reticulate::py_to_r(internal_input_idx)+1
                connect_info <- c(connect_info[idx],connect_info[-idx])
              }
              layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), bds = bds, name = name[l], prior_name = prior, scale = scale[k], scale_est = scale_est[k], nugget = nugget[l], nugget_est = nugget_est[k],
                                                   connect = reticulate::np_array(connect_info) )
            } else {
              if ( isTRUE(share) ){
                length_scale <- lengthscale[l]
              } else {
                length_scale <- rep(lengthscale[l], node)
              }
              layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), bds = bds, name = name[l], prior_name = prior, scale = scale[k], scale_est = scale_est[k], nugget = nugget[l], nugget_est = nugget_est[k])
            }
          }
        } else {
            if ( l == 1 ) {
              if ( isTRUE(share) ){
                length_scale <- lengthscale[l]
              } else {
                length_scale <- rep(lengthscale[l], n_dim_X)
              }
              layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), bds = bds, name = name[l], prior_name = prior, nugget = nugget[l],
                                                   input_dim = internal_input_idx, connect = external_input_idx)
            } else {
              if ( connect ) {
                if ( isTRUE(share) ){
                  length_scale <- lengthscale[l]
                } else {
                  length_scale <- rep(lengthscale[l], n_dim_X+node)
                }
                connect_info <- as.integer(1:n_dim_X - 1)
                if (!is.null(internal_input_idx)){
                  idx <- reticulate::py_to_r(internal_input_idx)+1
                  connect_info <- c(connect_info[idx],connect_info[-idx])
                }
                layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), bds = bds, name = name[l], prior_name = prior, nugget = nugget[l],
                                                     connect = reticulate::np_array(connect_info))
              } else {
                if ( isTRUE(share) ){
                  length_scale <- lengthscale[l]
                } else {
                  length_scale <- rep(lengthscale[l], node)
                }
                layer_l[[k]] <- pkg.env$dgpsi$kernel(length = reticulate::np_array(length_scale), bds = bds, name = name[l], prior_name = prior, nugget = nugget[l])
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
      } else if ( likelihood == 'Categorical' ) {
        struc[[depth]] <- c(pkg.env$dgpsi$Categorical(num_class))
      }
    }
    if ( isTRUE(verb) ) {
      message(" done")
      Sys.sleep(0.5)
    }


  if ( isTRUE(verb) ) message("Initializing the DGP emulator ...", appendLF = FALSE)

  obj <- pkg.env$dgpsi$dgp(X, Y, struc, check_rep, blocked_gibbs, vecchia, M, ord_wrapper)

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
    if ( identical(cores,as.integer(1)) ){
      obj$train(N, ess_burn, disable)
    } else {
      obj$ptrain(N, ess_burn, disable, cores)
    }
    est_obj <- obj$estimate(burnin)
  } else {
    est_obj <- obj$estimate(NULL)
  }

  if ( isTRUE(verb) ) message("Imputing ...", appendLF = FALSE)

  res <- list()
  res[['id']] <- if (is.null(id)) uuid::UUIDgenerate() else id
  res[['data']][['X']] <- unname(X)
  res[['data']][['Y']] <- unname(Y)
  res[['specs']] <- extract_specs(est_obj, "dgp")
  res[['specs']][['internal_dims']] <- if( is.null(internal_input_idx) ) 1:n_dim_X else as.integer(reticulate::py_to_r(internal_input_idx)+1)
  res[['specs']][['external_dims']] <- if( is.null(internal_input_idx) ) FALSE else as.integer(reticulate::py_to_r(external_input_idx)+1)
  res[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx
  res[['specs']][['vecchia']] <- vecchia
  res[['specs']][['M']] <- M
  res[['constructor_obj']] <- obj
  seed <- sample.int(100000, 1)
  set_seed(seed)
  res[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = blocked_gibbs)
  res[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_py, block = blocked_gibbs)
  res[['specs']][['seed']] <- seed
  res[['specs']][['B']] <- B

  class(res) <- "dgp"
  if ( isTRUE(verb) ) message(" done")
  return(res)
}


#' @title Continue training a DGP emulator
#'
#' @description This function implements additional training iterations for a DGP emulator.
#'
#' @param object an instance of the `dgp` class.
#' @param N additional number of iterations to train the DGP emulator. If set to `NULL`, the number of iterations is set to `500` if the DGP emulator
#'     was constructed without the Vecchia approximation, and is set to `200` if Vecchia approximation was used. Defaults to `NULL`.
#' @param cores the number of processes to be used to optimize GP components (in the same layer) at each M-step of the training. If set to `NULL`,
#'     the number of processes is set to `(max physical cores available - 1)` if the DGP emulator was constructed without the Vecchia approximation.
#'     Otherwise, the number of processes is set to `max physical cores available %/% 2`. Only use multiple processes when there is a large number of
#'     GP components in different layers and optimization of GP components is computationally expensive. Defaults to `1`.
#' @param ess_burn number of burnin steps for ESS-within-Gibbs
#'     at each I-step of the training. Defaults to `10`.
#' @param verb a bool indicating if a progress bar will be printed during training. Defaults to `TRUE`.
#' @param burnin the number of training iterations to be discarded for
#'     point estimates calculation. Must be smaller than the overall training iterations
#'     so-far implemented. If this is not specified, only the last 25% of iterations
#'     are used. This overrides the value of `burnin` set in [dgp()]. Defaults to `NULL`.
#' @param B the number of imputations to produce predictions. Increase the value to account for
#'     more imputation uncertainty. This overrides the value of `B` set in [dgp()] if `B` is not
#'     `NULL`. Defaults to `NULL`.
#'
#' @return An updated `object`.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @note
#' * One can also use this function to fit an untrained DGP emulator constructed by [dgp()] with `training = FALSE`.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()]
#'   in `object` will be removed and not contained in the returned object.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export

continue <- function(object, N = NULL, cores = 1, ess_burn = 10, verb = TRUE, burnin = NULL, B = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }

  if( is.null(N) ) {
    if (object[['specs']][['vecchia']]) {
      N <- 200
    } else {
      N <- 500
    }
  }
  N <- as.integer(N)

  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  if ( is.null(B) ){
    B <- as.integer(length(object$emulator_obj$all_layer_set))
  } else {
    B <- as.integer(B)
  }
  ess_burn <- as.integer(ess_burn)

  if( !is.null(burnin) ) {
    burnin <- as.integer(burnin)
  }

  if ( verb ) {
    disable <- FALSE
    message("Continue DGP training:")
  } else {
    disable <- TRUE
  }

  linked_idx <- object$container_obj$local_input_idx
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  isblock <- constructor_obj_cp$block
  if ( identical(cores,as.integer(1)) ){
    constructor_obj_cp$train(N, ess_burn, disable)
  } else {
    constructor_obj_cp$ptrain(N, ess_burn, disable, cores)
  }
  est_obj <- constructor_obj_cp$estimate(burnin)

  if ( isTRUE(verb) ) message("Imputing ...", appendLF = FALSE)
  new_object <- list()
  new_object[['id']] <- object$id
  new_object[['data']][['X']] <- object$data$X
  new_object[['data']][['Y']] <- object$data$Y
  new_object[['specs']] <- extract_specs(est_obj, "dgp")
  new_object[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
  new_object[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
  new_object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
  new_object[['specs']][['vecchia']] <- object[['specs']][['vecchia']]
  new_object[['specs']][['M']] <- object[['specs']][['M']]
  new_object[['constructor_obj']] <- constructor_obj_cp
  seed <- sample.int(100000, 1)
  set_seed(seed)
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, isblock)
  new_object[['specs']][['seed']] <- seed
  new_object[['specs']][['B']] <- B
  if ( "design" %in% names(object) ) new_object[['design']] <- object$design
  class(new_object) <- "dgp"
  if ( isTRUE(verb) ) message(" done")
  pkg.env$py_gc$collect()
  gc(full=T)
  return(new_object)
}

linked_idx_r_to_py <- function(linked_idx){
  if( !is.null(linked_idx) ) {
    if ( !is.list(linked_idx) ) {
      linked_idx <- reticulate::np_array(as.integer(linked_idx - 1))
    } else {
      for ( i in 1:length(linked_idx) ){
        if ( !is.null(linked_idx[[i]]) ) linked_idx[[i]] <- reticulate::np_array(as.integer(linked_idx[[i]] - 1))
      }
    }
  }
  return(linked_idx)
}

linked_idx_py_to_r <- function(linked_idx){
  if( !is.null(linked_idx) ) {
    if ( !is.list(linked_idx) ) {
      linked_idx <- as.integer(linked_idx + 1)
    } else {
      for ( i in 1:length(linked_idx) ){
        if ( !is.null(linked_idx[[i]]) ) linked_idx[[i]] <- as.integer(linked_idx[[i]] + 1)
      }
    }
  }
  return(linked_idx)
}
