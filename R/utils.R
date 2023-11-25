#' @title Combine layers
#'
#' @description This function combines customized layers into a DGP or linked (D)GP structure.
#'
#' @param ... a sequence of lists:
#' 1. For DGP emulations, each list represents a DGP layer and contains GP nodes (produced by [kernel()]), or
#'    likelihood nodes (produced by [Poisson()], [Hetero()], or [NegBin()]).
#' 2. For linked (D)GP emulations, each list represents a system layer and contains emulators (produced by [gp()] or
#'    [dgp()]) in that layer.
#'
#' @return A list defining a DGP structure (for `struc` of [dgp()]) or a linked (D)GP structure
#'     (for `struc` for [lgp()]).
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See lgp() for an example.
#' }
#' @md
#' @export
combine <- function(...) {
  res = list(...)
  return(res)
}

#' @title Pack GP and DGP emulators into a bundle
#'
#' @description This function packs GP emulators and DGP emulators into a `bundle` class for
#'     sequential designs if each emulator emulates one output dimension of the underlying simulator.
#'
#' @param ... a sequence or a list of emulators produced by [gp()] or [dgp()].
#' @param id an ID to be assigned to the bundle emulator. If an ID is not provided (i.e., `id = NULL`), a UUID (Universally Unique Identifier) will be automatically generated
#'    and assigned to the emulator. Default to `NULL`.
#'
#' @return An S3 class named `bundle` to be used by [design()] for sequential designs. It has:
#' - a slot called `id` that is assigned through the `id` argument.
#' - *N* slots named `emulator1,...,emulatorN`, each of which contains a GP or DGP emulator, where *N* is the number of emulators
#'   that are provided to the function.
#' - a slot called `data` which contains two elements `X` and `Y`. `X` contains *N* matrices named `emulator1,...,emulatorN` that are
#'   training input data for different emulators. `Y` contains *N* single-column matrices named `emulator1,...,emulatorN` that are
#'   training output data for different emulators.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # load packages and the Python env
#' library(lhs)
#' library(dgpsi)
#'
#' # construct a function with a two-dimensional output
#' f <- function(x) {
#'  y1 = sin(30*((2*x-1)/2-0.4)^5)*cos(20*((2*x-1)/2-0.4))
#'  y2 = 1/3*sin(2*(2*x - 1))+2/3*exp(-30*(2*(2*x-1))^2)+1/3
#'  return(cbind(y1,y2))
#' }
#'
#' # generate the initial design
#' X <- maximinLHS(10,1)
#' Y <- f(X)
#'
#' # generate the validation data
#' validate_x <- maximinLHS(30,1)
#' validate_y <- f(validate_x)
#'
#' # training a 2-layered DGP emulator with respect to each output with the global connection off
#' m1 <- dgp(X, Y[,1], connect=F)
#' m2 <- dgp(X, Y[,2], connect=F)
#'
#' # specify the range of the input dimension
#' lim <- c(0, 1)
#'
#' # pack emulators to form an emulator bundle
#' m <- pack(m1, m2)
#'
#' # 1st wave of the sequential design with 10 steps with target RMSE 0.01
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y, target = 0.01)
#'
#' # 2rd wave of the sequential design with 10 steps, the same target, and the aggregation
#' # function that takes the average of the criterion scores across the two outputs
#' g <- function(x){
#'   return(rowMeans(x))
#' }
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x,
#'                     y_test = validate_y, aggregate = g, target = 0.01)
#'
#' # draw sequential designs of the two packed emulators
#' draw(m, emulator = 1, type = 'design')
#' draw(m, emulator = 2, type = 'design')
#'
#' # inspect the traces of RMSEs of the two packed emulators
#' draw(m, emulator = 1, type = 'rmse')
#' draw(m, emulator = 2, type = 'rmse')
#'
#' # write and read the constructed emulator bundle
#' write(m, 'bundle_dgp')
#' m <- read('bundle_dgp')
#'
#' # unpack the bundle into individual emulators
#' m_unpacked <- unpack(m)
#'
#' # plot OOS validations of individual emulators
#' plot(m_unpacked[[1]], x_test = validate_x, y_test = validate_y[,1])
#' plot(m_unpacked[[2]], x_test = validate_x, y_test = validate_y[,2])
#' }
#' @md
#' @export
pack <- function(..., id = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  res <- list(...)
  if (length(res) == 1 && is.list(res[[1]])) {
    res <- res[[1]]
  }

  X_all <- list()
  Y_all <- list()
  if ( length(res)==1 ) stop("The function needs at least two emulators to pack.", call. = FALSE)
  training_input <- res[[1]]$data$X
  training_output <- c()
  for ( i in 1:length(res) ){
    if ( !inherits(res[[i]],"gp") & !inherits(res[[i]],"dgp") ) stop("The function only accepts GP or DGP emulators as inputs.", call. = FALSE)
    #if ( inherits(res[[i]],"dgp") ){
    #  if ( res[[i]]$constructor_obj$all_layer[[res[[i]]$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
    #    stop("The function can only pack DGP emulators without likelihood layers.", call. = FALSE)
    #  }
    #}
    if ( !identical(res[[i]]$data$X, training_input) ) stop("The function can only pack emulators with common training input data.", call. = FALSE)
    Y_dim <- ncol(res[[i]]$data$Y)
    if ( Y_dim!=1 ) stop(sprintf("The function is only applicable to emulators with 1D output. Your emulator %i has %i output dimensions.", i, Y_dim), call. = FALSE)
    X_all[[paste('emulator', i ,sep="")]] <- unname(training_input)
    Y_all[[paste('emulator', i ,sep="")]] <- unname(res[[i]]$data$Y)
    names(res)[i] <- paste('emulator', i, sep="")
  }
  res[['id']] <- if (is.null(id)) uuid::UUIDgenerate() else id
  res[['data']][['X']] <- X_all
  res[['data']][['Y']] <- Y_all
  class(res) <- "bundle"
  #pkg.env$py_gc$collect()
  #gc(full=T)
  return(res)
}


#' @title Unpack a bundle of (D)GP emulators
#'
#' @description This function unpacks a bundle of (D)GP emulators safely so any further manipulations of unpacked individual emulators
#'     will not impact the ones in the bundle.
#'
#' @param object an instance of the class `bundle`.
#'
#' @return A named list that contains individual emulators (named `emulator1,...,emulatorS`) packed in `object`,
#'    where `S` is the number of emulators in `object`.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See pack() for an example.
#' }
#' @md
#' @export
unpack <- function(object) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"bundle") ){
    stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)
  }
  n_emulators <- length(object) - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
  res <- list()
  for ( i in 1:n_emulators ){
    res[[paste('emulator', i, sep="")]] <- object[[paste('emulator', i, sep="")]]
    res[[paste('emulator', i, sep="")]]$constructor_obj <- pkg.env$copy$deepcopy(res[[paste('emulator', i, sep="")]]$constructor_obj)
    res[[paste('emulator', i, sep="")]]$container_obj <- pkg.env$copy$deepcopy(res[[paste('emulator', i, sep="")]]$container_obj)
    res[[paste('emulator', i, sep="")]]$emulator_obj <- pkg.env$copy$deepcopy(res[[paste('emulator', i, sep="")]]$emulator_obj)
  }
  pkg.env$py_gc$collect()
  gc(full=T)
  return(res)
}


#' @title Save the constructed emulator
#'
#' @description This function saves the constructed emulator to a `.pkl` file.
#'
#' @param object an instance of the S3 class `gp`, `dgp`, `lgp`, or `bundle`.
#' @param pkl_file the path to and the name of the `.pkl` file to which
#'     the emulator `object` is saved.
#' @param light a bool indicating if a light version of the constructed emulator (that requires a small storage) will be saved.
#'     This argument has no effects on GP or bundles of GP emulators. Defaults to `TRUE`.
#'
#' @return No return value. `object` will be save to a local `.pkl` file specified by `pkl_file`.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note Since the constructed emulators are 'python' objects, [save()] from R will not work as it is only for R objects.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), lgp(), or pack() for an example.
#' }
#' @md
#' @export
write <- function(object, pkl_file, light = TRUE) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  pkl_file <- tools::file_path_sans_ext(pkl_file)
  if (light) {
    if (inherits(object,"dgp")){
      if ( !"seed" %in% names(object$specs) ) stop("The supplied 'object' cannot be saved in light mode. To save, either set 'light = FALSE' or produce a new version of 'object' by set_imp().", call. = FALSE)
      object[['emulator_obj']] <- NULL
      object[['container_obj']] <- NULL
    } else if (inherits(object,"lgp")){
      if ( !"seed" %in% names(object$specs) ) stop("The supplied 'object' cannot be saved in light mode. To save, either set 'light = FALSE' or re-construct the 'object' by lgp().", call. = FALSE)
      object[['emulator_obj']] <- NULL
    } else if (inherits(object,"bundle")){
      N <- length(object) - 1
      if ( "id" %in% names(object) ) N <- N - 1
      if ( "design" %in% names(object) ) N <- N - 1
      for ( i in 1:N ){
        if ( inherits(object[[paste('emulator',i, sep='')]],"dgp") ) {
          if ( !"seed" %in% names(object[[paste('emulator',i, sep='')]][['specs']]) ) stop("The supplied 'object' cannot be saved in light mode. To save, either set 'light = FALSE' or produce a new version of 'object' by updating the included DGP emulators via set_imp().", call. = FALSE)
          object[[paste('emulator',i, sep='')]][['emulator_obj']] <- NULL
          object[[paste('emulator',i, sep='')]][['container_obj']] <- NULL
        }
      }
    }
  }
  label <- class(object)
  lst <- unclass(object)
  lst[['label']] <- label
  pkg.env$dgpsi$write(lst, pkl_file)
}

#' @title Random seed generator
#'
#' @description This function initializes a random number generator that sets the random seed in both R and Python
#'    to ensure reproducible results from the package.
#'
#' @param seed a single integer value.
#'
#' @return No return value.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export
set_seed <- function(seed) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  seed <- as.integer(seed)
  set.seed(seed)
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  pkg.env$dgpsi$nb_seed(seed)
}

#' @title Load the stored emulator
#'
#' @description This function loads the `.pkl` file that stores the emulator.
#'
#' @param pkl_file the path to and the name of the `.pkl` file where the emulator is stored.
#'
#' @return The S3 class of a GP emulator, a DGP emulator, a linked (D)GP emulator, or a bundle of (D)GP emulators.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), lgp(), or pack() for an example.
#' }
#' @md
#' @export
read <- function(pkl_file) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  pkl_file <- tools::file_path_sans_ext(pkl_file)
  res <- pkg.env$dgpsi$read(pkl_file)
  if ('label' %in% names(res)){
    label <- res$label
    res$label <- NULL
    if (label == "gp"){
      if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
      class(res) <- "gp"
    } else if (label == "dgp"){
      if ('emulator_obj' %in% names(res)){
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "dgp"
      } else {
        burnin <- res$constructor_obj$burnin
        est_obj <- res$constructor_obj$estimate(burnin)
        B <- res$specs$B
        isblock <- res$constructor_obj$block
        linked_idx <- if ( isFALSE( res[['specs']][['linked_idx']]) ) {NULL} else {res[['specs']][['linked_idx']]}
        set_seed(res$specs$seed)
        res[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
        res[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_r_to_py(linked_idx), isblock)
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "dgp"
      }
    } else if (label == "lgp"){
      if ('emulator_obj' %in% names(res)){
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "lgp"
      } else {
        B <- res$specs$B
        extracted_struc <- res$constructor_obj
        set_seed(res$specs$seed)
        obj <- pkg.env$dgpsi$lgp(all_layer = pkg.env$copy$deepcopy(extracted_struc), N = B)
        res[['emulator_obj']] <- obj
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "lgp"
      }
    } else if (label == 'bundle'){
      if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
      N <- length(res) - 2
      if ( "design" %in% names(res) ) N <- N - 1
      class(res) <- "bundle"
      for ( i in 1:N ){
        if ('emulator_obj' %in% names(res[[paste('emulator',i, sep='')]])) {
          type <- pkg.env$py_buildin$type(res[[paste('emulator',i, sep='')]]$emulator_obj)$'__name__'
          if ( type=='emulator' ) {
            class(res[[paste('emulator',i, sep='')]]) <- "dgp"
          } else if ( type=='gp' ) {
            class(res[[paste('emulator',i, sep='')]]) <- "gp"
          }
        } else {
          burnin <- res[[paste('emulator',i, sep='')]]$constructor_obj$burnin
          est_obj <- res[[paste('emulator',i, sep='')]]$constructor_obj$estimate(burnin)
          B <- res[[paste('emulator',i, sep='')]]$specs$B
          isblock <- res[[paste('emulator',i, sep='')]]$constructor_obj$block
          linked_idx <- if ( isFALSE( res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]) ) {NULL} else {res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]}
          set_seed(res[[paste('emulator',i, sep='')]]$specs$seed)
          res[[paste('emulator',i, sep='')]][['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
          res[[paste('emulator',i, sep='')]][['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_r_to_py(linked_idx), isblock)
          class(res[[paste('emulator',i, sep='')]]) <- "dgp"
        }
        if (!'id' %in% names(res[[paste('emulator',i, sep='')]])) res[[paste('emulator',i, sep='')]][['id']] <- uuid::UUIDgenerate()
      }
    }
  } else {
    if ('emulator_obj' %in% names(res)){
      type <- pkg.env$py_buildin$type(res$emulator_obj)$'__name__'
      if ( type=='emulator' ) {
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "dgp"
      } else if ( type=='gp' ) {
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "gp"
      } else if ( type=='lgp' ) {
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "lgp"
      }
    } else {
      N <- length(res) - 1
      if ( "design" %in% names(res) ) N <- N - 1
      if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
      class(res) <- "bundle"
      for ( i in 1:N ){
        type <- pkg.env$py_buildin$type(res[[paste('emulator',i, sep='')]]$emulator_obj)$'__name__'
        if ( type=='emulator' ) {
          if (!'id' %in% names(res[[paste('emulator',i, sep='')]])) res[[paste('emulator',i, sep='')]][['id']] <- uuid::UUIDgenerate()
          class(res[[paste('emulator',i, sep='')]]) <- "dgp"
        } else if ( type=='gp' ) {
          if (!'id' %in% names(res[[paste('emulator',i, sep='')]])) res[[paste('emulator',i, sep='')]][['id']] <- uuid::UUIDgenerate()
          class(res[[paste('emulator',i, sep='')]]) <- "gp"
        }
      }
    }
  }
  return(res)
}


#' @title Summary of a constructed GP, DGP, or linked (D)GP emulator
#'
#' @description This function summarizes key information of a GP, DGP or linked (D)GP emulator.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param ... N/A.
#'
#' @return A table summarizing key information contained in `object`.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name summary
NULL

#' @rdname summary
#' @method summary gp
#' @export
summary.gp <- function(object, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
}

#' @rdname summary
#' @method summary dgp
#' @export
summary.dgp <- function(object, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
}

#' @rdname summary
#' @method summary lgp
#' @export
summary.lgp <- function(object, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
}


#' @title Set linked indices
#'
#' @description This function adds the linked information to a GP or DGP emulator if the information is not provided
#'     when the emulator is constructed by [gp()] or [dgp()].
#'
#' @param object an instance of the S3 class `gp` or `dgp`.
#' @param idx same as the argument `linked_idx` of [gp()] and [dgp()].
#'
#' @return An updated `object` with the information of `idx` incorporated.
#'
#' @note This function is useful when different models are emulated by different teams. Each team can create their (D)GP emulator
#'     even without knowing how different emulators are connected together. When this information is available and
#'     different emulators are collected, the connection information between emulators can then be assigned to
#'     individual emulators with this function.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See lgp() for an example.
#' }
#' @md
#' @export
set_linked_idx <- function(object, idx) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  idx_py <- linked_idx_r_to_py(idx)
  object[['container_obj']] <- object$container_obj$set_local_input(idx_py, TRUE)
  object[['specs']][['linked_idx']] <- if ( is.null(idx) ) FALSE else idx
  return(object)
}

#' @title Reset number of imputations for a DGP emulator
#'
#' @description This function resets the number of imputations for predictions from a DGP emulator.
#'
#' @param object an instance of the S3 class `dgp`.
#' @param B the number of imputations to produce predictions from `object`. Increase the value to account for
#'     more imputation uncertainties with slower predictions. Decrease the value for lower imputation uncertainties
#'     but faster predictions. Defaults to `5`.
#'
#' @return An updated `object` with the information of `B` incorporated.
#'
#' @note
#' * This function is useful when a DGP emulator has been trained and one wants to make faster predictions by decreasing
#'    the number of imputations without rebuilding the emulator.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()]
#'   in `object` will be removed and not contained in the returned object.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See design() for an example.
#' }
#' @md
#' @export
set_imp <- function(object, B = 5) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  B <- as.integer(B)

  linked_idx <- object$container_obj$local_input_idx
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  burnin <- constructor_obj_cp$burnin
  isblock <- constructor_obj_cp$block
  est_obj <- constructor_obj_cp$estimate(burnin)

  new_object <- list()
  new_object[['id']] <- object$id
  new_object[['data']][['X']] <- object$data$X
  new_object[['data']][['Y']] <- object$data$Y
  new_object[['specs']] <- extract_specs(est_obj, "dgp")
  if ("internal_dims" %in% names(object[['specs']])){
    new_object[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
    new_object[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
  }
  new_object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
  new_object[['constructor_obj']] <- constructor_obj_cp
  id <- sample.int(100000, 1)
  set_seed(id)
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, isblock)
  new_object[['specs']][['seed']] <- id
  new_object[['specs']][['B']] <- B
  if ( "design" %in% names(object) ) new_object[['design']] <- object$design
  class(new_object) <- "dgp"
  pkg.env$py_gc$collect()
  gc(full=T)
  return(new_object)
}


#' @title Trim the sequences of model parameters of a DGP emulator
#'
#' @description This function trim the sequences of model parameters of a DGP emulator
#'     that are generated during the training.
#'
#' @param object an instance of the S3 class `dgp`.
#' @param start the first iteration before which all iterations are trimmed from the sequences.
#' @param end the last iteration after which all iterations are trimmed from the sequences.
#'     Set to `NULL` to keep all iterations after (including) `start`. Defaults to `NULL`.
#' @param thin the interval between the `start` and `end` iterations to thin out the sequences.
#'     Defaults to 1.
#'
#' @return An updated `object` with trimmed sequences of model parameters.
#'
#' @note
#' * This function is useful when a DGP emulator has been trained and one wants to trim
#'   the sequences of model parameters and use the trimmed sequences to generate the point estimates
#'   of DGP model parameters for predictions.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()]
#'   in `object` will be removed and not contained in the returned object.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export
window <- function(object, start, end = NULL, thin = 1) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  if (length(start) != 1) stop("'start' must be an integer.", call. = FALSE)
  if ( !is.null(end) ){
    if (length(end) != 1) stop("'end' must be an integer.", call. = FALSE)
  }
  if (length(thin) != 1) stop("'thin' must be an integer.", call. = FALSE)
  if ( !is.null(end) ){
    if (start > end) stop("'start' must be before 'end'", call. = FALSE)
  }

  linked_idx <- object$container_obj$local_input_idx
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  isblock <- constructor_obj_cp$block
  niter <- constructor_obj_cp$N + 1
  if ( is.null(end) ) end <- niter
  if (end > niter) end <- niter
  idx <- start:end
  idx <- idx[idx %% thin == 0]
  constructor_obj_cp$N <- as.integer(length(idx) - 1)
  for ( l in 1:constructor_obj_cp$n_layer ){
    n_kernel <- length(constructor_obj_cp$all_layer[[l]])
    for (k in 1:n_kernel){
      if (constructor_obj_cp$all_layer[[l]][[k]]$type == 'gp'){
        final_est <- constructor_obj_cp$all_layer[[l]][[k]]$para_path[idx[length(idx)],]
        constructor_obj_cp$all_layer[[l]][[k]]$scale <- reticulate::np_array(final_est[1])
        constructor_obj_cp$all_layer[[l]][[k]]$length <- reticulate::np_array(final_est[2:(length(final_est)-1)])
        constructor_obj_cp$all_layer[[l]][[k]]$nugget <- reticulate::np_array(final_est[length(final_est)])
        constructor_obj_cp$all_layer[[l]][[k]]$para_path <- constructor_obj_cp$all_layer[[l]][[k]]$para_path[idx,,drop=FALSE]
      }
    }
  }

  burnin <- 0L
  est_obj <- constructor_obj_cp$estimate(burnin)
  B <- length(object$emulator_obj$all_layer_set)

  new_object <- list()
  new_object[['id']] <- object$id
  new_object[['data']][['X']] <- object$data$X
  new_object[['data']][['Y']] <- object$data$Y
  new_object[['specs']] <- extract_specs(est_obj, "dgp")
  if ("internal_dims" %in% names(object[['specs']])){
    new_object[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
    new_object[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
  }
  new_object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
  new_object[['constructor_obj']] <- constructor_obj_cp
  id <- sample.int(100000, 1)
  set_seed(id)
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, isblock)
  new_object[['specs']][['seed']] <- id
  new_object[['specs']][['B']] <- B
  if ( "design" %in% names(object) ) new_object[['design']] <- object$design
  class(new_object) <- "dgp"
  pkg.env$py_gc$collect()
  gc(full=T)
  return(new_object)
}


#' @title Calculate negative predicted log-likelihood
#'
#' @description This function computes the negative predicted log-likelihood from a
#'     DGP emulator with a likelihood layer.
#'
#' @param object an instance of the `dgp` class and it should be produced by [dgp()] with one of the following two settings:
#' 1. if `struc = NULL`, `likelihood` is not `NULL`;
#' 2. if a customized structure is provided to `struc`, the final layer must be likelihood layer containing only one
#'    likelihood node produced by [Poisson()], [Hetero()], or [NegBin()].
#' @param x a matrix where each row is an input testing data point and each column is an input dimension.
#' @param y a matrix with only one column where each row is a scalar-valued testing output data point.
#'
#' @return An updated `object` is returned with an additional slot named `NLL` that contains two elements.
#'     The first one, named `meanNLL`, is a scalar that gives the average negative predicted log-likelihood
#'     across all testing data points. The second one, named `allNLL`, is a vector that gives the negative predicted
#'     log-likelihood for each testing data point.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note Any R vector detected in `x` and `y` will be treated as a column vector and automatically converted into a single-column
#'     R matrix. Thus, if `x` is a single testing data point with multiple dimensions, it must be given as a matrix.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to compute the negative predicted log-likelihood
#' # using nllik().
#' }
#' @md
#' @export
nllik <- function(object, x, y) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(y)&!is.vector(y) ) stop("'y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) x <- as.matrix(x)
  if ( is.vector(y) ) y <- as.matrix(y)

  if ( nrow(x)!=nrow(y) ) stop("'x' and 'y' have different number of data points.", call. = FALSE)
  n_dim_y <- ncol(y)
  if ( n_dim_y != 1 ) {
    stop("'y' must be a vector or a matrix with only one column.", call. = FALSE)
  }
  res <- object$emulator_obj$nllik(x, y)
  named_res <- list("meanNLL" = res[[1]], "allNLL" = res[[2]])
  object$NLL <- named_res
  return(object)
}


#' @title Plot of DGP model parameter traces
#'
#' @description This function plots the traces of model parameters of a chosen GP node
#'     in a DGP emulator.
#'
#' @param object an instance of the `dgp` class.
#' @param layer the index of a layer. Defaults to `NULL` for the final layer.
#' @param node the index of a GP node in the layer specified by `layer`. Defaults to `1` for the first GP node in the
#'     corresponding layer.
#'
#' @return A `ggplot` object.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export

trace_plot <- function(object, layer = NULL, node = 1) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)

  all_layer <- object$constructor_obj$all_layer

  layer_no <- length(all_layer)

  if ( is.null(layer) ){
    layer = layer_no
  }

  layer_l <- all_layer[[layer]]
  kernel_no <- length(layer_l)
  kernel <- layer_l[[node]]
  if (kernel$type == 'gp'){
    path <- kernel$para_path
    n_para <- ncol(path)
    dat <- as.data.frame(path)
    for ( i in 1:n_para){
      if ( i==1 ){
        colnames(dat)[i] <- 'Variance'
      } else if ( i==n_para ){
        colnames(dat)[i] <- 'Nugget'
      } else {
        colnames(dat)[i] <- paste('Lengthscale', i-1, sep = ' ')
      }
    }
    dat$Iteration <-  seq(nrow(path))
    mm <- reshape2::melt(dat, id.var="Iteration")

    color <- c("#E69F00", rep("#56B4E9",n_para-2), "#009E73")
    ggplot2::ggplot(mm, ggplot2::aes_(x = ~Iteration, y = ~value)) + ggplot2::geom_line(ggplot2::aes_(color = ~variable)) +
      ggplot2::facet_grid(variable ~ ., scales = "free_y") + ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title=paste("Node ", node, " (of ", kernel_no, ")", " in Layer ", layer, " (of ", layer_no, ")", sep = ""), x ="Iteration", y = "Parameter value") +
      ggplot2::scale_color_manual(values = color)
  } else {
    message('There is nothing to plot for a likelihood node, please choose a GP node instead.')
  }

}

#' @title Static pruning of a DGP emulator
#'
#' @description This function implements the static pruning of a DGP emulator.
#'
#' @param object an instance of the `dgp` class that is generated by `dgp()` with `struc = NULL`.
#' @param control a list that can supply the following two components to control the static pruning of the DGP emulator:
#' * `min_size`, the minimum number of design points required to trigger the pruning. Defaults to 10 times of the input dimensions.
#' * `threshold`, the R2 value above which a GP node is considered redundant and removable. Defaults to `0.97`.
#' @param verb a bool indicating if the trace information will be printed during the function execution. Defaults to `TRUE`.
#'
#' @return An updated `object` that could be an instance of `gp`, `dgp`, or `bundle` (of GP emulators) class.
#'
#' @note
#' * The function requires a DGP emulator that has been trained with a dataset comprising a minimum size equal to `min_size` in `control`.
#'    If the training dataset size is smaller than this, it is suggested to enrich the design of the DGP emulator and prune its
#'    structure dynamically using the `design()` function. Depending on the design of the DGP emulator, the static pruning may not be accurate.
#'    It is thus suggested to implement dynamic pruning as a part of the sequential design via `design()`.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()];
#'
#'   in `object` will be removed and not contained in the returned object.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # load the package and the Python env
#' library(dgpsi)
#'
#' # construct the borehole function over a hypercube
#' f <- function(x){
#'   x[,1] <- (0.15 - 0.5) * x[,1] + 0.5
#'   x[,2] <- exp((log(50000) - log(100)) * x[,2] + log(100))
#'   x[,3] <- (115600 - 63070) *x[,3] + 63070
#'   x[,4] <- (1110 - 990) * x[,4] + 990
#'   x[,5] <- (116 - 63.1) * x[,5] + 63.1
#'   x[,6] <- (820 - 700) * x[,6] + 700
#'   x[,7] <- (1680 - 1120) * x[,7] + 1120
#'   x[,8] <- (12045 - 9855) * x[,8] + 9855
#'   y <- apply(x, 1, RobustGaSP::borehole)
#' }
#'
#' # set a random seed
#' set_seed(999)
#'
#' # generate training data
#' X <- maximinLHS(80, 8)
#' Y <- f(X)
#'
#' # generate validation data
#' validate_x <- maximinLHS(500, 8)
#' validate_y <- f(validate_x)
#'
#' # training a DGP emulator with anisotropic squared exponential kernels
#' m <- dgp(X, Y, share = F)
#'
#' # OOS validation of the DGP emulator
#' plot(m, validate_x, validate_y)
#'
#' # prune the emulator until no more GP nodes are removable
#' m <- prune(m)
#'
#' # OOS validation of the resulting emulator
#' plot(m, validate_x, validate_y)
#' }
#' @md
#' @export

prune <- function(object, control = list(), verb = TRUE) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }

  n_dim_X <- ncol(object$data$X)
  con <- list(min_size = 10*n_dim_X, threshold = 0.97)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in 'control': ", paste(noNms,collapse=", "), call. = FALSE, immediate. = TRUE)
  }

  if (con$min_size<1) stop("'min_size' in 'control' must be at least 1.", call. = FALSE)
  if (con$threshold>1 || con$threshold<0) stop("'threshold' in 'control' must be between 0 and 1.", call. = FALSE)

  if (nrow(object$data$X) < con$min_size) {
    stop("To prune, 'object' needs to be trained with a dataset comprising a size at least equal to 'min_size' in 'control'. Use design() to enrich the design size.", call. = FALSE)
  }

  if (!"internal_dims" %in% names(object[['specs']])) {
    stop("'object' must be an instance of the 'dgp' class generated by dgp() with 'struc = NULL'.", call. = FALSE)
  } else {
    n_layer <- object$constructor_obj$n_layer
    if (object$constructor_obj$all_layer[[n_layer]][[1]]$type!='gp') {
      n_layer <- n_layer - 1
      if (n_layer == 1) {
        stop("To prune when a likelihood layer is included, the 'object' must be generated by dgp() with the argument 'depth' of at least 3.", call. = FALSE)
      }
    }
    for (i in 2:n_layer){
      for (ker in object$constructor_obj$all_layer[[i]]){
        if (is.null(ker$global_input)) {
          stop("To prune, 'object' must be generated by dgp() with 'connect = TRUE'.", call. = FALSE)
        }
      }
    }
  }

  is.finish <- FALSE
  cropping_times <- 0
  while (!is.finish){
    crop_id_list <- create_drop_list(object)
    r2 <- object$constructor_obj$aggregate_r2()
    N_cropped <- 0
    for (l in length(crop_id_list):1){
      crop_id <- (r2[[l+1]][[1]]>con$threshold)
      crop_id_list[[l]] <- crop_id
      N_cropped <- N_cropped + sum(crop_id)
    }
    if (N_cropped!=0) {
      object <- copy_in_design(object)
      object <- crop(object, crop_id_list, refit_cores = as.integer(1), verb = verb)
      if ( !inherits(object,"dgp") ) {
        is.finish <- TRUE
      } else {
        n_layer <- object$constructor_obj$n_layer
        if (object$constructor_obj$all_layer[[n_layer]][[1]]$type!='gp') {
          n_layer <- n_layer - 1
          if (n_layer == 1) is.finish <- TRUE
        }
      }
      cropping_times <- cropping_times + 1
    } else {
      is.finish <- TRUE
    }
  }
  if (cropping_times == 0) {
    if (verb) message("No GP nodes can be pruned.", appendLF = FALSE)
  } else {
    if ('loo' %in% names(object)) object[['loo']] <- NULL
    if ('oos' %in% names(object)) object[['oos']] <- NULL
    if ('results' %in% names(object)) object[['results']] <- NULL
    if (verb) message(" * No more GP nodes can be pruned.", appendLF = FALSE)
  }
  return(object)
}

extract_specs <- function(obj, type) {
  res <- list()
  if ( type=='gp' ){
    res[['kernel']] = obj$kernel$name
    res[['lengthscales']] = as.vector(obj$kernel$length)
    res[['scale']] = as.numeric(obj$kernel$scale)
    res[['nugget']] = as.numeric(obj$kernel$nugget)
  } else if ( type=='dgp' ){
    no_layer <- length(obj)
    for (l in 1:no_layer){
      no_node <- length(obj[[l]])
      for (k in 1:no_node)
        if ( obj[[l]][[k]]$type == 'likelihood' ){
          res[[paste('layer', l, sep="")]][[paste('node', k, sep="")]][['type']] <- obj[[l]][[k]]$name
        } else {
          res[[paste('layer', l, sep="")]][[paste('node', k, sep="")]][['kernel']] <- obj[[l]][[k]]$name
          res[[paste('layer', l, sep="")]][[paste('node', k, sep="")]][['lengthscales']] <- as.vector(obj[[l]][[k]]$length)
          res[[paste('layer', l, sep="")]][[paste('node', k, sep="")]][['scale']] <- as.numeric(obj[[l]][[k]]$scale)
          res[[paste('layer', l, sep="")]][[paste('node', k, sep="")]][['nugget']] <- as.numeric(obj[[l]][[k]]$nugget)
        }
    }
  }
  return(res)
}

crop <- function(object, crop_id_list, refit_cores, verb) {
  total_layer <- object$constructor_obj$n_layer
  all_layer <- object$constructor_obj$all_layer
  blocked_gibbs <- object$constructor_obj$block
  B <- as.integer(length(object$emulator_obj$all_layer_set))
  X <- object[['data']][['X']]
  Y <- object[['data']][['Y']]
  linked_idx <- object$container_obj$local_input_idx
  n_layer <- length(crop_id_list)
  for (i in n_layer:1){
    if (all(crop_id_list[[i]])){
      if ( (total_layer-i)==1 ){
        nodes <- all_layer[[total_layer]]
        if (length(nodes)==1){
          if ( verb ) message(paste(" - Transiting to a GP emulator ...", collapse=" "), appendLF = FALSE)
          struc <- nodes[[1]]
          struc$input_dim <- reticulate::np_array(as.integer(object[['specs']][['internal_dims']] - 1))
          struc$connect <- if( isFALSE(object[['specs']][['external_dims']]) ) NULL else reticulate::np_array(as.integer(object[['specs']][['external_dims']]-1))
          struc$input <- X[,struc$input_dim+1,drop=F]
          if ( is.null(struc$connect) ){
            struc$global_input <- NULL
          } else {
            struc$global_input <- X[,struc$connect+1,drop=F]
          }
          if (length(struc$length)!=1) {
            length_dim <- ncol(X)
            struc$length <- utils::tail(struc$length, length_dim)
            struc$para_path <- matrix(c(struc$scale, struc$length, struc$nugget), nrow = 1, byrow=T)
          }
          obj <- pkg.env$dgpsi$gp(X, Y, struc)
          with(pkg.env$np$errstate(divide = 'ignore'), obj$train())
          res <- list()
          res[['id']] <- object$id
          res[['data']][['X']] <- X
          res[['data']][['Y']] <- Y
          res[['specs']] <- extract_specs(obj, "gp")
          res[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
          res[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
          res[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
          res[['constructor_obj']] <- obj
          res[['container_obj']] <- pkg.env$dgpsi$container(obj$export(), linked_idx)
          res[['emulator_obj']] <- obj
          class(res) <- "gp"
          if ( verb ) message(" done")
          return(res)
        } else {
          if ( verb ) message(paste(" - Transiting to a bundle of GP emulators ...", collapse=" "), appendLF = FALSE)
          gp_list <- vector('list', length(nodes))
          for (j in 1:length(nodes)){
            struc <- nodes[[j]]
            struc$input_dim <- reticulate::np_array(as.integer(object[['specs']][['internal_dims']] - 1))
            struc$connect <- if( isFALSE(object[['specs']][['external_dims']]) ) NULL else reticulate::np_array(as.integer(object[['specs']][['external_dims']]-1))
            struc$input <- X[,struc$input_dim+1,drop=F]
            if ( is.null(struc$connect) ){
              struc$global_input <- NULL
            } else {
              struc$global_input <- X[,struc$connect+1,drop=F]
            }
            if (length(struc$length)!=1) {
              length_dim <- ncol(X)
              struc$length <- utils::tail(struc$length, length_dim)
              struc$para_path <- matrix(c(struc$scale, struc$length, struc$nugget), nrow = 1, byrow = T)
            }
            obj <- pkg.env$dgpsi$gp(X, Y[,j,drop=F], struc)
            with(pkg.env$np$errstate(divide = 'ignore'), obj$train())
            res_j <- list()
            res_j[['id']] <- uuid::UUIDgenerate()
            res_j[['data']][['X']] <- X
            res_j[['data']][['Y']] <- Y[,j,drop=F]
            res_j[['specs']] <- extract_specs(obj, "gp")
            res_j[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
            res_j[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
            res_j[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
            res_j[['constructor_obj']] <- obj
            res_j[['container_obj']] <- pkg.env$dgpsi$container(obj$export(), linked_idx)
            res_j[['emulator_obj']] <- obj
            class(res_j) <- "gp"
            gp_list[[j]] <- res_j
          }
          res <- pack(gp_list, id = object$id)
          if ( verb ) message(" done")
          return(res)
        }
      } else {
        if ( verb ) message(paste(c(" - Pruning all nodes below layer", sprintf("%i", i+1), "..."), collapse=" "), appendLF = FALSE)
        all_layer <- all_layer[-c(1:i)]
        for ( j in 1:length(all_layer[[1]]) ){
          all_layer[[1]][[j]]$input_dim <- reticulate::np_array(as.integer(object[['specs']][['internal_dims']] - 1))
          all_layer[[1]][[j]]$connect <- if( isFALSE(object[['specs']][['external_dims']]) ) NULL else reticulate::np_array(as.integer(object[['specs']][['external_dims']]-1))
          all_layer[[1]][[j]]$input <- object$constructor_obj$X[,all_layer[[1]][[j]]$input_dim+1,drop=F]
          all_layer[[1]][[j]]$R2 <- NULL
          if ( is.null(all_layer[[1]][[j]]$connect) ){
            all_layer[[1]][[j]]$global_input <- NULL
          } else {
            all_layer[[1]][[j]]$global_input <- object$constructor_obj$X[,all_layer[[1]][[j]]$connect+1,drop=F]
          }
          if (length(all_layer[[1]][[j]]$length)!=1) {
            length_dim <- ncol(all_layer[[1]][[j]]$input)+ifelse( is.null(all_layer[[1]][[j]]$global_input), 0, ncol(all_layer[[1]][[j]]$global_input) )
            all_layer[[1]][[j]]$length <- utils::tail(all_layer[[1]][[j]]$length, length_dim)
          }
        }
        if ( verb ) message(" done")
        object$constructor_obj$update_all_layer(all_layer)
        if ( verb ) message(" - Re-fitting ...", appendLF = FALSE)
        if ( identical(refit_cores, as.integer(1)) ){
          with(pkg.env$np$errstate(divide = 'ignore'), object$constructor_obj$train(as.integer(100), as.integer(10), TRUE))
        } else {
          with(pkg.env$np$errstate(divide = 'ignore'), object$constructor_obj$ptrain(as.integer(100), as.integer(10), TRUE, refit_cores))
        }
        est_obj <- object$constructor_obj$estimate()
        internal_dims <- object[['specs']][['internal_dims']]
        external_dims <- object[['specs']][['external_dims']]
        object[['specs']] <- extract_specs(est_obj, "dgp")
        object[['specs']][['internal_dims']] <- internal_dims
        object[['specs']][['external_dims']] <- external_dims
        object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
        id <- sample.int(100000, 1)
        set_seed(id)
        object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = blocked_gibbs)
        object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, blocked_gibbs)
        object[['specs']][['seed']] <- id
        if ( verb ) message(" done")
        return(object)
      }
    } else {
      if (any(crop_id_list[[i]])){
        if ( verb ) message(paste(c(" - Pruning", sprintf("%i", sum(crop_id_list[[i]])), "node(s) in layer", sprintf("%i", i), "..."), collapse=" "), appendLF = FALSE)
        all_layer[[i]] <- all_layer[[i]][!crop_id_list[[i]]]
        for ( j in 1:length(all_layer[[i+1]]) ){
          all_layer[[i+1]][[j]]$input_dim <- reticulate::np_array(as.integer((1:sum(!crop_id_list[[i]])) - 1))
          all_layer[[i+1]][[j]]$input <- all_layer[[i+1]][[j]]$input[,!crop_id_list[[i]],drop=F]
          if (length(all_layer[[i+1]][[j]]$length)!=1) {
            all_layer[[i+1]][[j]]$length <- all_layer[[i+1]][[j]]$length[c(!crop_id_list[[i]], rep(T, ncol(all_layer[[i+1]][[j]]$global_input)))]
          }
        }
        if ( verb ) message(" done")
      }
    }
  }

  object$constructor_obj$update_all_layer(all_layer)
  if ( verb ) message(" - Re-fitting ...", appendLF = FALSE)
  if ( identical(refit_cores, as.integer(1)) ){
    with(pkg.env$np$errstate(divide = 'ignore'), object$constructor_obj$train(as.integer(100), as.integer(10), TRUE))
  } else {
    with(pkg.env$np$errstate(divide = 'ignore'), object$constructor_obj$ptrain(as.integer(100), as.integer(10), TRUE, refit_cores))
  }
  est_obj <- object$constructor_obj$estimate()
  internal_dims <- object[['specs']][['internal_dims']]
  external_dims <- object[['specs']][['external_dims']]
  object[['specs']] <- extract_specs(est_obj, "dgp")
  object[['specs']][['internal_dims']] <- internal_dims
  object[['specs']][['external_dims']] <- external_dims
  object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
  id <- sample.int(100000, 1)
  set_seed(id)
  object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = blocked_gibbs)
  object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, blocked_gibbs)
  object[['specs']][['seed']] <- id
  if ( verb ) message(" done")
  return(object)
}
