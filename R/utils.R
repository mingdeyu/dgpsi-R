#' @title Combine layers
#'
#' @description
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated and will be removed in the next release, as it is
#' simply a wrapper for the [list()] function. To construct linked (D)GP structures,
#' please use the updated [lgp()] function, which provides a simpler and more efficient
#' approach to building (D)GP emulators.
#'
#' @param ... a sequence of lists. Each list represents a system layer and contains emulators (produced by [gp()] or
#'    [dgp()]) in that layer.
#'
#' @return A list defining a linked (D)GP structure to be passed to `struc` of [lgp()].
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @md
#' @keywords internal
#' @export
combine <- function(...) {
  lifecycle::deprecate_warn(
    when = "2.5.0",
    what = "combine()",
    details = c(i = "The function will be removed in the next release.",
                i = "To construct linked (D)GP structure, please use the updated `lgp()` function instead."
    )
  )
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
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # load packages
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
#' m1 <- dgp(X, Y[,1], connect = F)
#' m2 <- dgp(X, Y[,2], connect = F)
#'
#' # specify the range of the input dimension
#' lim <- c(0, 1)
#'
#' # pack emulators to form an emulator bundle
#' m <- pack(m1, m2)
#'
#' # 1st wave of the sequential design with 10 iterations and the target RMSE of 0.01
#' m <- design(m, N = 10, limits = lim, f = f, x_test = validate_x, y_test = validate_y, target = 0.01)
#'
#' # 2rd wave of the sequential design with additional 10 iterations and the same target
#' m <- design(m, N = 10, limits = lim, f = f, x_test = validate_x, y_test = validate_y, target = 0.01)
#'
#' # draw sequential designs of the two packed emulators
#' draw(m, type = 'design')
#'
#' # inspect the traces of RMSEs of the two packed emulators
#' draw(m, type = 'rmse')
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
#' @description This function unpacks a bundle of (D)GP emulators safely so that any further manipulations of unpacked individual emulators
#'     will not impact those in the bundle.
#'
#' @param object an instance of the class `bundle`.
#'
#' @return A named list that contains individual emulators (named `emulator1,...,emulatorS`) packed in `object`,
#'    where `S` is the number of emulators in `object`.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
#' @param light a bool indicating if a light version of the constructed emulator
#'     (that requires less disk space to store) will be saved. Defaults to `TRUE`.
#'
#' @return No return value. `object` will be saved to a local `.pkl` file specified by `pkl_file`.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @note Since emulators built from the package are 'python' objects, [save()] from R will not work as it would for R objects. If `object`
#'     was processed by [set_vecchia()] to add or remove the Vecchia approximation, `light` should be set to `FALSE` to ensure
#'     reproducibility after the saved emulator is reloaded by [read()].
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
    if (inherits(object,"gp")){
      object[['container_obj']] <- NULL
    } else if (inherits(object,"dgp")){
      if ( !"seed" %in% names(object$specs) ) stop("The supplied 'object' cannot be saved in light mode. To save, either set 'light = FALSE' or produce a new version of 'object' by set_imp().", call. = FALSE)
      object[['emulator_obj']] <- NULL
      object[['container_obj']] <- NULL
    } else if (inherits(object,"lgp")){
      if ( !"seed" %in% names(object$specs) ) stop("The supplied 'object' cannot be saved in light mode. To save, either set 'light = FALSE' or re-construct and activate the 'object' by lgp().", call. = FALSE)
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
        } else {
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
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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

#' @title Set Emulator ID
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function assigns a unique identifier to an emulator.
#'
#' @param object an emulator object to which the ID will be assigned.
#' @param id a unique identifier for the emulator as either a numeric or character
#'     string. Ensure this ID does not conflict with other emulator IDs, especially
#'     when used in linked emulations.
#'
#' @return The updated `object`, with the assigned ID stored in its `id` slot.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See lgp() for an example.
#' }
#' @md
#' @export
set_id <- function(object, id) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  object[['id']] <- id
  return(object)
}

#' @title Set the number of threads
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function sets the number of threads for parallel computations involved
#'    in the package.
#'
#' @param num the number of threads. If it is greater than the maximum number of threads available, the
#'     number of threads will be set to the maximum value.
#'
#' @return No return value.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @md
#' @export
set_thread_num <- function(num) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if (num >= pkg.env$thread_num) num <- pkg.env$thread_num
  num <- as.integer(num)
  pkg.env$dgpsi$set_thread(num)
}

#' @title Get the number of threads
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function gets the number of threads used for parallel computations involved
#'    in the package.
#'
#' @return the number of threads.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @md
#' @export
get_thread_num <- function() {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  thread_num <- pkg.env$dgpsi$get_thread()
  return(thread_num)
}

#' @title Load the stored emulator
#'
#' @description This function loads the `.pkl` file that stores the emulator.
#'
#' @param pkl_file the path to and the name of the `.pkl` file where the emulator is stored.
#'
#' @return The S3 class of a GP emulator, a DGP emulator, a linked (D)GP emulator, or a bundle of (D)GP emulators.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), lgp(), or pack() for an example.
#' }
#' @md
#' @export
read <- function(pkl_file) {
  # in next release remove the loading of linked_idx as this will be no longer needed.
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
      if ('container_obj' %in% names(res)){
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        class(res) <- "gp"
      } else {
        est_obj <- res$constructor_obj$export()
        if (is.null(res[['specs']][['linked_idx']])){
          linked_idx <- NULL
        } else {
          linked_idx <- if ( isFALSE( res[['specs']][['linked_idx']]) ) {NULL} else {res[['specs']][['linked_idx']]}
        }
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        res[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_r_to_py(linked_idx))
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "gp"
      }
    } else if (label == "dgp"){
      if ('emulator_obj' %in% names(res)){
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        class(res) <- "dgp"
      } else {
        burnin <- res$constructor_obj$burnin
        est_obj <- res$constructor_obj$estimate(burnin)
        B <- res$specs$B
        isblock <- res$constructor_obj$block
        if (is.null(res[['specs']][['linked_idx']])){
          linked_idx <- NULL
        } else {
          linked_idx <- if ( isFALSE( res[['specs']][['linked_idx']]) ) {NULL} else {res[['specs']][['linked_idx']]}
        }
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
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
        if ('seed' %in% names(res$specs)){
          B <- res$specs$B
          extracted_struc <- res$constructor_obj
          set_seed(res$specs$seed)
          obj <- pkg.env$dgpsi$lgp(all_layer = extracted_struc, N = B)
          res[['emulator_obj']] <- obj
          if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        } else {
          if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        }
        if ('metadata' %in% names(res$specs)){
          res$specs$metadata <- as.data.frame(res$specs$metadata)
          res$specs$struc <- as.data.frame(res$specs$struc)
        }
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
            if ('container_obj' %in% names(res[[paste('emulator',i, sep='')]])){
              class(res[[paste('emulator',i, sep='')]]) <- "gp"
            } else {
              est_obj <- res[[paste('emulator',i, sep='')]]$constructor_obj$export()
              if (is.null(res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']])){
                linked_idx <- NULL
              } else {
                linked_idx <- if ( isFALSE( res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]) ) {NULL} else {res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]}
              }
              res[[paste('emulator',i, sep='')]][['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_r_to_py(linked_idx))
              class(res[[paste('emulator',i, sep='')]]) <- "gp"
            }
          }
        } else {
          burnin <- res[[paste('emulator',i, sep='')]]$constructor_obj$burnin
          est_obj <- res[[paste('emulator',i, sep='')]]$constructor_obj$estimate(burnin)
          B <- res[[paste('emulator',i, sep='')]]$specs$B
          isblock <- res[[paste('emulator',i, sep='')]]$constructor_obj$block
          if (is.null(res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']])){
            linked_idx <- NULL
          } else {
            linked_idx <- if ( isFALSE( res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]) ) {NULL} else {res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]}
          }
          set_seed(res[[paste('emulator',i, sep='')]]$specs$seed)
          res[[paste('emulator',i, sep='')]][['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
          res[[paste('emulator',i, sep='')]][['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_r_to_py(linked_idx), isblock)
          class(res[[paste('emulator',i, sep='')]]) <- "dgp"
        }
        if (!'id' %in% names(res[[paste('emulator',i, sep='')]])) res[[paste('emulator',i, sep='')]][['id']] <- uuid::UUIDgenerate()
        if (!'vecchia' %in% names(res[[paste('emulator',i, sep='')]]$specs)) {
          res[[paste('emulator',i, sep='')]][['specs']][['vecchia']] <- FALSE
          res[[paste('emulator',i, sep='')]][['specs']][['M']] <- 25
        }
      }
    }
  } else {
    if ('emulator_obj' %in% names(res)){
      type <- pkg.env$py_buildin$type(res$emulator_obj)$'__name__'
      if ( type=='emulator' ) {
        if (!'specs' %in% names(res)) {
          est_obj <- res$emulator_obj$all_layer
          res[['specs']] <- extract_specs(est_obj, "dgp")
        }
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        if (!'B' %in% names(res$specs)) {
          res[['specs']][['B']] <- length(res$emulator_obj$all_layer_set)
        }
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "dgp"
      } else if ( type=='gp' ) {
        if (!'specs' %in% names(res)) {
          res[['specs']] <- extract_specs(res[['constructor_obj']], "gp")
        }
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        class(res) <- "gp"
      } else if ( type=='lgp' ) {
        if (!'specs' %in% names(res)) {
          res[['specs']][['B']] <- length(res$emulator_obj$all_layer_set)
        }
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
          if (!'specs' %in% names(res[[paste('emulator',i, sep='')]])) {
            est_obj <- res[[paste('emulator',i, sep='')]]$emulator_obj$all_layer
            res[[paste('emulator',i, sep='')]][['specs']] <- extract_specs(est_obj, "dgp")
          }
          if (!'vecchia' %in% names(res[[paste('emulator',i, sep='')]]$specs)) {
            res[[paste('emulator',i, sep='')]][['specs']][['vecchia']] <- FALSE
            res[[paste('emulator',i, sep='')]][['specs']][['M']] <- 25
          }
          if (!'B' %in% names(res[[paste('emulator',i, sep='')]]$specs)) {
            res[[paste('emulator',i, sep='')]][['specs']][['B']] <- length(res[[paste('emulator',i, sep='')]]$emulator_obj$all_layer_set)
          }
          if (!'id' %in% names(res[[paste('emulator',i, sep='')]])) res[[paste('emulator',i, sep='')]][['id']] <- uuid::UUIDgenerate()
          class(res[[paste('emulator',i, sep='')]]) <- "dgp"
        } else if ( type=='gp' ) {
          if (!'specs' %in% names(res[[paste('emulator',i, sep='')]])) {
            res[[paste('emulator',i, sep='')]][['specs']] <- extract_specs(res[[paste('emulator',i, sep='')]][['constructor_obj']], "gp")
          }
          if (!'vecchia' %in% names(res[[paste('emulator',i, sep='')]]$specs)) {
            res[[paste('emulator',i, sep='')]][['specs']][['vecchia']] <- FALSE
            res[[paste('emulator',i, sep='')]][['specs']][['M']] <- 25
          }
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
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function provides a summary of key information for a GP, DGP, or linked (D)GP emulator
#' by generating either a table or an interactive plot of the emulator’s structure.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param type a character string, either `"table"` or `"plot"`, indicating the format of the output.
#'     If set to `"table"`, the function returns a summary in table. If set to `"plot"`, the function
#'     returns an interactive visualization. Defaults to `"plot"`. If the `object` was created with
#'     [lgp()] where `struc` is not a data frame, `type` will automatically default to `"table"`.
#' @param group_size an integer specifying the number of consecutive layers to be grouped together
#'     in the interactive visualization of linked emulators when `type = "plot"`.
#'     This argument is only applicable if `object` is an instance of the `lgp` class.
#'     Defaults to `1`.
#' @param ... Any arguments that can be passed to [kableExtra::kbl()] when `type = "table"`.
#'
#' @return Either a summary table (returned as `kableExtra` object) or an interactive visualization
#' (returned as a `visNetwork` object) of the emulator. The visualization is compatible with R Markdown
#' documents and the RStudio Viewer. The summary table can be further customized by [kableExtra::kableExtra] package.
#' The resulting `visNetwork` object can be saved as an HTML file using [visNetwork::visSave()] from the [visNetwork::visNetwork] package.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
summary.gp <- function(object, type = "plot", ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  if ( type!='plot' & type!='table' ) stop("'type' can only be either 'plot' or 'table'.", call. = FALSE)

  if (type == "table"){
    nodes <- data.frame(
      KernelType = ifelse(object$specs$kernel=='sexp', "Squared Exp", "Mat\u00e9rn-2.5"),
      InputDims = length(object$emulator_obj$kernel$input_dim),
      OutputDims = 1,
      LengthScales = paste(format(object$specs$lengthscales, digits = 3, nsmall = 3), collapse = ", "),
      Variance = format(object$specs$scale, digits = 3, nsmall = 3),
      Nugget = format(object$specs$nugget, digits = 3, nsmall = 3, scientific = TRUE),
      Vecchia = ifelse(object$specs$vecchia, "ON", "OFF"),
      stringsAsFactors = FALSE
    )
    defaults <- list(
      col.names = c(
        "Kernel", "Input Dim(s)", "Output Dim",
        "Length-scale(s)", "Scale (Prior Var)", "Nugget", "Vecchia"
      ),
      format = "simple",
      row.names = FALSE,
      align = rep("c", ncol(nodes)),
      caption = "Summary of GP Emulator",
      escape = FALSE
    )

    # Capture user-supplied arguments as a list
    user_args <- list(...)

    # Merge defaults with user-supplied arguments (only if not already provided)
    final_args <- defaults
    for (name in names(user_args)) {
      final_args[[name]] <- user_args[[name]]
    }

    # Call knitr::kable with the merged arguments
    do.call(kableExtra::kbl, c(list(nodes), final_args))
  } else {

    c24_rgba_lighter <- c(
      "rgba(115, 184, 255, 1)",   # lighter dodgerblue2
      "rgba(240, 115, 115, 1)",   # lighter #E31A1C
      "rgba(102, 180, 102, 1)",   # lighter green4
      "rgba(161, 115, 186, 1)",   # lighter #6A3D9A
      "rgba(255, 180, 100, 1)",   # lighter #FF7F00
      "rgba(128, 128, 128, 1)",   # lighter black
      "rgba(255, 233, 100, 1)",   # lighter gold1
      "rgba(180, 225, 245, 1)",   # lighter skyblue2
      "rgba(255, 190, 190, 1)",   # lighter #FB9A99
      #"rgba(180, 255, 180, 1)",   # lighter palegreen2
      "rgba(225, 205, 235, 1)",   # lighter #CAB2D6
      "rgba(255, 220, 170, 1)",   # lighter #FDBF6F
      #"rgba(245, 240, 190, 1)",   # lighter khaki2
      "rgba(160, 60, 60, 1)",     # lighter maroon
      "rgba(240, 150, 240, 1)",   # lighter orchid1
      "rgba(255, 130, 185, 1)",   # lighter deeppink1
      "rgba(100, 100, 255, 1)",   # lighter blue1
      "rgba(120, 150, 180, 1)",   # lighter steelblue4
      "rgba(100, 230, 230, 1)",   # lighter darkturquoise
      #"rgba(170, 255, 170, 1)",   # lighter green1
      "rgba(200, 200, 100, 1)",   # lighter yellow4
      #"rgba(255, 255, 150, 1)",   # lighter yellow3
      "rgba(210, 130, 100, 1)",   # lighter darkorange4
      "rgba(200, 110, 110, 1)"    # lighter brown
    )

    c24_rgba <- c(
      "rgba(30, 144, 255, 1)",    # dodgerblue2
      "rgba(227, 26, 28, 1)",     # #E31A1C
      "rgba(0, 139, 0, 1)",       # green4
      "rgba(106, 61, 154, 1)",    # #6A3D9A
      "rgba(255, 127, 0, 1)",     # #FF7F00
      "rgba(0, 0, 0, 1)",         # black
      "rgba(255, 215, 0, 1)",     # gold1
      "rgba(135, 206, 235, 1)",   # skyblue2
      "rgba(251, 154, 153, 1)",   # #FB9A99
      #"rgba(144, 238, 144, 1)",   # palegreen2
      "rgba(202, 178, 214, 1)",   # #CAB2D6
      "rgba(253, 191, 111, 1)",   # #FDBF6F
      #"rgba(240, 230, 140, 1)",   # khaki2
      "rgba(128, 0, 0, 1)",       # maroon
      "rgba(218, 112, 214, 1)",   # orchid1
      "rgba(255, 20, 147, 1)",    # deeppink1
      "rgba(0, 0, 255, 1)",       # blue1
      "rgba(54, 100, 139, 1)",    # steelblue4
      "rgba(0, 206, 209, 1)",     # darkturquoise
      #"rgba(0, 255, 0, 1)",       # green1
      "rgba(139, 139, 0, 1)",     # yellow4
      #"rgba(255, 255, 0, 1)",     # yellow3
      "rgba(139, 69, 19, 1)",     # darkorange4
      "rgba(165, 42, 42, 1)"      # brown
    )

    nodes <- data.frame(
      id = c("Global", "gp"),
      label = c("", "GP"),
      title = c(paste0(
        "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333; text-align: center;'>",
        "<b style='font-size: 12px; color:", c24_rgba[1] ,";'>Global Input</b> ", "<br>",
        "<span style='font-size: 12px; color:", c24_rgba[1] ,";'>Total Dim(s):</span> ", ncol(object$data$X),
        "</div>"
      ),
      paste0(
        "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333; white-space: normal; max-width: 200px;'>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Total Input Dim(s):</b> ", length(object$emulator_obj$kernel$input_dim), "<br>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Total Output Dim:</b> ", 1, "<br>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Kernel Type:</b> ", ifelse(object$specs$kernel=='sexp', "Squared Exp", "Mat&eacute;rn-2.5"), "<br>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Length-scales:</b> ", paste(format(object$specs$lengthscales, digits = 3, nsmall = 3), collapse = ", "), "<br>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Scale (Prior Variance):</b> ",format(object$specs$scale, digits = 3, nsmall = 3), "<br>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Nugget:</b> ", format(object$specs$nugget, digits = 3, nsmall = 3, scientific = TRUE), "<br>",
        "<b style='font-size: 12px; color: ", c24_rgba[-1][1], ";'>Vecchia Mode:</b> ", ifelse(object$specs$vecchia, "ON", "OFF"),
        "</div>"
      )
      ),
      level = c(0, 1),
      group = c(as.character(0), as.character(1)),
      stringsAsFactors = FALSE,
      shape = c('database', "circle")
    )

    edge <- data.frame(
      from = "Global",
      to = "gp",
      label = "",
      title = paste0(
        "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
        "<b style='font-size: 12px;'>From Dim(s):</b> ", paste(object$emulator_obj$kernel$input_dim + 1, collapse=", "), "<br>",
        "<b style='font-size: 12px;'>To Dim(s):</b> ", paste(1:length(object$emulator_obj$kernel$input_dim), collapse=", "),
        "</div>"
      ),
      arrows = "to",
      stringsAsFactors = FALSE
    )


    network <- visNetwork::visNetwork(nodes, edge, width = "100%", height = "300px",
                                      main = "Gaussian Process Emulator",
                                      submain = "Graphical Summary") %>%
      visNetwork::visEdges(
        arrows = list(to = list(enable = T, scaleFactor = 0.5)),
        color = list(color = "darkgrey", highlight = "black", hover = "black"),
        selectionWidth = 0.8,
        hoverWidt = 0.5
      )

    # Apply color for the first layer explicitly if needed
    network <- network %>%
      visNetwork::visGroups(
        groupname = as.character(0),
        borderWidthSelected = 1.5,
        color = list(
          background = c24_rgba_lighter[1],
          border = c24_rgba[1],
          hover = list(background = c24_rgba[1], border = c24_rgba_lighter[1]),
          highlight = list(background = c24_rgba[1], border = c24_rgba_lighter[1])
        )
      ) %>%
      visNetwork::visGroups(
        groupname = as.character(1),
        borderWidthSelected = 3,
        color = list(
          background = c24_rgba_lighter[-1][1],
          border = c24_rgba_lighter[-1][1],
          hover = list(background = c24_rgba[-1][1], border = c24_rgba[-1][1]),
          highlight = list(background = c24_rgba[-1][1], border = c24_rgba[-1][1])
        )
      )


    network %>%
      visNetwork::visNodes(size=15,
                           font = list(
                             size = 14,              # Font size
                             color = "white",     # Font color
                             face = "helvetica",
                             ital = TRUE
                           )) %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = FALSE),
        nodesIdSelection = FALSE) %>%
      visNetwork::visInteraction(hover = TRUE, hideEdgesOnDrag = FALSE)%>%
      visNetwork::visHierarchicalLayout(direction = "LR")
  }
}

#' @rdname summary
#' @method summary dgp
#' @importFrom magrittr %>%
#' @export
summary.dgp <- function(object, type = "plot", ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  if ( type!='plot' & type!='table' ) stop("'type' can only be either 'plot' or 'table'.", call. = FALSE)

  if (type == "table"){
    nodes <- data.frame()

    # Iterate through layers and nodes to populate the dataframe
    n_layer <- object$constructor_obj$n_layer
    for (i in 1:n_layer) {
      layer_id <- object$specs[[paste0("layer", i)]]

      for (j in 1:length(layer_id)) {
        node_id <- layer_id[[paste0("node", j)]]

        # Default values
        formatted_lengthscales <- NA
        total_input_dim <- NA
        kernel_display <- NA

        if (is.null(node_id[["type"]])) {
          # Format length scales
          formatted_lengthscales <- paste(format(node_id$lengthscales, digits = 3, nsmall = 3), collapse = ", ")

          # Calculate total input dimension
          total_input_dim <- length(object$emulator_obj$all_layer[[i]][[j]]$input_dim)
          global_connect_dim <- object$emulator_obj$all_layer[[i]][[j]]$connect
          if (!is.null(global_connect_dim)) {
            total_input_dim <- total_input_dim + length(global_connect_dim)
          }

          # Determine kernel display
          kernel_display <- ifelse(node_id$kernel == "sexp", "Squared Exp", "Mat\u00e9rn-2.5")
        }

        # Combine Node Type and Kernel Type
        node_type_display <- ifelse(
          is.null(node_id[["type"]]),
          paste0("GP (", kernel_display, ")"),
          paste0("Likelihood (", node_id$type, ")")
        )

        # Append to the nodes dataframe
        nodes <- rbind(nodes, data.frame(
          NodeType = node_type_display,
          Layer = i,
          NodeNo = j,
          InputDims = ifelse(is.null(node_id[["type"]]), total_input_dim, length(object$emulator_obj$all_layer[[i]][[j]]$input_dim)),
          OutputDims = 1,
          LengthScales = formatted_lengthscales,
          Variance = ifelse(is.null(node_id[["type"]]), format(node_id$scale, digits = 3, nsmall = 3), NA),
          Nugget = ifelse(is.null(node_id[["type"]]), format(node_id$nugget, digits = 3, nsmall = 3, scientific = TRUE), NA),
          stringsAsFactors = FALSE
        ))
      }
    }

    # Rename columns to add spaces and make them more readable
    defaults <- list(
      col.names = c(
        "Node Type", "Layer", "Node No", "Input Dim(s)", "Output Dim",
        "Length-scale(s)", "Scale (Prior Var)", "Nugget"
      ),
      format = "simple",
      row.names = FALSE,
      align = rep("c", ncol(nodes)),
      caption = "Summary of DGP Emulator Nodes",
      escape = FALSE
    )

    # Capture user-supplied arguments as a list
    user_args <- list(...)

    # Merge defaults with user-supplied arguments (only if not already provided)
    final_args <- defaults
    for (name in names(user_args)) {
      final_args[[name]] <- user_args[[name]]
    }

    # Call knitr::kable with the merged arguments
    do.call(kableExtra::kbl, c(list(nodes), final_args))
  } else {
    c24_rgba_lighter <- c(
      "rgba(115, 184, 255, 1)",   # lighter dodgerblue2
      "rgba(240, 115, 115, 1)",   # lighter #E31A1C
      "rgba(102, 180, 102, 1)",   # lighter green4
      "rgba(161, 115, 186, 1)",   # lighter #6A3D9A
      "rgba(255, 180, 100, 1)",   # lighter #FF7F00
      "rgba(128, 128, 128, 1)",   # lighter black
      "rgba(255, 233, 100, 1)",   # lighter gold1
      "rgba(180, 225, 245, 1)",   # lighter skyblue2
      "rgba(255, 190, 190, 1)",   # lighter #FB9A99
      #"rgba(180, 255, 180, 1)",   # lighter palegreen2
      "rgba(225, 205, 235, 1)",   # lighter #CAB2D6
      "rgba(255, 220, 170, 1)",   # lighter #FDBF6F
      #"rgba(245, 240, 190, 1)",   # lighter khaki2
      "rgba(160, 60, 60, 1)",     # lighter maroon
      "rgba(240, 150, 240, 1)",   # lighter orchid1
      "rgba(255, 130, 185, 1)",   # lighter deeppink1
      "rgba(100, 100, 255, 1)",   # lighter blue1
      "rgba(120, 150, 180, 1)",   # lighter steelblue4
      "rgba(100, 230, 230, 1)",   # lighter darkturquoise
      #"rgba(170, 255, 170, 1)",   # lighter green1
      "rgba(200, 200, 100, 1)",   # lighter yellow4
      #"rgba(255, 255, 150, 1)",   # lighter yellow3
      "rgba(210, 130, 100, 1)",   # lighter darkorange4
      "rgba(200, 110, 110, 1)"    # lighter brown
    )

    c24_rgba <- c(
      "rgba(30, 144, 255, 1)",    # dodgerblue2
      "rgba(227, 26, 28, 1)",     # #E31A1C
      "rgba(0, 139, 0, 1)",       # green4
      "rgba(106, 61, 154, 1)",    # #6A3D9A
      "rgba(255, 127, 0, 1)",     # #FF7F00
      "rgba(0, 0, 0, 1)",         # black
      "rgba(255, 215, 0, 1)",     # gold1
      "rgba(135, 206, 235, 1)",   # skyblue2
      "rgba(251, 154, 153, 1)",   # #FB9A99
      #"rgba(144, 238, 144, 1)",   # palegreen2
      "rgba(202, 178, 214, 1)",   # #CAB2D6
      "rgba(253, 191, 111, 1)",   # #FDBF6F
      #"rgba(240, 230, 140, 1)",   # khaki2
      "rgba(128, 0, 0, 1)",       # maroon
      "rgba(218, 112, 214, 1)",   # orchid1
      "rgba(255, 20, 147, 1)",    # deeppink1
      "rgba(0, 0, 255, 1)",       # blue1
      "rgba(54, 100, 139, 1)",    # steelblue4
      "rgba(0, 206, 209, 1)",     # darkturquoise
      #"rgba(0, 255, 0, 1)",       # green1
      "rgba(139, 139, 0, 1)",     # yellow4
      #"rgba(255, 255, 0, 1)",     # yellow3
      "rgba(139, 69, 19, 1)",     # darkorange4
      "rgba(165, 42, 42, 1)"      # brown
    )

    nodes <- data.frame(
      id = "Global",
      label = "",
      title = paste0(
        "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333; text-align: center;'>",
        "<b style='font-size: 12px; color:", c24_rgba[1] ,";'>Global Input</b> ", "<br>",
        "<span style='font-size: 12px; color:", c24_rgba[1] ,";'>Total Dim(s):</span> ", ncol(object$data$X),
        "</div>"
      ),
      level = 0,
      group = as.character(0),
      stringsAsFactors = FALSE,
      shape = 'database'
    )

    edge <- data.frame(
      from = character(),
      to = character(),
      label = character(),
      title = character(),
      arrows = character(),
      stringsAsFactors = FALSE
    )

    n_layer <- object$constructor_obj$n_layer
    for (i in 1:n_layer){
      layer_id <- object$specs[[paste0("layer",i)]]
      for (j in 1:length(layer_id)){
        node_id <- layer_id[[paste0("node",j)]]

        if (is.null(node_id[["type"]])){
          formatted_lengthscales <- paste(format(node_id$lengthscales, digits = 3, nsmall = 3), collapse = ", ")
          total_input_dim <- length(object$emulator_obj$all_layer[[i]][[j]]$input_dim)
          global_connect_dim <- object$emulator_obj$all_layer[[i]][[j]]$connect
          if (!is.null(global_connect_dim)) total_input_dim <- total_input_dim + length(global_connect_dim)
        }

        nodes <- rbind(nodes, data.frame(
          id = paste0("e",i,j),
          label = ifelse(is.null(node_id[["type"]]), "GP", "LH"), # Leave the label empty to avoid duplication with the hover title
          title = if (is.null(node_id[["type"]])){
            paste0(
              "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333; white-space: normal; max-width: 200px;'>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Node Type:</b> ", "GP", "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Total Input Dim(s):</b> ", total_input_dim, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Total Output Dim:</b> ", 1, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Kernel Type:</b> ", ifelse(node_id$kernel=='sexp', "Squared Exp", "Mat&eacute;rn-2.5"), "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Length-scales:</b> ", formatted_lengthscales, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Scale (Prior Variance):</b> ",format(node_id$scale, digits = 3, nsmall = 3), "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Nugget:</b> ", format(node_id$nugget, digits = 3, nsmall = 3, scientific = TRUE), "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Layer in Hierarchy:</b> ", i, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Position in Layer:</b> ", j, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Vecchia Mode:</b> ", ifelse(object$specs$vecchia, "ON", "OFF"),
              "</div>"
            )
          } else {
            paste0(
              "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Node Type:</b> ", node_id$type, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Total Input Dim(s):</b> ", length(object$emulator_obj$all_layer[[i]][[j]]$input_dim), "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Total Output Dim:</b> ", 1, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Layer in Hierarchy:</b> ", i, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Position in Layer:</b> ", j, "<br>",
              "<b style='font-size: 12px; color: ", c24_rgba[-1][(i - 1) %% length(c24_rgba[-1]) + 1], ";'>Vecchia Mode:</b> ", ifelse(object$specs$vecchia, "ON", "OFF"),
              "</div>"
            )
          },
          level = i,
          group = as.character(i),
          stringsAsFactors = FALSE,
          shape = ifelse(is.null(node_id[["type"]]), 'circle', 'box')
        )
        )

        if (i==1){
          edge <- rbind(edge, data.frame(
            from = "Global",
            to = paste0("e",i,j),
            label = "",
            title = paste0(
              "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
              "<b style='font-size: 12px;'>From Dim(s):</b> ", paste(object$emulator_obj$all_layer[[i]][[j]]$input_dim + 1, collapse=", "), "<br>",
              "<b style='font-size: 12px;'>To Dim(s):</b> ", paste(1:length(object$emulator_obj$all_layer[[i]][[j]]$input_dim), collapse=", "),
              "</div>"
            ),
            arrows = "to",
            stringsAsFactors = FALSE
          )
          )
        } else {
          connection <- object$emulator_obj$all_layer[[i]][[j]]$input_dim + 1
          if (is.null(node_id[["type"]])){
            global_connection <- object$emulator_obj$all_layer[[i]][[j]]$connect
          } else {
            global_connection <- NULL
          }
          for (k in 1:length(connection)){
            edge <- rbind(edge, data.frame(
              from = paste0("e",i-1,connection[k]),
              to = paste0("e",i,j),
              label = "",
              title = paste0(
                "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
                "<b style='font-size: 12px;'>From Ouput Dim:</b> ", 1, "<br>",
                "<b style='font-size: 12px;'>To Input Dim:</b> ", k,
                "</div>"
              ),
              arrows = "to",
              stringsAsFactors = FALSE
            )
            )
          }

          if (!is.null(global_connection)){
            edge <- rbind(edge, data.frame(
              from = "Global",
              to = paste0("e",i,j),
              label = "",
              title = paste0(
                "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
                "<b style='font-size: 12px;'>From Dim(s):</b> ", paste(global_connection + 1, collapse=", "), "<br>",
                "<b style='font-size: 12px;'>To Dim(s):</b> ", paste((1:length(global_connection)) + length(object$emulator_obj$all_layer[[i]][[j]]$input_dim), collapse=", "),
                "</div>"
              ),
              arrows = "to",
              stringsAsFactors = FALSE
            )
            )
          }
        }
      }
    }

    network <- visNetwork::visNetwork(nodes, edge, , width = "100%", height = "300px",
                                      main = "Deep Gaussian Process Emulator",
                                      submain = "Graphical Summary",
                                      footer = paste0(
                                        "<b>GP</b> = Gaussian Process Node &nbsp;&nbsp;&nbsp;&nbsp; <b>LH</b> = Likelihood Node"
                                      )) %>%
      visNetwork::visEdges(
        arrows = list(to = list(enable = T, scaleFactor = 0.5)),
        color = list(color = "darkgrey", highlight = "black", hover = "black"),
        smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.3),
        selectionWidth = 0.8,
        hoverWidt = 0.5
      )

    # Apply color for the first layer explicitly if needed
    network <- network %>%
      visNetwork::visGroups(
        groupname = as.character(0),
        borderWidthSelected = 1.5,
        color = list(
          background = c24_rgba_lighter[1],
          border = c24_rgba[1],
          hover = list(background = c24_rgba[1], border = c24_rgba_lighter[1]),
          highlight = list(background = c24_rgba[1], border = c24_rgba_lighter[1])
        )
      )

    # Loop over remaining layers, cycling through colors if needed
    for (i in 1:max(unique(nodes$group))) {
      color_index <- (i - 1) %% (length(c24_rgba) - 1) + 1  # Cycle through color indices

      network <- network %>%
        visNetwork::visGroups(
          groupname = as.character(i),
          borderWidthSelected = 3,
          color = list(
            background = c24_rgba_lighter[-1][color_index],
            border = c24_rgba_lighter[-1][color_index],
            hover = list(background = c24_rgba[-1][color_index], border = c24_rgba[-1][color_index]),
            highlight = list(background = c24_rgba[-1][color_index], border = c24_rgba[-1][color_index])
          )
        )
    }

    network %>%
      visNetwork::visNodes(size=15,
                           font = list(
                             size = 14,              # Font size
                             color = "white",     # Font color
                             face = "helvetica",
                             ital = TRUE
                           )) %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = FALSE),
        nodesIdSelection = FALSE) %>%
      visNetwork::visInteraction(hover = TRUE, hideEdgesOnDrag = FALSE)%>%
      visNetwork::visHierarchicalLayout(direction = "LR")
  }
}

#' @rdname summary
#' @method summary lgp
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
summary.lgp <- function(object, type = "plot", group_size = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  if ( type!='plot' & type!='table' ) stop("'type' can only be either 'plot' or 'table'.", call. = FALSE)

  if ( "metadata" %in% names(object$specs) ){
    if (type == 'table') {
      struc <- object$specs$struc
      metadata <- object$specs$metadata
      metadata$Global_Output_Dims <- sapply(metadata$Emulator, function(emulator_id) {
        # Get the total output dimensions for this emulator
        total_out_dim <- metadata$Total_Output_Dims[metadata$Emulator == emulator_id]

        # Full set of output dimensions for this emulator
        all_output_dims <- 1:total_out_dim

        # Find connected output dims in struc
        connected_outputs <- struc$From_Output[struc$From_Emulator == emulator_id]

        # Global output dims are those not in connected outputs
        global_outputs <- setdiff(all_output_dims, connected_outputs)

        # Convert to string if there are multiple or no global outputs
        if (length(global_outputs) == 0) {
          return(NA)  # No global outputs
        } else {
          return(paste(global_outputs, collapse = ", "))
        }
      })

      nodes <- data.frame(
        Emulator = metadata$Emulator,
        Type = toupper(metadata$Type),
        Layer = metadata$Layer,
        Position = metadata$Pos_in_Layer,
        Input_Dims = metadata$Total_Input_Dims,
        Output_Dims = metadata$Total_Output_Dims,
        Global_Output_Indices = ifelse(is.na(metadata$Global_Output_Dims), 'None', metadata$Global_Output_Dims),
        Vecchia = ifelse(metadata$Vecchia, "ON", "OFF"),
        stringsAsFactors = FALSE
      )

      defaults <- list(
        col.names = c(
          "Emulator ID", "Type", "Layer", "Position", "Input Dim(s)",
          "Output Dim(s)", "Global Output Idx", "Vecchia"
        ),
        format = "simple",
        row.names = FALSE,
        align = rep("c", ncol(nodes)),
        caption = "Summary of Linked Emulators"
      )

      # Capture user-supplied arguments as a list
      user_args <- list(...)

      # Merge defaults with user-supplied arguments (only if not already provided)
      final_args <- defaults
      for (name in names(user_args)) {
        final_args[[name]] <- user_args[[name]]
      }

      # Call knitr::kable with the merged arguments
      do.call(kableExtra::kbl, c(list(nodes), final_args))
    } else {
      N <- as.integer(group_size)
      c24_rgba_lighter <- c(
        "rgba(115, 184, 255, 1)",   # lighter dodgerblue2
        "rgba(240, 115, 115, 1)",   # lighter #E31A1C
        "rgba(102, 180, 102, 1)",   # lighter green4
        "rgba(161, 115, 186, 1)",   # lighter #6A3D9A
        "rgba(255, 180, 100, 1)",   # lighter #FF7F00
        "rgba(128, 128, 128, 1)",   # lighter black
        "rgba(255, 233, 100, 1)",   # lighter gold1
        "rgba(180, 225, 245, 1)",   # lighter skyblue2
        "rgba(255, 190, 190, 1)",   # lighter #FB9A99
        #"rgba(180, 255, 180, 1)",   # lighter palegreen2
        "rgba(225, 205, 235, 1)",   # lighter #CAB2D6
        "rgba(255, 220, 170, 1)",   # lighter #FDBF6F
        #"rgba(245, 240, 190, 1)",   # lighter khaki2
        "rgba(160, 60, 60, 1)",     # lighter maroon
        "rgba(240, 150, 240, 1)",   # lighter orchid1
        "rgba(255, 130, 185, 1)",   # lighter deeppink1
        "rgba(100, 100, 255, 1)",   # lighter blue1
        "rgba(120, 150, 180, 1)",   # lighter steelblue4
        "rgba(100, 230, 230, 1)",   # lighter darkturquoise
        #"rgba(170, 255, 170, 1)",   # lighter green1
        "rgba(200, 200, 100, 1)",   # lighter yellow4
        #"rgba(255, 255, 150, 1)",   # lighter yellow3
        "rgba(210, 130, 100, 1)",   # lighter darkorange4
        "rgba(200, 110, 110, 1)"    # lighter brown
      )

      c24_rgba <- c(
        "rgba(30, 144, 255, 1)",    # dodgerblue2
        "rgba(227, 26, 28, 1)",     # #E31A1C
        "rgba(0, 139, 0, 1)",       # green4
        "rgba(106, 61, 154, 1)",    # #6A3D9A
        "rgba(255, 127, 0, 1)",     # #FF7F00
        "rgba(0, 0, 0, 1)",         # black
        "rgba(255, 215, 0, 1)",     # gold1
        "rgba(135, 206, 235, 1)",   # skyblue2
        "rgba(251, 154, 153, 1)",   # #FB9A99
        #"rgba(144, 238, 144, 1)",   # palegreen2
        "rgba(202, 178, 214, 1)",   # #CAB2D6
        "rgba(253, 191, 111, 1)",   # #FDBF6F
        #"rgba(240, 230, 140, 1)",   # khaki2
        "rgba(128, 0, 0, 1)",       # maroon
        "rgba(218, 112, 214, 1)",   # orchid1
        "rgba(255, 20, 147, 1)",    # deeppink1
        "rgba(0, 0, 255, 1)",       # blue1
        "rgba(54, 100, 139, 1)",    # steelblue4
        "rgba(0, 206, 209, 1)",     # darkturquoise
        #"rgba(0, 255, 0, 1)",       # green1
        "rgba(139, 139, 0, 1)",     # yellow4
        #"rgba(255, 255, 0, 1)",     # yellow3
        "rgba(139, 69, 19, 1)",     # darkorange4
        "rgba(165, 42, 42, 1)"      # brown
      )

      struc <- object$specs$struc
      metadata <- object$specs$metadata
      metadata$Global_Output_Dims <- sapply(metadata$Emulator, function(emulator_id) {
        # Get the total output dimensions for this emulator
        total_out_dim <- metadata$Total_Output_Dims[metadata$Emulator == emulator_id]

        # Full set of output dimensions for this emulator
        all_output_dims <- 1:total_out_dim

        # Find connected output dims in struc
        connected_outputs <- struc$From_Output[struc$From_Emulator == emulator_id]

        # Global output dims are those not in connected outputs
        global_outputs <- setdiff(all_output_dims, connected_outputs)

        # Convert to string if there are multiple or no global outputs
        if (length(global_outputs) == 0) {
          return(NA)  # No global outputs
        } else {
          return(paste(global_outputs, collapse = ", "))
        }
      })
      metadata <- rbind(
        data.frame(Emulator = "Global", Layer = 0, Pos_in_Layer = NA,
                   Total_Input_Dims = NA, Total_Output_Dims = NA, Global_Output_Dims = NA,
                   Type = NA, Vecchia = FALSE),
        metadata
      )

      max_layer <- max(metadata$Layer)

      nodes <- data.frame(
        id = metadata$Emulator,
        label = toupper(metadata$Type), # Leave the label empty to avoid duplication with the hover title
        title = paste0(
          "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Emulator ID:</b> ", metadata$Emulator, "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Emulator Type:</b> ", toupper(metadata$Type), "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Total Input Dim(s):</b> ", metadata$Total_Input_Dims, "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Total Output Dim(s):</b> ", metadata$Total_Output_Dims, "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Global Output Indices:</b> ", ifelse(is.na(metadata$Global_Output_Dims), "NA", metadata$Global_Output_Dims), "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Layer in Network:</b> ", metadata$Layer, "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Position in Layer:</b> ", metadata$Pos_in_Layer, "<br>",
          "<b style='font-size: 12px; color: ", c24_rgba[-1][((metadata$Layer - 1) %/% N) %% length(c24_rgba[-1]) + 1], ";'>Vecchia Mode:</b> ", ifelse(metadata$Vecchia, "ON", "OFF"),
          "</div>"
        ),
        level = metadata$Layer,
        group = paste0(
          ((metadata$Layer - 1) %/% N) * N + 1,
          "-",
          pmin(((metadata$Layer - 1) %/% N + 1) * N, max_layer)
        ),
        stringsAsFactors = FALSE,
        shape = 'circle'
      )

      nodes[nodes$level == 0,]$shape <- "database"
      nodes[nodes$level == 0,]$title <- paste0(
        "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
        "<b style='font-size: 12px; color:", c24_rgba[1] ,";'>Global Input</b> ",
        "</div>"
      )
      nodes[nodes$level == 0,]$group <- as.character(0)

      #edges <- data.frame(
      #  from = struc$From_Emulator,
      #  to = struc$To_Emulator,
      #  label = paste(struc$From_Output, "to", struc$To_Input),
      #  arrows = "to",
      #  stringsAsFactors = FALSE
      #)

      edges <- struc %>%
        dplyr::group_by(.data$From_Emulator, .data$To_Emulator) %>%
        dplyr::summarize(
          title = paste0(
            "<div style='font-family: Arial, sans-serif; font-size: 12px; line-height: 1.5em; color: #333;'>",
            "<b style='font-size: 12px;'>From Dim(s):</b> ", paste(.data$From_Output, collapse = ", "), "<br>",
            "<b style='font-size: 12px;'>To Dim(s):</b> ", paste(.data$To_Input, collapse = ", "),
            "</div>"
          ),
          arrows = "to",
          .groups = 'drop'  # This ensures that the result is ungrouped
        ) %>%
        dplyr::rename(from = .data$From_Emulator, to = .data$To_Emulator)  # Rename columns after summarizing

      edges <- as.data.frame(edges, stringsAsFactors = FALSE)

      network <- visNetwork::visNetwork(nodes, edges, , width = "100%", height = "300px",
                                        main = "Network of Linked Emulators",
                                        submain = "Graphical Summary",
                                        footer = paste0(
                                          "<b>GP</b> = Gaussian Process Emulator &nbsp;&nbsp;&nbsp; <b>DGP</b> = Deep Gaussian Process Emulator"
                                        )) %>%
        visNetwork::visEdges(
          arrows = list(to = list(enable = T, scaleFactor = 0.5)),
          color = list(color = "darkgrey", highlight = "black", hover = "black"),
          smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.3),
          selectionWidth = 0.8,
          hoverWidt = 0.5
        )

      # Apply color for the first layer explicitly if needed
      network <- network %>%
        visNetwork::visGroups(
          groupname = as.character(0),
          borderWidthSelected = 1.5,
          color = list(
            background = c24_rgba_lighter[1],
            border = c24_rgba[1],
            hover = list(background = c24_rgba[1], border = c24_rgba_lighter[1]),
            highlight = list(background = c24_rgba[1], border = c24_rgba_lighter[1])
          )
        )

      # Loop over remaining layers, cycling through colors if needed
      unique_groups <- unique(nodes$group[nodes$group != "0"])

      # Sort unique groups in the correct numeric order
      unique_groups <- unique_groups[order(as.numeric(sub("-.*", "", unique_groups)))]

      for (i in 1:length(unique_groups)) {
        color_index <- (i - 1) %% (length(c24_rgba) - 1) + 1  # Cycle through color indices

        network <- network %>%
          visNetwork::visGroups(
            groupname = unique_groups[i],
            borderWidthSelected = 3,
            color = list(
              background = c24_rgba_lighter[-1][color_index],
              border = c24_rgba_lighter[-1][color_index],
              hover = list(background = c24_rgba[-1][color_index], border = c24_rgba[-1][color_index]),
              highlight = list(background = c24_rgba[-1][color_index], border = c24_rgba[-1][color_index])
            )
          )
      }

      network <- network %>%
        visNetwork::visNodes(size=15,
                             font = list(
                               size = 14,              # Font size
                               color = "white",     # Font color
                               face = "helvetica",
                               ital = TRUE
                             )) %>%
        visNetwork::visOptions(
          highlightNearest = list(enabled = FALSE),
          nodesIdSelection = list(enabled = TRUE, useLabels = FALSE, main = "Select by ID")) %>%
        visNetwork::visInteraction(hover = TRUE, hideEdgesOnDrag = FALSE) %>%
        visNetwork::visHierarchicalLayout(direction = "LR")
      if (N!=1){
        network <- network %>%
        visNetwork::visClusteringByGroup(groups = unique(nodes$group), label = "Layer ")
      }
      network
    }
  } else {
    pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
    if (type == 'plot') {
      warning("'type' has been automatically set to 'table' because 'object' was created by lgp() with 'struc' not in data frame format.", call. = FALSE)
    }
  }
}

#' @title Add or remove the Vecchia approximation
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function adds or removes the Vecchia approximation from a GP, DGP or linked (D)GP emulator
#'     constructed by [gp()], [dgp()] or [lgp()].
#'
#' @param object an instance of the S3 class `gp`, `dgp`, or `lgp`.
#' @param vecchia a bool or a list of bools to indicate the addition or removal of the Vecchia approximation:
#' * if `object` is an instance of the `gp` or `dgp` class, `vecchia` is a bool that indicates
#'   either addition (`vecchia = TRUE`) or removal (`vecchia = FALSE`) of the Vecchia approximation from `object`.
#' * if `object` is an instance of the `lgp` class, `x` can be a bool or a list of bools:
#'   - if `vecchia` is a bool, it indicates either addition (`vecchia = TRUE`) or removal (`vecchia = FALSE`) of
#'     the Vecchia approximation from all individual (D)GP emulators contained in `object`.
#'   - if `vecchia` is a list of bools, it should have same shape as `struc` that was supplied to [lgp()]. Each bool
#'     in the list indicates if the corresponding (D)GP emulator contained in `object` shall have the Vecchia approximation
#'     added or removed.
#' @param M the size of the conditioning set for the Vecchia approximation in the (D)GP emulator training. Defaults to `25`.
#' @param ord an R function that returns the ordering of the input to the (D)GP emulator for the Vecchia approximation. The
#'    function must satisfy the following basic rules:
#' * the first argument represents the lengthscale-scaled input to the GP emulator or the lengthscale-scaled input to a GP node
#'   of the DGP emulator.
#' * the output of the function is a vector of indices that gives the ordering of the input to the GP emulator or the input to
#'   the GP nodes of the DGP emulator.
#'
#' If `ord = NULL`, the default random ordering is used. Defaults to `NULL`.
#' @return An updated `object` with the Vecchia approximation either added or removed.
#'
#' @note This function is useful for quickly switching between Vecchia and non-Vecchia approximations for an existing emulator
#'     without the need to reconstruct the emulator. If the emulator was built without the Vecchia approximation, the function
#'     can add it, and if the emulator was built with the Vecchia approximation, the function can remove it. If the current
#'     state already matches the requested state, the emulator remains unchanged.
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @md
#' @export
set_vecchia <- function(object, vecchia = TRUE, M = 25, ord = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
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

  if ( inherits(object,"gp") ){
    if (vecchia){
      if (object$specs$vecchia) {
        return(object)
      } else {
        object$constructor_obj$to_vecchia(m = M, ord_fun = ord_wrapper)
        object$container_obj$to_vecchia()
        object[['specs']][['vecchia']] <- TRUE
        object[['specs']][['M']] <- M
        return(object)
      }
    } else {
      if (object$specs$vecchia) {
        object$constructor_obj$remove_vecchia()
        object$container_obj$remove_vecchia()
        object[['specs']][['vecchia']] <- FALSE
        return(object)
      } else {
        return(object)
      }
    }
  } else if ( inherits(object,"dgp") ) {
    if (vecchia){
      if (object$specs$vecchia) {
        return(object)
      } else {
        object$constructor_obj$to_vecchia(m = M, ord_fun = ord_wrapper)
        object$emulator_obj$to_vecchia()
        object$container_obj$to_vecchia()
        object[['specs']][['vecchia']] <- TRUE
        object[['specs']][['M']] <- M
        return(object)
      }
    } else {
      if (object$specs$vecchia) {
        object$constructor_obj$remove_vecchia()
        object$emulator_obj$remove_vecchia()
        object$container_obj$remove_vecchia()
        object[['specs']][['vecchia']] <- FALSE
        return(object)
      } else {
        return(object)
      }
    }
  } else if ( inherits(object,"lgp") ) {
    tryCatch({
      object$emulator_obj$set_vecchia(mode = vecchia)
    }, error = function(e) {
      if(grepl("mode has a different shape as all_layer", e$message)) {
        stop("'vecchia' has a different shape as 'struc' supplied to lgp().", call. = FALSE)
      }
    })
    return(object)
  }
}


#' @title Set linked indices
#'
#' @description
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated and will be removed in the next release. The updated
#' [lgp()] function now offers a simpler, more efficient way to specify linked information
#' for (D)GP emulators.
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
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See lgp() for an example.
#' }
#' @md
#' @keywords internal
#' @export
set_linked_idx <- function(object, idx) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  lifecycle::deprecate_warn(
    when = "2.5.0",
    what = "set_linked_idx()",
    details = c(i = "The function will be removed in the next release.",
                i = "Please use the updated `lgp()` function to specify linked information for (D)GP emulators."
    )
  )

  idx_py <- linked_idx_r_to_py(idx)
  object[['container_obj']] <- object$container_obj$set_local_input(idx_py, TRUE)
  object[['specs']][['linked_idx']] <- if ( is.null(idx) ) FALSE else idx
  return(object)
}

#' @title Reset number of imputations for a DGP emulator
#'
#' @description This function resets the number of imputations for prediction from a DGP emulator.
#'
#' @param object an instance of the S3 class `dgp`.
#' @param B the number of imputations to produce predictions from `object`. Increase the value to improve imputation uncertainty quantification. Decrease the value to improve speed of prediction. Defaults to `5`.
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
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
  new_object[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
  new_object[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
  new_object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
  new_object[['specs']][['vecchia']] <- object[['specs']][['vecchia']]
  new_object[['specs']][['M']] <- object[['specs']][['M']]
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


#' @title Trim the sequence of hyperparameter estimates within a DGP emulator
#'
#' @description This function trims the sequence of hyperparameter estimates within a DGP emulator
#'     generated during training.
#'
#' @param object an instance of the S3 class `dgp`.
#' @param start the first iteration before which all iterations are trimmed from the sequence.
#' @param end the last iteration after which all iterations are trimmed from the sequence.
#'     Set to `NULL` to keep all iterations after (including) `start`. Defaults to `NULL`.
#' @param thin the interval between the `start` and `end` iterations to thin out the sequence.
#'     Defaults to 1.
#'
#' @return An updated `object` with a trimmed sequence of hyperparameters.
#'
#' @note
#' * This function is useful when a DGP emulator has been trained and one wants to trim
#'   the sequence of hyperparameters estimated and to use the trimmed sequence to generate point estimates
#'   of the DGP model parameters for prediction.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()]
#'   in `object` will be removed and not contained in the returned object.
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
  new_object[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
  new_object[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
  new_object[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
  new_object[['specs']][['vecchia']] <- object[['specs']][['vecchia']]
  new_object[['specs']][['M']] <- object[['specs']][['M']]
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


#' @title Calculate the predictive negative log-likelihood
#'
#' @description This function computes the predictive negative log-likelihood from a
#'     DGP emulator with a likelihood layer.
#'
#' @param object an instance of the `dgp` class and it should be produced by [dgp()] with `likelihood` not being `NULL`;
#' @param x a matrix where each row is an input testing data point and each column is an input dimension.
#' @param y a matrix with only one column where each row is a scalar-valued testing output data point.
#'
#' @return An updated `object` is returned with an additional slot named `NLL` that contains two elements.
#'     The first one, named `meanNLL`, is a scalar that gives the average negative predicted log-likelihood
#'     across all testing data points. The second one, named `allNLL`, is a vector that gives the negative predicted
#'     log-likelihood for each testing data point.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
  if ( is.vector(x) ) {
    if ( ncol(object$data$X)!=1 ){
      x <- matrix(x, nrow = 1)
    } else {
      x <- as.matrix(x)
    }
  }
  if ( is.vector(y) ) {
    if ( ncol(object$data$Y)!=1 ){
      y <- matrix(y, nrow = 1)
    } else {
      y <- as.matrix(y)
    }
  }

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


#' @title Trace plot for DGP hyperparameters
#'
#' @description This function draws trace plots for the hyperparameters of a chosen GP node
#'     in a DGP emulator.
#'
#' @param object an instance of the `dgp` class.
#' @param layer the index of a layer. Defaults to `NULL` for the final layer.
#' @param node the index of a GP node in the layer specified by `layer`. Defaults to `1` for the first GP node in the
#'     corresponding layer.
#'
#' @return A `ggplot` object.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
#' @description This function implements static pruning for a DGP emulator.
#'
#' @param object an instance of the `dgp` class that is generated by `dgp()`.
#' @param control a list that can supply the following two components to control static pruning of the DGP emulator:
#' * `min_size`, the minimum number of design points required to trigger pruning. Defaults to 10 times of the input dimensions.
#' * `threshold`, the \eqn{R^2} value above which a GP node is considered redundant and removable. Defaults to `0.97`.
#' @param verb a bool indicating if trace information will be printed during the function execution. Defaults to `TRUE`.
#'
#' @return An updated `object` that could be an instance of `gp`, `dgp`, or `bundle` (of GP emulators) class.
#'
#' @note
#' * The function requires a DGP emulator that has been trained with a dataset comprising a minimum size equal to `min_size` in `control`.
#'    If the training dataset size is smaller than this, it is recommended that the design of the DGP emulator is enriched and its
#'    structure pruned dynamically using the `design()` function. Depending on the design of the DGP emulator, static pruning may not be accurate.
#'    It is thus recommended that dynamic pruning is implemented as a part of a sequential design via `design()`.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()];
#'
#'   in `object` will be removed and not contained in the returned object.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
    stop("To prune, 'object' needs to be trained with a dataset comprising a size at least equal to 'min_size' in 'control'. Use design() to enrich the training set.", call. = FALSE)
  }

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
  vecchia <- object[['specs']][['vecchia']]
  M <- object[['specs']][['M']]
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
          obj <- pkg.env$dgpsi$gp(X, Y, struc, vecchia, M)
          with(pkg.env$np$errstate(divide = 'ignore'), obj$train())
          res <- list()
          res[['id']] <- object$id
          res[['data']][['X']] <- X
          res[['data']][['Y']] <- Y
          res[['specs']] <- extract_specs(obj, "gp")
          res[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
          res[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
          res[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
          res[['specs']][['vecchia']] <- vecchia
          res[['specs']][['M']] <- M
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
            obj <- pkg.env$dgpsi$gp(X, Y[,j,drop=F], struc, vecchia, M)
            with(pkg.env$np$errstate(divide = 'ignore'), obj$train())
            res_j <- list()
            res_j[['id']] <- uuid::UUIDgenerate()
            res_j[['data']][['X']] <- X
            res_j[['data']][['Y']] <- Y[,j,drop=F]
            res_j[['specs']] <- extract_specs(obj, "gp")
            res_j[['specs']][['internal_dims']] <- object[['specs']][['internal_dims']]
            res_j[['specs']][['external_dims']] <- object[['specs']][['external_dims']]
            res_j[['specs']][['linked_idx']] <- if ( is.null(linked_idx) ) FALSE else linked_idx_py_to_r(linked_idx)
            res_j[['specs']][['vecchia']] <- vecchia
            res_j[['specs']][['M']] <- M
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
        object[['specs']][['vecchia']] <- vecchia
        object[['specs']][['M']] <- M
        id <- sample.int(100000, 1)
        set_seed(id)
        object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = blocked_gibbs)
        object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, blocked_gibbs)
        object[['specs']][['seed']] <- id
        object[['specs']][['B']] <- B
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
  object[['specs']][['vecchia']] <- vecchia
  object[['specs']][['M']] <- M
  id <- sample.int(100000, 1)
  set_seed(id)
  object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = blocked_gibbs)
  object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx, blocked_gibbs)
  object[['specs']][['seed']] <- id
  object[['specs']][['B']] <- B
  if ( verb ) message(" done")
  return(object)
}

#Badge function to render customised badges

new_badge <- function(stage) {

  url <- paste0("https://lifecycle.r-lib.org/articles/stages.html#", stage)
  html <- sprintf(
    "\\href{%s}{\\figure{%s}{options: alt='[%s]'}}",
    url,
    file.path(sprintf("lifecycle-%s.svg", stage)),
    upcase2(stage)
  )
  text <- sprintf("\\strong{[%s]}", upcase2(stage))

  sprintf("\\ifelse{html}{%s}{%s}", html, text)
}

upcase2 <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

get_docs_url <- function() {
  pkg_version <- as.character(utils::packageVersion("dgpsi"))

  is_dev <- grepl("\\.9000$", pkg_version)

  if (is_dev) {
    "https://mingdeyu.github.io/dgpsi-R/dev/"
  } else {
    "https://mingdeyu.github.io/dgpsi-R/"
  }
}

