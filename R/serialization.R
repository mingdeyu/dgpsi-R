#' @title Serialize the constructed emulator
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function serializes the constructed emulator.
#'
#' @param object an instance of the S3 class `gp`, `dgp`, `lgp`, or `bundle`.
#' @param light a bool indicating if a light version of the constructed emulator (that requires a small storage) will be serialized.
#'     Defaults to `TRUE`.
#'
#' @return A serialized version of `object`.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @note Since the constructed emulators are 'python' objects, they cannot be directly exported to other R processes for parallel
#'    processing. This function provides a solution by converting the emulators into serialized objects, which can be restored
#'    using [deserialize()] for multi-process parallel implementation.
#' @examples
#' \dontrun{
#'
#' library(parallel)
#' library(dgpsi)
#'
#' # model
#' f <- function(x) {
#'  (sin(7.5*x)+1)/2
#' }
#'
#' # training data
#' X <- seq(0, 1, length = 10)
#' Y <- sapply(X, f)
#'
#' # train a DGP emulator
#' m <- dgp(X, Y, name = "matern2.5")
#'
#' # testing input data
#' X_dgp <- seq(0, 1, length = 100)
#'
#' # serialize the DGP emulator
#' m_serialized <- serialize(m)
#'
#' # create a cluster with 8 workers for parallel predictions
#' cl <- makeCluster(8)
#'
#' # export objects to the cluster
#' clusterExport(cl, varlist = c("m_serialized", "X_dgp"))
#'
#' # initialize deserialized object on each worker
#' res <- clusterEvalQ(cl, {
#'   library(dgpsi)
#'   assign("m_deserialized", deserialize(m_serialized), envir = .GlobalEnv)
#' })
#'
#' # perform parallel predictions
#' results <- parLapply(cl, 1:length(X_dgp), function(i) {
#'   mean_i <- predict(m_deserialized, X_dgp[i])$results$mean
#' })
#'
#' # reset the cluster
#' stopCluster(cl)
#'
#' # combine mean predictions
#' pred_mean <- do.call(rbind, results)
#' }
#' @md
#' @export
serialize <- function(object, light = TRUE) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if (light) {
    if (inherits(object,"gp")){
      if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
      object[['container_obj']] <- NULL
    } else if (inherits(object,"dgp")){
      if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
      if ( !"seed" %in% names(object$specs) ) stop("The supplied 'object' cannot be serialized in light mode. To serialize, either set 'light = FALSE' or produce a new version of 'object' by set_imp().", call. = FALSE)
      object[['emulator_obj']] <- NULL
      object[['container_obj']] <- NULL
    } else if (inherits(object,"lgp")){
      if ( "metadata" %in% names(object$specs) ){
        if ( !("emulator_obj" %in% names(object)) ){
          stop("'object' must be activated for serialization. Please set `activate = TRUE` in `lgp()` to activate the emulator.", call. = FALSE)
        }
      }
      if ( reticulate::py_is_null_xptr(object$emulator_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
      if ( !"seed" %in% names(object$specs) ) stop("The supplied 'object' cannot be serialized in light mode. To serialize, either set 'light = FALSE' or re-construct and activate the 'object' by lgp().", call. = FALSE)
      object[['emulator_obj']] <- NULL
    } else if (inherits(object,"bundle")){
      if ( reticulate::py_is_null_xptr(object$emulator1$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulators in the bundle or, if the bundle was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
      N <- length(object) - 1
      if ( "id" %in% names(object) ) N <- N - 1
      if ( "design" %in% names(object) ) N <- N - 1
      for ( i in 1:N ){
        if ( inherits(object[[paste('emulator',i, sep='')]],"dgp") ) {
          if ( !"seed" %in% names(object[[paste('emulator',i, sep='')]][['specs']]) ) stop("The supplied 'object' cannot be serialized in light mode. To serialize, either set 'light = FALSE' or produce a new version of 'object' by updating the included DGP emulators via set_imp().", call. = FALSE)
          object[[paste('emulator',i, sep='')]][['emulator_obj']] <- NULL
          object[[paste('emulator',i, sep='')]][['container_obj']] <- NULL
        } else {
          object[[paste('emulator',i, sep='')]][['container_obj']] <- NULL
        }
      }
    }
  } else {
    if (inherits(object,"bundle")){
      if ( reticulate::py_is_null_xptr(object$emulator1$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulators in the bundle or, if the bundle was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
    } else if (inherits(object,"lgp")) {
      if ( "metadata" %in% names(object$specs) ){
        if ( !("emulator_obj" %in% names(object)) ){
          stop("'object' must be activated for serialization. Please set `activate = TRUE` in `lgp()` to activate the emulator.", call. = FALSE)
        }
      }
      if ( reticulate::py_is_null_xptr(object$emulator_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
    } else {
      if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
    }
  }
  label <- class(object)
  if (inherits(object,"bundle")) {
    X_train_name <- colnames(object$data$X[[1]])
    Y_train_name <- colnames(object$data$Y[[1]])
  } else if (inherits(object,"gp") || inherits(object,"dgp")) {
    X_train_name <- colnames(object$data$X)
    Y_train_name <- colnames(object$data$Y)
  } else {
    X_train_name <- NULL
    Y_train_name <- NULL
  }

  if (is.null(X_train_name)){
    X_train_name <- 'None'
  }

  if (is.null(Y_train_name)){
    Y_train_name <- 'None'
  }

  if (inherits(object, "dgp")){
    L = object$constructor_obj$n_layer
    final_node <- object$specs[[paste('layer', L, sep="")]][['node1']]
    if ("type" %in% names(final_node) && final_node$type == "Categorical") {
      class_name <- as.character(object$constructor_obj$all_layer[[L]][[1]]$class_encoder$classes_)
    } else {
      class_name <- NULL
    }
  } else {
    class_name <- NULL
  }

  if (is.null(class_name)){
    class_name <- 'None'
  }

  lst <- unclass(object)
  lst[['label']] <- label
  lst[['X_train_name']] <- X_train_name
  lst[['Y_train_name']] <- Y_train_name
  lst[['class_name']] <- class_name
  serialized_binary <- pkg.env$dill$dumps(lst)

  # Encode the binary string as a Base64 string
  serialized_obj <- pkg.env$base64$b64encode(serialized_binary)$decode("utf-8")
  return(serialized_obj)
}


#' @title Restore the serialized emulator
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function restores the serialized emulator created by [serialize()].
#'
#' @param object the serialized object of an emulator.
#'
#' @return The S3 class of a GP emulator, a DGP emulator, a linked (D)GP emulator, or a bundle of (D)GP emulators.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @note See the *Note* section in [serialize()].
#' @examples
#' \dontrun{
#'
#' # See serialize() for an example.
#' }
#' @md
#' @export
deserialize <- function(object) {
  # in next release remove the loading of linked_idx as this will be no longer needed.
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  serialized_binary <- pkg.env$base64$b64decode(object)
  res <- pkg.env$dill$loads(serialized_binary)
  if ('label' %in% names(res)){
    if ('X_train_name' %in% names(res)){
      X_train_name <- if (identical(res$X_train_name, "None")) NULL else res$X_train_name
      Y_train_name <- if (identical(res$Y_train_name, "None")) NULL else res$Y_train_name
      class_name <- if (identical(res$class_name, "None")) NULL else res$class_name
      res$X_train_name <- NULL
      res$Y_train_name <- NULL
      res$class_name <- NULL
    } else {
      X_train_name <- NULL
      Y_train_name <- NULL
      class_name <- NULL
    }
    label <- res$label
    res$label <- NULL
    if (label == "gp"){
      if ('container_obj' %in% names(res)){
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        colnames(res[['data']][['X']]) <- X_train_name
        colnames(res[['data']][['Y']]) <- Y_train_name
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
        colnames(res[['data']][['X']]) <- X_train_name
        colnames(res[['data']][['Y']]) <- Y_train_name
        class(res) <- "gp"
      }
    } else if (label == "dgp"){
      if ('emulator_obj' %in% names(res)){
        if (!'id' %in% names(res)) res[['id']] <- uuid::UUIDgenerate()
        if (!'vecchia' %in% names(res$specs)) {
          res[['specs']][['vecchia']] <- FALSE
          res[['specs']][['M']] <- 25
        }
        colnames(res[['data']][['X']]) <- X_train_name
        colnames(res[['data']][['Y']]) <- Y_train_name
        if (!is.null(class_name)){
          if ("loo" %in% names(res)) {
            colnames(res$loo$probability) <- class_name
          }
          if ("oos" %in% names(res)){
            colnames(res$oos$probability) <- class_name
          }
          if ("results" %in% names(res)){
            if (is.list(res$results$mean)){
              colnames(res[['results']][['mean']][[paste0('layer',length(res[['results']][['mean']]))]]) <- class_name
              colnames(res[['results']][['var']][[paste0('layer',length(res[['results']][['mean']]))]]) <- class_name
            } else {
              colnames(res$results$mean) <- class_name
              colnames(res$results$var) <- class_name
            }
          }
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
        colnames(res[['data']][['X']]) <- X_train_name
        colnames(res[['data']][['Y']]) <- Y_train_name
        if (!is.null(class_name)){
          if ("loo" %in% names(res)) {
            colnames(res$loo$probability) <- class_name
          }
          if ("oos" %in% names(res)){
            colnames(res$oos$probability) <- class_name
          }
          if ("results" %in% names(res)){
            if (is.list(res$results$mean)){
              colnames(res[['results']][['mean']][[paste0('layer',length(res[['results']][['mean']]))]]) <- class_name
              colnames(res[['results']][['var']][[paste0('layer',length(res[['results']][['mean']]))]]) <- class_name
            } else {
              colnames(res$results$mean) <- class_name
              colnames(res$results$var) <- class_name
            }
          }
        }
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
        colnames(res[['data']][['X']][[paste('emulator',i, sep='')]]) <- X_train_name
        colnames(res[['data']][['Y']][[paste('emulator',i, sep='')]]) <- Y_train_name
        if ('emulator_obj' %in% names(res[[paste('emulator',i, sep='')]])) {
          type <- pkg.env$py_buildin$type(res[[paste('emulator',i, sep='')]]$emulator_obj)$'__name__'
          if ( type=='emulator' ) {
            colnames(res[[paste('emulator',i, sep='')]][['data']][['X']]) <- X_train_name
            colnames(res[[paste('emulator',i, sep='')]][['data']][['Y']]) <- Y_train_name
            class(res[[paste('emulator',i, sep='')]]) <- "dgp"
          } else if ( type=='gp' ) {
            if ('container_obj' %in% names(res[[paste('emulator',i, sep='')]])){
              colnames(res[[paste('emulator',i, sep='')]][['data']][['X']]) <- X_train_name
              colnames(res[[paste('emulator',i, sep='')]][['data']][['Y']]) <- Y_train_name
              class(res[[paste('emulator',i, sep='')]]) <- "gp"
            } else {
              est_obj <- res[[paste('emulator',i, sep='')]]$constructor_obj$export()
              if (is.null(res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']])){
                linked_idx <- NULL
              } else {
                linked_idx <- if ( isFALSE( res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]) ) {NULL} else {res[[paste('emulator',i, sep='')]][['specs']][['linked_idx']]}
              }
              res[[paste('emulator',i, sep='')]][['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx_r_to_py(linked_idx))
              colnames(res[[paste('emulator',i, sep='')]][['data']][['X']]) <- X_train_name
              colnames(res[[paste('emulator',i, sep='')]][['data']][['Y']]) <- Y_train_name
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
          colnames(res[[paste('emulator',i, sep='')]][['data']][['X']]) <- X_train_name
          colnames(res[[paste('emulator',i, sep='')]][['data']][['Y']]) <- Y_train_name
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
