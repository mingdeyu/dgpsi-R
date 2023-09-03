#' @title Validate a constructed GP, DGP, or linked (D)GP emulator
#'
#' @description This function validate a constructed GP, DGP, or linked (D)GP emulator via the Leave-One-Out (LOO)
#'    cross validation or Out-Of-Sample (OOS) validation.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param x_test the OOS testing input data:
#' * if `x` is an instance of the `gp` or `dgp` class, `x_test` is a matrix where each row is an input testing data point and each column is an input dimension.
#' * if `x` is an instance of the `lgp` class, `x_test` can be a matrix or a list:
#'    - if `x_test` is a matrix, it is the global testing input data that feed into the emulators in the first layer of a system.
#'      The rows of `x_test` represent different input data points and the columns represent input dimensions across all emulators in
#'      the first layer of the system. In this case, it is assumed that the only global input to the system is the input to the
#'      emulators in the first layer and there is no global input to emulators in other layers.
#'    - if `x_test` is a list, it should have *L* (the number of layers in an emulator system) elements. The first element
#'      is a matrix that represents the global testing input data that feed into the emulators in the first layer of the system. The
#'      remaining *L-1* elements are *L-1* sub-lists, each of which contains a number (the same number of emulators in
#'      the corresponding layer) of matrices (rows being testing input data points and columns being input dimensions) that represent the
#'      global testing input data to the emulators in the corresponding layer. The matrices must be placed in the sub-lists based on how
#'      their corresponding emulators are placed in `struc` argument of [lgp()]. If there is no global input data to a certain emulator,
#'      set `NULL` in the corresponding sub-list of `x_test`.
#'
#' `x_test` must be provided for the validation if `x` is an instance of the `lgp`. Defaults to `NULL`.
#' @param y_test the OOS testing output data that correspond to `x_test`:
#' * if `x` is an instance of the `gp` class, `y_test` is a matrix with only one column and each row being an testing output data point.
#' * if `x` is an instance of the `dgp` class, `y_test` is a matrix with its rows being testing output data points and columns being
#'   output dimensions.
#' * if `x` is an instance of the `lgp` class, `y_test` can be a single matrix or a list of matrices:
#'   - if `y_test` is a single matrix, then there is only one emulator in the final layer of the linked emulator system and `y_test`
#'     represents the emulator's output with rows being testing positions and columns being output dimensions.
#'   - if `y_test` is a list, then `y_test` should have *M* number (the same number of emulators in the final layer of the system) of matrices.
#'     Each matrix has its rows corresponding to testing positions and columns corresponding to output dimensions of the associated emulator
#'     in the final layer.
#'
#' `y_test` must be provided for the validation if `x` is an instance of the `lgp`. Defaults to `NULL`.
#' @param method the prediction approach in validations: mean-variance (`"mean_var"`) or sampling (`"sampling"`) approach. Defaults to `"mean_var"`.
#' @param verb a bool indicating if the trace information on validations will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param force a bool indicating whether to force the LOO or OOS re-evaluation when `loo` or `oos` slot already exists in `object`. When `force = FALSE`,
#'     [validate()] will try to determine automatically if the LOO or OOS re-evaluation is needed. Set `force` to `TRUE` when LOO or OOS re-evaluation
#'     is required. Defaults to `FALSE`.
#' @param cores the number of cores/workers to be used for the LOO or OOS validation. If set to `NULL`,
#'     the number of cores is set to `(max physical cores available - 1)`. Defaults to `1`.
#' @param threading a bool indicating whether to use the multi-threading to accelerate the LOO or OOS.
#'     Turning this option on could improve the speed of validations when the emulator is built with a moderately large number of
#'     training data points and the Matérn-2.5 kernel.
#' @param ... N/A.
#'
#' @return
#' * If `object` is an instance of the `gp` class, an updated `object` is returned with an additional slot called `loo` (for LOO cross validation) or
#'   `oos` (for OOS validation) that contains:
#'   - two slots called `x_train` (or `x_test`) and `y_train` (or `y_test`) that contain the validation data points for LOO (or OOS).
#'   - a column matrix called `mean`, if `method = "mean_var"`, or `median`, if `method = "sampling"`, that contains the predictive means or medians of the
#'     GP emulator at validation positions.
#'   - three column matrices called `std`, `lower`, and `upper` that contain the predictive standard deviations and credible intervals of the
#'     GP emulator at validation positions. If `method = "mean_var"`, the upper and lower bounds of a credible interval are two standard deviations above
#'     and below the predictive mean. If `method = "sampling"`, the upper and lower bounds of a credible interval are 2.5th and 97.5th percentiles.
#'   - a numeric value called `rmse` that contains the root mean/median squared error of the GP emulator.
#'   - a numeric value called `nrmse` that contains the (min-max) normalized root mean/median squared error of the GP emulator. The min-max normalization
#'     is based on the maximum and minimum values of the validation outputs contained in `y_train` (or `y_test`).
#'
#'   The rows of matrices (`mean`, `median`, `std`, `lower`, and `upper`) correspond to the validation positions.
#' * If `object` is an instance of the `dgp` class, an updated `object` is returned with an additional slot called `loo` (for LOO cross validation) or
#'   `oos` (for OOS validation) that contains:
#'   - two slots called `x_train` (or `x_test`) and `y_train` (or `y_test`) that contain the validation data points for LOO (or OOS).
#'   - a matrix called `mean`, if `method = "mean_var"`, or `median`, if `method = "sampling"`, that contains the predictive means or medians of the
#'     DGP emulator at validation positions.
#'   - three matrices called `std`, `lower`, and `upper` that contain the predictive standard deviations and credible intervals of the
#'     DGP emulator at validation positions. If `method = "mean_var"`, the upper and lower bounds of a credible interval are two standard deviations above
#'     and below the predictive mean. If `method = "sampling"`, the upper and lower bounds of a credible interval are 2.5th and 97.5th percentiles.
#'   - a vector called `rmse` that contains the root mean/median squared errors of the DGP emulator across different output
#'     dimensions.
#'   - a vector called `nrmse` that contains the (min-max) normalized root mean/median squared errors of the DGP emulator across different output
#'     dimensions. The min-max normalization is based on the maximum and minimum values of the validation outputs contained in `y_train` (or `y_test`).
#'
#'   The rows and columns of matrices (`mean`, `median`, `std`, `lower`, and `upper`) correspond to the validation positions and DGP emulator output
#' dimensions, respectively.
#' * If `object` is an instance of the `lgp` class, an updated `object` is returned with an additional slot called `oos` (for OOS validation) that contains:
#'   - two slots called `x_test` and `y_test` that contain the validation data points for OOS.
#'   - a list called `mean`, if `method = "mean_var"`, or `median`, if `method = "sampling"`, that contains the predictive means or medians of
#'     the linked (D)GP emulator at validation positions.
#'   - three lists called `std`, `lower`, and `upper` that contain the predictive standard deviations and credible intervals of
#'     the linked (D)GP emulator at validation positions. If `method = "mean_var"`, the upper and lower bounds of a credible interval are two standard
#'     deviations above and below the predictive mean. If `method = "sampling"`, the upper and lower bounds of a credible interval are 2.5th and 97.5th percentiles.
#'   - a list called `rmse` that contains the root mean/median squared errors of the linked (D)GP emulator.
#'   - a list called `nrmse` that contains the (min-max) normalized root mean/median squared errors of the linked (D)GP emulator. The min-max normalization
#'     is based on the maximum and minimum values of the validation outputs contained in `y_test`.
#'
#'   Each element in `mean`, `median`, `std`, `lower`, `upper`, `rmse`, and `nrmse` corresponds to a (D)GP emulator in the final layer of the linked (D)GP
#' emulator.
#'
#' @note
#' * When both `x_test` and `y_test` are `NULL`, the LOO cross validation will be implemented. Otherwise, OOS validation will
#'   be implemented. The LOO validation is only applicable to a GP or DGP emulator (i.e., `x` is an instance of the `gp` or `dgp`
#'   class). If a linked (D)GP emulator (i.e., `x` is an instance of the `lgp` class) is provided, `x_test` and `y_test` must
#'   also be provided for OOS validation.
#' * Any R vector detected in `x_test` and `y_test` will be treated as a column vector and automatically converted into a single-column
#'   R matrix. Thus, if `x_test` or `y_test` is a single testing data point with multiple dimensions, it must be given as a matrix.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name validate
#' @export
validate <- function(object, x_test, y_test, method, verb, force, cores, ...){
  UseMethod("validate")
}

#' @rdname validate
#' @method validate gp
#' @export
validate.gp <- function(object, x_test = NULL, y_test = NULL, method = 'mean_var', verb = TRUE, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)

  #For LOO
  if (is.null(x_test) & is.null(y_test)){
    #check existing LOO
    if ( isFALSE(force) ){
      if ( "loo" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" LOO results found in the gp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( ((method == 'mean_var')&("mean" %in% names(object$loo)))|((method == 'sampling')&("median" %in% names(object$loo))) ){
          if ( isTRUE(verb) ) message(" LOO re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting gp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else {
          if ( isTRUE(verb) ) message(" LOO re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the LOO ...", appendLF = FALSE)
    x_train <- object$constructor_obj$X
    y_train <- object$constructor_obj$Y
    dat <- list('x_train' = x_train,'y_train' = y_train)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the LOO ...", appendLF = FALSE)
    res <- object$emulator_obj$loo(method = method, sample_size = 500L)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Saving results to the slot 'loo' in the gp object ...", appendLF = FALSE)
    if ( method=='mean_var' ){
      dat[["mean"]] <- res[[1]]
      dat[["std"]] <- sqrt(res[[2]])
      dat[["lower"]] <- dat$mean-2*dat$std
      dat[["upper"]] <- dat$mean+2*dat$std
      dat[["rmse"]] <- sqrt(mean((dat$mean-dat$y_train)^2))
      dat[["nrmse"]] <- dat$rmse/(max(dat$y_train)-min(dat$y_train))
    } else if ( method=='sampling' ){
      quant <- t(pkg.env$np$quantile(res, c(0.025, 0.5, 0.975), axis=1L))
      std <- pkg.env$np$std(res, axis=1L, keepdims=TRUE)
      dat[["median"]] <- quant[,2,drop=F]
      dat[["std"]] <- std
      dat[["lower"]] <- quant[,1,drop=F]
      dat[["upper"]] <- quant[,3,drop=F]
      dat[["rmse"]] <- sqrt(mean((dat$median-dat$y_train)^2))
      dat[["nrmse"]] <- dat$rmse/(max(dat$y_train)-min(dat$y_train))
    }
    object$loo <- dat
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    return(object)
    #For OOS
  } else if (!is.null(x_test) & !is.null(y_test)) {
    x_test <- unname(x_test)
    y_test <- unname(y_test)
    if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
    if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)

    #check existing OOS
    if ( isFALSE(force) ){
      if ( "oos" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" OOS results found in the gp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (((method == 'mean_var')&("mean" %in% names(object$oos)))|((method == 'sampling')&("median" %in% names(object$oos)))) ){
          if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting gp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else {
          if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the OOS ...", appendLF = FALSE)
    dat <- list('x_test' = x_test,'y_test' = y_test)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the OOS ...", appendLF = FALSE)
    rep_x <- pkg.env$np$unique(x_test, return_inverse=TRUE, axis=0L)
    x_test_unique <- rep_x[[1]]
    rep <- rep_x[[2]] + 1

    if ( identical(cores,as.integer(1)) ){
      res <- object$emulator_obj$predict(x_test_unique, method = method, sample_size=500L)
    } else {
      res <- object$emulator_obj$ppredict(x_test_unique, method = method, sample_size=500L, core_num = cores)
    }
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Saving results to the slot 'oos' in the gp object ...", appendLF = FALSE)
    if ( method == 'mean_var' ){
      dat[["mean"]] <- res[[1]][rep,,drop=FALSE]
      dat[["std"]] <- sqrt(res[[2]][rep,,drop=FALSE])
      dat[["lower"]] <- dat$mean-2*dat$std
      dat[["upper"]] <- dat$mean+2*dat$std
      dat[["rmse"]] <- sqrt(mean((dat$mean-dat$y_test)^2))
      dat[["nrmse"]] <- dat$rmse/(max(dat$y_test)-min(dat$y_test))
    } else if ( method == 'sampling' ){
      quant <- t(pkg.env$np$quantile(res, c(0.025, 0.5, 0.975), axis=1L))[rep,,drop=FALSE]
      std <- pkg.env$np$std(res, axis=1L, keepdims=TRUE)[rep,,drop=FALSE]
      dat[["median"]] <- quant[,2,drop=F]
      dat[["std"]] <- std
      dat[["lower"]] <- quant[,1,drop=F]
      dat[["upper"]] <- quant[,3,drop=F]
      dat[["rmse"]] <- sqrt(mean((dat$median-dat$y_test)^2))
      dat[["nrmse"]] <- dat$rmse/(max(dat$y_test)-min(dat$y_test))
    }
    object$oos <- dat
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    return(object)
    #For other cases
  } else {
    stop("Either 'x_test' or 'y_test' is not given.", call. = FALSE)
  }
}

#' @rdname validate
#' @method validate dgp
#' @export
validate.dgp <- function(object, x_test = NULL, y_test = NULL, method = 'mean_var', verb = TRUE, force = FALSE, cores = 1, threading = FALSE, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }
  object$emulator_obj$set_nb_parallel(threading)

  if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)

  #For LOO
  if (is.null(x_test) & is.null(y_test)){
    #check existing LOO
    if ( isFALSE(force) ){
      if ( "loo" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" LOO results found in the dgp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( ((method == 'mean_var')&("mean" %in% names(object$loo)))|((method == 'sampling')&("median" %in% names(object$loo))) ){
          if ( isTRUE(verb) ) message(" LOO re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting dgp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else {
          if ( isTRUE(verb) ) message(" LOO re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the LOO ...", appendLF = FALSE)
    x_train <- object$constructor_obj$X
    y_train <- object$constructor_obj$Y
    rep <- object$constructor_obj$indices
    if ( !is.null(rep) ){
      rep <- rep + 1
      x_train <- x_train[rep,,drop=FALSE]
    }
    dat <- list('x_train' = x_train,'y_train' = y_train)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the LOO ...", appendLF = FALSE)
    if ( identical(cores,as.integer(1)) ){
      res <- object$emulator_obj$loo(X = x_train, method = method)
    } else {
      res <- object$emulator_obj$ploo(X = x_train, method = method, core_num = cores)
    }
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Saving results to the slot 'loo' in the dgp object ...", appendLF = FALSE)
    if ( method == 'sampling' ){
      res_np <- pkg.env$np$array(res)
      quant <- pkg.env$np$transpose(pkg.env$np$quantile(res_np, c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
      std <- pkg.env$np$std(res_np, axis=2L)
      dat[["median"]] <- as.matrix(quant[2,,])
      dat[["std"]] <- t(std)
      dat[["lower"]] <- as.matrix(quant[1,,])
      dat[["upper"]] <- as.matrix(quant[3,,])
      dat[["rmse"]] <- sqrt(colMeans((dat$median-dat$y_train)^2))
      dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_train, axis=0L)-pkg.env$np$amin(dat$y_train, axis=0L))
    } else if ( method == 'mean_var' ) {
      dat[["mean"]] <- res[[1]]
      dat[["std"]] <- sqrt(res[[2]])
      dat[["lower"]] <- dat$mean-2*dat$std
      dat[["upper"]] <- dat$mean+2*dat$std
      dat[["rmse"]] <- sqrt(colMeans((dat$mean-dat$y_train)^2))
      dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_train, axis=0L)-pkg.env$np$amin(dat$y_train, axis=0L))
    }
    object$loo <- dat
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    return(object)
    #For OOS
  } else if (!is.null(x_test) & !is.null(y_test)) {
    x_test <- unname(x_test)
    y_test <- unname(y_test)
    if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
    if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)

    #check existing OOS
    if ( isFALSE(force) ){
      if ( "oos" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" OOS results found in the dgp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (((method == 'mean_var')&("mean" %in% names(object$oos)))|((method == 'sampling')&("median" %in% names(object$oos)))) ){
          if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting dgp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else {
          if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the OOS ...", appendLF = FALSE)
    dat <- list('x_test' = x_test,'y_test' = y_test)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the OOS ...", appendLF = FALSE)
    rep_x <- pkg.env$np$unique(x_test, return_inverse=TRUE, axis=0L)
    x_test_unique <- rep_x[[1]]
    rep <- rep_x[[2]] + 1

    if ( identical(cores,as.integer(1)) ){
      res <- object$emulator_obj$predict(x = x_test_unique, method = method)
    } else {
      res <- object$emulator_obj$ppredict(x = x_test_unique, method = method, core_num = cores)
    }
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Saving results to the slot 'oos' in the dgp object ...", appendLF = FALSE)
    if ( method == 'sampling' ){
      res_np <- pkg.env$np$array(res)
      quant <- pkg.env$np$transpose(pkg.env$np$quantile(res_np, c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
      std <- pkg.env$np$std(res_np, axis=2L)
      dat[["median"]] <- as.matrix(quant[2,,])[rep,,drop=F]
      dat[["std"]] <- t(std)[rep,,drop=F]
      dat[["lower"]] <- as.matrix(quant[1,,])[rep,,drop=F]
      dat[["upper"]] <- as.matrix(quant[3,,])[rep,,drop=F]
      dat[["rmse"]] <- sqrt(colMeans((dat$median-dat$y_test)^2))
      dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_test, axis=0L)-pkg.env$np$amin(dat$y_test, axis=0L))
    } else if ( method == 'mean_var' ) {
      dat[["mean"]] <- res[[1]][rep,,drop=F]
      dat[["std"]] <- sqrt(res[[2]][rep,,drop=F])
      dat[["lower"]] <- dat$mean-2*dat$std
      dat[["upper"]] <- dat$mean+2*dat$std
      dat[["rmse"]] <- sqrt(colMeans((dat$mean-dat$y_test)^2))
      dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_test, axis=0L)-pkg.env$np$amin(dat$y_test, axis=0L))
    }
    object$oos <- dat
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    return(object)
    #For other cases
  } else {
    stop("Either 'x_test' or 'y_test' is not given.", call. = FALSE)
  }
}


#' @rdname validate
#' @method validate lgp
#' @export
validate.lgp <- function(object, x_test = NULL, y_test = NULL, method = 'mean_var', verb = TRUE, force = FALSE, cores = 1, threading = FALSE, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"lgp") ) stop("'object' must be an instance of the 'lgp' class.", call. = FALSE)
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }
  object$emulator_obj$set_nb_parallel(threading)

  if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)

  #For OOS
  if (!is.null(x_test) & !is.null(y_test)) {
    #check testing input
    if ( !is.list(x_test) ) {
      if ( !is.matrix(x_test)&!is.vector(x_test) ) {
        stop("'x_test' must be a vector or a matrix.", call. = FALSE)
      } else {
        x_test <- unname(x_test)
        if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
        nrow_x <- nrow(x_test)
      }
    } else {
      for ( l in 1:length(x_test) ){
        if ( l==1 ){
          if ( !is.matrix(x_test[[l]])&!is.vector(x_test[[l]]) ) {
            stop("The first element of 'x_test' must be a vector or a matrix.", call. = FALSE)
          } else {
            x_test[[l]] <- unname(x_test[[l]])
            if ( is.vector(x_test[[l]]) ) x_test[[l]] <- as.matrix(x_test[[l]])
            nrow_x <- nrow(x_test[[l]])
          }
        } else {
          for ( k in 1:length(x_test[[l]]) ){
            if ( !is.matrix(x_test[[l]][[k]])&!is.null(x_test[[l]][[k]])&!is.vector(x_test[[l]][[k]]) ) stop(sprintf("The element %i in the sublist %i of 'x_test' must be a vector, a matrix, or 'NULL'.", k, l), call. = FALSE)
            if ( is.matrix(x_test[[l]][[k]])|is.vector(x_test[[l]][[k]]) ){
              x_test[[l]][[k]] <- unname(x_test[[l]][[k]])
              if (is.vector(x_test[[l]][[k]])) x_test[[l]][[k]] <- as.matrix(x_test[[l]][[k]])
              if ( nrow(x_test[[l]][[k]])!=nrow_x ) {
                stop(sprintf("The element %i in the sublist %i of 'x_test' has inconsistent number of data points with the first element of 'x_test'.", k, l), call. = FALSE)
              }
            }
          }
        }
      }
    }

    #check testing output
    if ( !is.list(y_test) ) {
      if ( !is.matrix(y_test)&!is.vector(y_test) ) {
        stop("'y_test' must be a vector or a matrix.", call. = FALSE)
      } else {
        y_test <- unname(y_test)
        if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
        nrow_y <- nrow(y_test)
        if ( nrow_y!=nrow_x ) stop("The number of data points are inconsistent between 'x_test' and 'y_test'.", call. = FALSE)
        }
    } else {
      for ( l in 1:length(y_test) ){
          if ( !is.matrix(y_test[[l]])&!is.vector(y_test[[l]]) ) {
            stop(sprintf("The element %i of 'y_test' must be a vector or a matrix.", l), call. = FALSE)
          } else {
            y_test[[l]] <- unname(y_test[[l]])
            if ( is.vector(y_test[[l]]) ) y_test[[l]] <- as.matrix(y_test[[l]])
            nrow_y <- nrow(y_test[[l]])
            if ( nrow_y!=nrow_x ) stop(sprintf("The number of data points are inconsistent between 'x_test' and the element %i of 'y_test'.", l), call. = FALSE)
          }
      }
    }

    #check existing OOS
    if ( isFALSE(force) ){
      if ( "oos" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" OOS results found in the lgp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (((method == 'mean_var')&("mean" %in% names(object$oos)))|((method == 'sampling')&("median" %in% names(object$oos)))) ){
          if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting lgp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else {
          if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the OOS ...", appendLF = FALSE)
    dat <- list('x_test' = x_test,'y_test' = y_test)
    if ( !is.list(y_test) ) y_test <- list(y_test)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the OOS ...", appendLF = FALSE)
    if ( identical(cores,as.integer(1)) ){
      res <- object$emulator_obj$predict(x = x_test, method = method)
    } else {
      res <- object$emulator_obj$ppredict(x = x_test, method = method, core_num = cores)
    }
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Saving results to the slot 'oos' in the lgp object ...", appendLF = FALSE)
    if ( method == 'sampling' ){
      median_lst <- list()
      std_lst <- list()
      lower_lst <- list()
      upper_lst <- list()
      rmse_lst <- list()
      nrmse_lst <- list()
      for ( l in 1:length(res) ) {
        quant <- pkg.env$np$transpose(pkg.env$np$quantile(res[[l]], c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
        std <- pkg.env$np$std(res[[l]], axis=2L)
        median_lst[[l]] <- as.matrix(quant[2,,])
        std_lst[[l]] <- t(std)
        lower_lst[[l]] <- as.matrix(quant[1,,])
        upper_lst[[l]] <- as.matrix(quant[3,,])
        rmse_lst[[l]] <- sqrt(colMeans((median_lst[[l]]-y_test[[l]])^2))
        nrmse_lst[[l]] <- rmse_lst[[l]]/(pkg.env$np$amax(y_test[[l]], axis=0L)-pkg.env$np$amin(y_test[[l]], axis=0L))
      }
      dat[["median"]] <- median_lst
      dat[["std"]] <- std_lst
      dat[["lower"]] <- lower_lst
      dat[["upper"]] <- upper_lst
      dat[["rmse"]] <- rmse_lst
      dat[["nrmse"]] <- nrmse_lst
    } else if ( method == 'mean_var' ) {
      mean_lst <- list()
      std_lst <- list()
      lower_lst <- list()
      upper_lst <- list()
      rmse_lst <- list()
      nrmse_lst <- list()
      for ( l in 1:length(res[[1]]) ) {
        mean_lst[[l]] <- res[[1]][[l]]
        std_lst[[l]] <- sqrt(res[[2]][[l]])
        lower_lst[[l]] <- mean_lst[[l]]-2*std_lst[[l]]
        upper_lst[[l]] <- mean_lst[[l]]+2*std_lst[[l]]
        rmse_lst[[l]] <- sqrt(colMeans((mean_lst[[l]]-y_test[[l]])^2))
        nrmse_lst[[l]] <- rmse_lst[[l]]/(pkg.env$np$amax(y_test[[l]], axis=0L)-pkg.env$np$amin(y_test[[l]], axis=0L))
      }
      dat[["mean"]] <- mean_lst
      dat[["std"]] <- std_lst
      dat[["lower"]] <- lower_lst
      dat[["upper"]] <- upper_lst
      dat[["rmse"]] <- rmse_lst
      dat[["nrmse"]] <- nrmse_lst
    }
    object$oos <- dat
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    return(object)
    #For other cases
  } else {
    stop("Both 'x_test' and 'y_test' must be provided for the validation of a linked (D)GP emulator.", call. = FALSE)
  }
}
