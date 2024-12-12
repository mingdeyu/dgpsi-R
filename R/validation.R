#' @title Validate a constructed GP, DGP, or linked (D)GP emulator
#'
#' @description This function calculates Leave-One-Out (LOO) cross validation or Out-Of-Sample (OOS) validation statistics for a constructed GP, DGP, or linked (D)GP emulator.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param x_test OOS testing input data:
#' * if `object` is an instance of the `gp` or `dgp` class, `x_test` is a matrix where each row is a new input location to be used for validating the emulator and each column is an input dimension.
#' * `r lifecycle::badge("deprecated")` if `object` is an instance of the `lgp` class, `x_test` can be a matrix or a list:
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
#'    **This option for linked (D)GP emulators is deprecated and will be removed in the next release.**
#' * `r new_badge("new")` If `object` is an instance of the `lgp` class created by [lgp()] with argument `struc` in data frame form,
#'   `x_test` must be a matrix representing the global input, where each row corresponds to a test data point and each column represents a global input dimension.
#'   The column indices in `x_test` must align with the indices specified in the `From_Output` column of the `struc` data frame (used in [lgp()]),
#'   corresponding to rows where the `From_Emulator` column is `"Global"`.
#'
#' `x_test` must be provided if `object` is an instance of the `lgp`. `x_test` must also be provided if `y_test` is provided. Defaults to `NULL`, in which case LOO validation is performed.
#' @param y_test the OOS output data corresponding to `x_test`:
#' * if `object` is an instance of the `gp` class, `y_test` is a matrix with only one column where each row represents the output corresponding to the matching row of `x_test`.
#' * if `object` is an instance of the `dgp` class, `y_test` is a matrix where each row represents the output corresponding to the matching row of `x_test` and with columns representing output dimensions.
#' * if `object` is an instance of the `lgp` class, `y_test` can be a single matrix or a list of matrices:
#'   - if `y_test` is a single matrix, then there should be only one emulator in the final layer of the linked emulator system and `y_test`
#'     represents the emulator's output with rows being testing positions and columns being output dimensions.
#'   - if `y_test` is a list, then `y_test` should have *L* matrices, where *L* is the number of emulators in the final layer of the system.
#'     Each matrix has its rows corresponding to testing positions and columns corresponding to output dimensions of the associated emulator
#'     in the final layer.
#'
#' `y_test` must be provided if `object` is an instance of the `lgp`. `y_test` must also be provided if `x_test` is provided. Defaults to `NULL`, in which case LOO validation is performed.
#' @param method `r new_badge("updated")` the prediction approach to use for validation: either the mean-variance approach (`"mean_var"`) or the sampling approach (`"sampling"`). For details see [predict()].
#'      For DGP emulators with a categorical likelihood (`likelihood = "Categorical"` in [dgp()]), only the sampling approach is supported.
#'      By default, the method is set to `"sampling"` for DGP emulators with Poisson, Negative Binomial, and Categorical likelihoods and `"mean_var"` otherwise.
#' @param sample_size the number of samples to draw for each given imputation if `method = "sampling"`. Defaults to `50`.
#' @param verb a bool indicating if trace information for validation should be printed during function execution.
#'     Defaults to `TRUE`.
#' @param M `r new_badge("new")` the size of the conditioning set for the Vecchia approximation in emulator validation. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation. Defaults to `50`.
#' @param force a bool indicating whether to force LOO or OOS re-evaluation when the `loo` or `oos` slot already exists in `object`. When `force = FALSE`,
#'     [validate()] will only re-evaluate the emulators if the `x_test` and `y_test` are not identical to the values in the `oos` slot. If the existing `loo` or `oos` validation used a different `M` in a Vecchia approximation or a different `method` to the one prescribed in this call, the emulator will be re-evaluated. Set `force` to `TRUE` when LOO or OOS re-evaluation
#'     is required. Defaults to `FALSE`.
#' @param cores the number of processes to be used for validation. If set to `NULL`, the number of processes is set to `max physical cores available %/% 2`.
#'     Defaults to `1`.
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
#'   - a numeric value called `nrmse` that contains the (max-min) normalized root mean/median squared error of the GP emulator. The max-min normalization
#'     uses the maximum and minimum values of the validation outputs contained in `y_train` (or `y_test`).
#'   - `r new_badge("new")` an integer called `M` that contains the size of the conditioning set used for the Vecchia approximation, if used, for emulator validation.
#'   - an integer called `sample_size` that contains the number of samples used for validation if `method = "sampling"`.
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
#'   - a vector called `nrmse` that contains the (max-min) normalized root mean/median squared errors of the DGP emulator across different output
#'     dimensions. The max-min normalization uses the maximum and minimum values of the validation outputs contained in `y_train` (or `y_test`).
#'   - `r new_badge("new")` an integer called `M` that contains size of the conditioning set used for the Vecchia approximation, if used, for emulator validation.
#'   - an integer called `sample_size` that contains the number of samples used for validation if `method = "sampling"`.
#'
#'   The rows and columns of matrices (`mean`, `median`, `std`, `lower`, and `upper`) correspond to the validation positions and DGP emulator output
#' dimensions, respectively.
#' * `r new_badge("new")` If `object` is an instance of the `dgp` class with a categorical likelihood, an updated `object` is returned with an additional slot called `loo`
#'   (for LOO cross validation) or `oos` (for OOS validation) that contains:
#'   - two slots called `x_train` (or `x_test`) and `y_train` (or `y_test`) that contain the validation data points for LOO (or OOS).
#'   - a matrix called `label` that contains predictive samples of labels from the DGP emulator at validation positions. The matrix has its rows
#'     corresponding to validation positions and columns corresponding to samples of labels.
#'   - a list called `probability` that contains predictive samples of probabilities for each class from the DGP emulator at validation positions. The element in the list
#'     is a matrix that has its rows corresponding to validation positions and columns corresponding to samples of probabilities.
#'   - a scalar called `log_loss` that represents the average log loss of the predicted labels in the DGP emulator across all validation positions. Log loss measures the
#'     accuracy of probabilistic predictions, with lower values indicating better classification performance. `log_loss` ranges from `0` to positive infinity, where a
#'     value closer to `0` suggests more confident and accurate predictions.
#'   - an integer called `M` that contains size of the conditioning set used for the Vecchia approximation, if used, in emulator validation.
#'   - an integer called `sample_size` that contains the number of samples used for validation.
#' * If `object` is an instance of the `lgp` class, an updated `object` is returned with an additional slot called `oos` (for OOS validation) that contains:
#'   - two slots called `x_test` and `y_test` that contain the validation data points for OOS.
#'   - a list called `mean`, if `method = "mean_var"`, or `median`, if `method = "sampling"`, that contains the predictive means or medians of
#'     the linked (D)GP emulator at validation positions.
#'   - three lists called `std`, `lower`, and `upper` that contain the predictive standard deviations and credible intervals of
#'     the linked (D)GP emulator at validation positions. If `method = "mean_var"`, the upper and lower bounds of a credible interval are two standard
#'     deviations above and below the predictive mean. If `method = "sampling"`, the upper and lower bounds of a credible interval are 2.5th and 97.5th percentiles.
#'   - a list called `rmse` that contains the root mean/median squared errors of the linked (D)GP emulator.
#'   - a list called `nrmse` that contains the (max-min) normalized root mean/median squared errors of the linked (D)GP emulator. The max-min normalization
#'     uses the maximum and minimum values of the validation outputs contained in `y_test`.
#'   - `r new_badge("new")` an integer called `M` that contains size of the conditioning set used for the Vecchia approximation, if used, in emulator validation.
#'   - an integer called `sample_size` that contains the number of samples used for validation if `method = "sampling"`.
#'
#'   Each element in `mean`, `median`, `std`, `lower`, `upper`, `rmse`, and `nrmse` corresponds to a (D)GP emulator in the final layer of the linked (D)GP
#' emulator.
#'
#' @note
#' * When both `x_test` and `y_test` are `NULL`, LOO cross validation will be implemented. Otherwise, OOS validation will
#'   be implemented. LOO validation is only applicable to a GP or DGP emulator (i.e., `object` is an instance of the `gp` or `dgp`
#'   class). If a linked (D)GP emulator (i.e., `object` is an instance of the `lgp` class) is provided, `x_test` and `y_test` must
#'   also be provided for OOS validation.
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name validate
#' @export
validate <- function(object, x_test, y_test, method, sample_size, verb, M, force, cores, ...){
  UseMethod("validate")
}

#' @rdname validate
#' @method validate gp
#' @export
validate.gp <- function(object, x_test = NULL, y_test = NULL, method = NULL, sample_size = 50, verb = TRUE, M = 50, force = FALSE, cores = 1, ...) {
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
  M <- as.integer(M)

  if ( is.null(method) ){
    method = 'mean_var'
  } else {
    if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)
  }
  sample_size <- as.integer(sample_size)
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
        if ( (method == 'mean_var')&("mean" %in% names(object$loo)) & (M == object$loo$M) ){
          if ( isTRUE(verb) ) message(" LOO re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting gp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else if ( (method == 'sampling')&("median" %in% names(object$loo)) & (M == object$loo$M) ){
          if (sample_size == object$loo$sample_size){
            if ( isTRUE(verb) ) message(" LOO re-evaluation not needed.")
            if ( isTRUE(verb) ) message("Exporting gp object without re-evaluation ...", appendLF = FALSE)
            if ( isTRUE(verb) ) Sys.sleep(0.5)
            if ( isTRUE(verb) ) message(" done")
            return(object)
          } else {
            if ( isTRUE(verb) ) message(" LOO re-evaluation needed.")
            if ( isTRUE(verb) ) message("Start re-evaluation: ")
          }
        } else {
          if ( isTRUE(verb) ) message(" LOO re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the LOO ...", appendLF = FALSE)
    x_train <- object$constructor_obj$X
    y_train <- object$constructor_obj$Y
    dat <- list('x_train' = x_train, 'y_train' = y_train)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the LOO ...", appendLF = FALSE)
    res <- object$emulator_obj$loo(method = method, sample_size = sample_size, m = M)
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
    dat[["M"]] <- M
    if (method == "sampling"){
      dat[["sample_size"]] <- sample_size
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
    if ( is.vector(x_test) ) {
      if ( ncol(object$data$X)!=1 ){
        x_test <- matrix(x_test, nrow = 1)
      } else {
        x_test <- as.matrix(x_test)
      }
    }
    if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)
    if ( ncol(x_test) != ncol(object$data$X) ) stop("'x_test' must have the same number of dimensions as the training input.", call. = FALSE)
    #check existing OOS
    if ( isFALSE(force) ){
      if ( "oos" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" OOS results found in the gp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (method == 'mean_var')&("mean" %in% names(object$oos)) & (M == object$oos$M) ){
          if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting gp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else if ( ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (method == 'sampling')&("median" %in% names(object$oos)) & (M == object$oos$M) ) ){
          if ( sample_size == object$oos$sample_size ){
            if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
            if ( isTRUE(verb) ) message("Exporting gp object without re-evaluation ...", appendLF = FALSE)
            if ( isTRUE(verb) ) Sys.sleep(0.5)
            if ( isTRUE(verb) ) message(" done")
            return(object)
          } else {
            if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
            if ( isTRUE(verb) ) message("Start re-evaluation: ")
          }
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
      res <- object$emulator_obj$predict(x_test_unique, method = method, sample_size = sample_size, m = M)
    } else {
      res <- object$emulator_obj$ppredict(x_test_unique, method = method, sample_size = sample_size, m = M, core_num = cores)
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
    dat[["M"]] <- M
    if (method == "sampling"){
      dat[["sample_size"]] <- sample_size
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
validate.dgp <- function(object, x_test = NULL, y_test = NULL, method = NULL, sample_size = 50, verb = TRUE, M = 50, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }
  M <- as.integer(M)

  L = object$constructor_obj$n_layer
  final_node <- object$specs[[paste('layer', L, sep="")]][['node1']]
  if ("type" %in% names(final_node) && final_node$type == "Categorical") {
    is.categorical <- TRUE
    is.Poisson <- FALSE
    is.NegBin <- FALSE
  } else if ("type" %in% names(final_node) && final_node$type == "Poisson") {
    is.categorical <- FALSE
    is.Poisson <- TRUE
    is.NegBin <- FALSE
  } else if ("type" %in% names(final_node) && final_node$type == "NegBin") {
    is.categorical <- FALSE
    is.Poisson <- FALSE
    is.NegBin <- TRUE
  } else {
    is.categorical <- FALSE
    is.Poisson <- FALSE
    is.NegBin <- FALSE
  }

  if ( is.null(method) ){
    if (is.categorical|is.Poisson|is.NegBin) {
      method = 'sampling'
    } else {
      method = 'mean_var'
    }
  } else {
    if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)
    if ( method=='mean_var' && is.categorical){
      stop("'method' can only be 'sampling' for DGP emulators with categorical likelihoods.", call. = FALSE)
    }
  }

  sample_size <- as.integer(sample_size)
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
        if ( (method == 'mean_var')&("mean" %in% names(object$loo)) & (M == object$loo$M) ){
          if ( isTRUE(verb) ) message(" LOO re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting dgp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else if ( (method == 'sampling') & (any(c("median", "label") %in% names(object$loo)))  & (M == object$loo$M) ) {
          if ( sample_size == object$loo$sample_size ) {
            if ( isTRUE(verb) ) message(" LOO re-evaluation not needed.")
            if ( isTRUE(verb) ) message("Exporting dgp object without re-evaluation ...", appendLF = FALSE)
            if ( isTRUE(verb) ) Sys.sleep(0.5)
            if ( isTRUE(verb) ) message(" done")
            return(object)
          } else {
            if ( isTRUE(verb) ) message(" LOO re-evaluation needed.")
            if ( isTRUE(verb) ) message("Start re-evaluation: ")
          }
        } else {
          if ( isTRUE(verb) ) message(" LOO re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the LOO ...", appendLF = FALSE)
    x_train <- object$constructor_obj$X
    y_train <- object$data$Y
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
      res <- object$emulator_obj$loo(X = reticulate::np_array(x_train), method = method, sample_size = sample_size, m = M)
    } else {
      res <- object$emulator_obj$ploo(X = reticulate::np_array(x_train), method = method, sample_size = sample_size, m = M, core_num = cores)
    }
    if ( isTRUE(verb) ) message(" done")
    if ( isTRUE(verb) ) message("Saving results to the slot 'loo' in the dgp object ...", appendLF = FALSE)
    if ( method == 'sampling' ){
      if ( is.categorical ) {
        encoder <- object$constructor_obj$all_layer[[L]][[1]]$class_encoder
        prob_samp <- pkg.env$np$transpose(pkg.env$np$asarray(res), c(2L,1L,0L))
        label_sample <- pkg.env$dgpsi$functions$categorical_sampler_3d(prob_samp)
        original_label <- encoder$inverse_transform(pkg.env$np$ravel(label_sample))
        original_label = pkg.env$np$reshape(original_label, pkg.env$np$shape(label_sample) )
        dat[["label"]] <- original_label
        index_to_label <- as.character(encoder$classes_)
        names(res) <- index_to_label
        dat[["probability"]] <- res
        y_train_encode <- encoder$transform(pkg.env$np$ravel(dat$y_train))
        dat[["log_loss"]] <- pkg.env$dgpsi$functions$logloss(prob_samp, as.integer(y_train_encode))
        dat[["accuracy"]] <- mean(rowMeans(pkg.env$np$equal(original_label, dat$y_train)))
      } else {
        res_np <- pkg.env$np$array(res)
        quant <- pkg.env$np$transpose(pkg.env$np$quantile(res_np, c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
        std <- pkg.env$np$std(res_np, axis=2L)
        dat[["median"]] <- as.matrix(quant[2,,])
        dat[["std"]] <- t(std)
        dat[["lower"]] <- as.matrix(quant[1,,])
        dat[["upper"]] <- as.matrix(quant[3,,])
        dat[["rmse"]] <- sqrt(colMeans((dat$median-dat$y_train)^2))
        dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_train, axis=0L)-pkg.env$np$amin(dat$y_train, axis=0L))
      }
    } else if ( method == 'mean_var' ) {
      dat[["mean"]] <- res[[1]]
      dat[["std"]] <- sqrt(res[[2]])
      dat[["lower"]] <- dat$mean-2*dat$std
      dat[["upper"]] <- dat$mean+2*dat$std
      dat[["rmse"]] <- sqrt(colMeans((dat$mean-dat$y_train)^2))
      dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_train, axis=0L)-pkg.env$np$amin(dat$y_train, axis=0L))
    }
    dat[["M"]] <- M
    if (method == "sampling"){
      dat[["sample_size"]] <- sample_size
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
    if ( is.vector(x_test) ) {
      if ( ncol(object$data$X)!=1 ){
        x_test <- matrix(x_test, nrow = 1)
      } else {
        x_test <- as.matrix(x_test)
      }
    }
    if ( is.vector(y_test) ) {
      if ( ncol(object$data$Y)!=1 ){
        y_test <- matrix(y_test, nrow = 1)
      } else {
        y_test <- as.matrix(y_test)
      }
    }
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)
    if ( ncol(x_test) != ncol(object$data$X) ) stop("'x_test' must have the same number of dimensions as the training input.", call. = FALSE)
    #check existing OOS
    if ( isFALSE(force) ){
      if ( "oos" %in% names(object) ){
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( isTRUE(verb) ) message(" OOS results found in the dgp object.")
        if ( isTRUE(verb) ) message("Checking ...", appendLF = FALSE)
        if ( isTRUE(verb) ) Sys.sleep(0.5)
        if ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (method == 'mean_var')&("mean" %in% names(object$oos)) & (M == object$oos$M) ){
          if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting dgp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else if ( identical(object$oos$x_test, x_test) & identical(object$oos$y_test, y_test) & (method == 'sampling')&(any(c("median", "label") %in% names(object$oos))) & (M == object$oos$M) ){
          if ( sample_size == object$oos$sample_size ){
            if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
            if ( isTRUE(verb) ) message("Exporting dgp object without re-evaluation ...", appendLF = FALSE)
            if ( isTRUE(verb) ) Sys.sleep(0.5)
            if ( isTRUE(verb) ) message(" done")
            return(object)
          } else {
            if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
            if ( isTRUE(verb) ) message("Start re-evaluation: ")
          }
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
      if (is.categorical) {
        res <- object$emulator_obj$classify(x = reticulate::np_array(x_test_unique), mode = 'prob', sample_size = sample_size, m = M)
      } else {
        res <- object$emulator_obj$predict(x = x_test_unique, method = method, sample_size = sample_size, m = M)
      }
    } else {
      if (is.categorical) {
        res <- object$emulator_obj$pclassify(x = reticulate::np_array(x_test_unique), mode = 'prob', sample_size = sample_size, m = M, core_num = cores)
      } else {
        res <- object$emulator_obj$ppredict(x = x_test_unique, method = method, sample_size = sample_size, m = M, core_num = cores)
      }
    }
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Saving results to the slot 'oos' in the dgp object ...", appendLF = FALSE)
    if ( method == 'sampling' ){
      if ( is.categorical ) {
        encoder <- object$constructor_obj$all_layer[[L]][[1]]$class_encoder
        prob_samp <- pkg.env$np$transpose(pkg.env$np$asarray(res), c(2L,1L,0L))
        label_sample <- pkg.env$dgpsi$functions$categorical_sampler_3d(prob_samp)
        original_label <- encoder$inverse_transform(pkg.env$np$ravel(label_sample))
        original_label = pkg.env$np$reshape(original_label, pkg.env$np$shape(label_sample) )
        dat[["label"]] <- original_label[rep,,drop=F]
        index_to_label <- as.character(encoder$classes_)
        for (k in 1:length(res)){
          res[[k]] <- res[[k]][rep,,drop=F]
        }
        names(res) <- index_to_label
        dat[["probability"]] <- res
        y_test_encode <- encoder$transform(pkg.env$np$ravel(dat$y_test))
        dat[["log_loss"]] <- pkg.env$dgpsi$functions$logloss(prob_samp, as.integer(y_test_encode), as.integer(rep-1))
        dat[["accuracy"]] <- mean(rowMeans(pkg.env$np$equal(dat[["label"]], dat$y_test)))
      } else {
        res_np <- pkg.env$np$array(res)
        quant <- pkg.env$np$transpose(pkg.env$np$quantile(res_np, c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
        std <- pkg.env$np$std(res_np, axis=2L)
        dat[["median"]] <- as.matrix(quant[2,,])[rep,,drop=F]
        dat[["std"]] <- t(std)[rep,,drop=F]
        dat[["lower"]] <- as.matrix(quant[1,,])[rep,,drop=F]
        dat[["upper"]] <- as.matrix(quant[3,,])[rep,,drop=F]
        dat[["rmse"]] <- sqrt(colMeans((dat$median-dat$y_test)^2))
        dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_test, axis=0L)-pkg.env$np$amin(dat$y_test, axis=0L))
      }
    } else if ( method == 'mean_var' ) {
      dat[["mean"]] <- res[[1]][rep,,drop=F]
      dat[["std"]] <- sqrt(res[[2]][rep,,drop=F])
      dat[["lower"]] <- dat$mean-2*dat$std
      dat[["upper"]] <- dat$mean+2*dat$std
      dat[["rmse"]] <- sqrt(colMeans((dat$mean-dat$y_test)^2))
      dat[["nrmse"]] <- dat$rmse/(pkg.env$np$amax(dat$y_test, axis=0L)-pkg.env$np$amin(dat$y_test, axis=0L))
    }
    dat[["M"]] <- M
    if (method == "sampling"){
      dat[["sample_size"]] <- sample_size
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
validate.lgp <- function(object, x_test = NULL, y_test = NULL, method = NULL, sample_size = 50, verb = TRUE, M = 50, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"lgp") ) stop("'object' must be an instance of the 'lgp' class.", call. = FALSE)

  if ( "metadata" %in% names(object$specs) ){
    if ( !("emulator_obj" %in% names(object)) ){
      stop("'object' is not activated for predictions. Please set `activate = TRUE` in `lgp()` to activate the emulator.", call. = FALSE)
    }
  }
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }
  M <- as.integer(M)

  if ( is.null(method) ){
    method = 'mean_var'
  } else {
    if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)
  }
  sample_size <- as.integer(sample_size)
  #For OOS
  if (!is.null(x_test) & !is.null(y_test)) {
    #check testing input
    if ( "metadata" %in% names(object$specs) ){
      if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
      x_test <- unname(x_test)
      global_dim <- unique(subset(object$specs$struc, object$specs$struc[["From_Emulator"]] == "Global")$From_Output)
      if ( is.vector(x_test) ) {
        if ( global_dim!=1 ){
          x_test <- matrix(x_test, nrow = 1)
          if ( ncol(x_test)<global_dim ) stop(sprintf("'x_test' has missing dimensions. Expected %d, found %d in 'x'.",
                                                 global_dim, ncol(x_test)), call. = FALSE)
        } else {
          x_test <- as.matrix(x_test)
        }
      }
      nrow_x <- nrow(x_test)
      num_layers <- length(object$constructor)
      x_list <- vector("list", num_layers)
      x_list[[1]] <- x_test
      for (l in 2:num_layers) {
        layer_metadata <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == l)
        layer_metadata <- layer_metadata[order(layer_metadata$Pos_in_Layer), ]
        layer_matrices <- vector("list", nrow(layer_metadata))
        for (k in seq_len(nrow(layer_metadata))) {
          emulator_id <- layer_metadata$Emulator[k]
          input_connections <- subset(object$specs$struc, object$specs$struc[["To_Emulator"]] == emulator_id)
          global_outputs <- input_connections$From_Output[input_connections$From_Emulator == "Global"]
          if ( length(global_outputs)>0 ) {
            layer_matrices[[k]] <- x_test[,global_outputs,drop=F]
          }
        }
        x_list[[l]] <- layer_matrices
      }
      x_test <- x_list
    } else {
      lifecycle::deprecate_warn(
        when = "2.5.0",
        what = I("The `object` created by `lgp()` without specifying `struc` as a data frame"),
        details = c(
          i = "Support for `object` structures created without `struc` specified as a data frame will be removed in the next release.",
          i = "To ensure compatibility with future versions, please recreate the `object` by calling the updated `lgp()` with `struc` provided as a data frame."
        )
      )
      # To be deprecated
      if ( !is.list(x_test) ) {
        if ( !is.matrix(x_test)&!is.vector(x_test) ) {
          stop("'x_test' must be a vector or a matrix.", call. = FALSE)
        } else {
          x_test <- unname(x_test)
          if ( is.vector(x_test) ) {
            is.1d <- TRUE
            for (item in object$constructor_obj[[1]]){
              if (item$type == 'gp'){
                if (length(item$structure$input_dim) != 1){
                  is.1d <- FALSE
                  break
                }
              } else {
                for (ker in item$structure[[1]]){
                  if (length(ker$input_dim) != 1){
                    is.1d <- FALSE
                    break
                  }
                }
                if (isFALSE(is.1d)) break
              }
            }
            if ( is.1d ){
              x_test <- as.matrix(x_test)
            } else {
              x_test <- matrix(x_test, nrow = 1)
            }
          }
          nrow_x <- nrow(x_test)
        }
      } else {
        for ( l in 1:length(x_test) ){
          if ( l==1 ){
            if ( !is.matrix(x_test[[l]])&!is.vector(x_test[[l]]) ) {
              stop("The first element of 'x_test' must be a vector or a matrix.", call. = FALSE)
            } else {
              x_test[[l]] <- unname(x_test[[l]])
              if ( is.vector(x_test[[l]]) ) {
                is.1d <- TRUE
                for (item in object$constructor_obj[[l]]){
                  if (item$type == 'gp'){
                    if (length(item$structure$input_dim) != 1){
                      is.1d <- FALSE
                      break
                    }
                  } else {
                    for (ker in item$structure[[1]]){
                      if (length(ker$input_dim) != 1){
                        is.1d <- FALSE
                        break
                      }
                    }
                    if (isFALSE(is.1d)) break
                  }
                }
                if ( is.1d ){
                  x_test[[l]] <- as.matrix(x_test[[l]])
                } else {
                  x_test[[l]] <- matrix(x_test[[l]], nrow = 1)
                }
              }
              nrow_x <- nrow(x_test[[l]])
            }
          } else {
            for ( k in 1:length(x_test[[l]]) ){
              if ( !is.matrix(x_test[[l]][[k]])&!is.null(x_test[[l]][[k]])&!is.vector(x_test[[l]][[k]]) ) stop(sprintf("The element %i in the sublist %i of 'x_test' must be a vector, a matrix, or 'NULL'.", k, l), call. = FALSE)
              if ( is.matrix(x_test[[l]][[k]])|is.vector(x_test[[l]][[k]]) ){
                x_test[[l]][[k]] <- unname(x_test[[l]][[k]])
                if (is.vector(x_test[[l]][[k]])) {
                  is.1d <- TRUE
                  item_cont <- object$constructor_obj[[l]][[k]]
                  if (item_cont$type == 'gp'){
                    if (length(item_cont$structure$input_dim) != 1){
                      is.1d <- FALSE
                    }
                  } else {
                    for (ker in item_cont$structure[[1]]){
                      if (length(ker$input_dim) != 1){
                        is.1d <- FALSE
                        break
                      }
                    }
                  }
                  if ( is.1d ){
                    x_test[[l]][[k]] <- as.matrix(x_test[[l]][[k]])
                  } else {
                    x_test[[l]][[k]] <- matrix(x_test[[l]][[k]], nrow = 1)
                  }
                }
                if ( nrow(x_test[[l]][[k]])!=nrow_x ) {
                  stop(sprintf("The element %i in the sublist %i of 'x_test' has an inconsistent number of data points with the first element of 'x_test'.", k, l), call. = FALSE)
                }
              }
            }
          }
        }
      }
    }

    #check testing output
    if ( "metadata" %in% names(object$specs) ){
      if ( !is.list(y_test) ) {
        if ( !is.matrix(y_test)&!is.vector(y_test) ) {
          stop("'y_test' must be a vector or a matrix.", call. = FALSE)
        } else {
          final_layer <- max(object$specs$metadata$Layer)
          # Filter metadata to get emulators in the final layer
          final_layer_emulators <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == final_layer)
          emu_num <- nrow(final_layer_emulators)
          if ( emu_num!=1 ) stop(sprintf("The linked system contains %d emulators in its final layer. 'y_test' must be a list of %d vectors or matrices.",
                                         emu_num, emu_num), call. = FALSE)
          y_test <- unname(y_test)
          if ( is.vector(y_test) ) {
            output_dim <- final_layer_emulators$Total_Output_Dims
            if (output_dim==1){
              y_test <- as.matrix(y_test)
            } else {
              y_test <- matrix(y_test, nrow = 1)
            }
          }
        }
        nrow_y <- nrow(y_test)
        if ( nrow_y!=nrow_x ) stop("The number of data points are inconsistent between 'x_test' and 'y_test'.", call. = FALSE)
      } else {
        final_layer <- max(object$specs$metadata$Layer)
        # Filter metadata to get emulators in the final layer
        final_layer_emulators <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == final_layer)
        emu_num <- nrow(final_layer_emulators)
        if ( emu_num!=length(y_test) ) stop(sprintf("The linked system's final layer contains %d emulators. 'y_test' should contain %d vectors or matrices, but found %d.",
                                                    emu_num, emu_num, length(y_test)), call. = FALSE)
        for ( l in 1:length(y_test) ){
          if ( !is.matrix(y_test[[l]])&!is.vector(y_test[[l]]) ) {
            stop(sprintf("The element %i of 'y_test' must be a vector or a matrix.", l), call. = FALSE)
          } else {
            y_test[[l]] <- unname(y_test[[l]])
            if ( is.vector(y_test[[l]]) ) {
              emu_cont <- subset(final_layer_emulators, final_layer_emulators[["Pos_in_Layer"]] == l)
              emu_cont_output_dim <- emu_cont$Total_Output_Dims
              if (emu_cont_output_dim==1){
                y_test[[l]] <- as.matrix(y_test[[l]])
              } else {
                y_test[[l]] <- matrix(y_test[[l]], nrow = 1)
              }
            }
            nrow_y <- nrow(y_test[[l]])
            if ( nrow_y!=nrow_x ) stop(sprintf("The number of data points is inconsistent between 'x_test' and the element %i of 'y_test'.", l), call. = FALSE)
          }
        }
      }
    } else {
      if ( !is.list(y_test) ) {
        if ( !is.matrix(y_test)&!is.vector(y_test) ) {
          stop("'y_test' must be a vector or a matrix.", call. = FALSE)
        } else {
          total_layer <- length(object$constructor_obj)
          emu_num <- length(object$constructor_obj[[total_layer]])
          if ( emu_num!=1 ) stop(sprintf("The linked system contains %d emulators in its final layer. 'y_test' must be a list of %d vectors or matrices.",
                                          emu_num, emu_num), call. = FALSE)
          y_test <- unname(y_test)
          if ( is.vector(y_test) ) {
            emu_cont <- object$constructor_obj[[total_layer]][[emu_num]]
            if (emu_cont$type == 'gp'){
              y_test <- as.matrix(y_test)
            } else {
              dgp_layer <- length(emu_cont$structure)
              if (length(emu_cont$structure[[dgp_layer]])==1){
                y_test <- as.matrix(y_test)
              } else {
                y_test <- matrix(y_test, nrow = 1)
              }
            }
          }
          nrow_y <- nrow(y_test)
          if ( nrow_y!=nrow_x ) stop("The number of data points is inconsistent between 'x_test' and 'y_test'.", call. = FALSE)
        }
      } else {
        total_layer <- length(object$constructor_obj)
        emu_num <- length(object$constructor_obj[[total_layer]])
        if ( emu_num!=length(y_test) ) stop(sprintf("The linked system's final layer contains %d emulators. 'y_test' should contain %d vectors or matrices, but found %d.",
                                                    emu_num, emu_num, length(y_test)), call. = FALSE)
        for ( l in 1:length(y_test) ){
          if ( !is.matrix(y_test[[l]])&!is.vector(y_test[[l]]) ) {
            stop(sprintf("The element %i of 'y_test' must be a vector or a matrix.", l), call. = FALSE)
          } else {
            y_test[[l]] <- unname(y_test[[l]])
            if ( is.vector(y_test[[l]]) ) {
              emu_cont <- object$constructor_obj[[total_layer]][[l]]
              if (emu_cont$type == 'gp'){
                y_test[[l]] <- as.matrix(y_test[[l]])
              } else {
                dgp_layer <- length(emu_cont$structure)
                if (length(emu_cont$structure[[dgp_layer]])==1){
                  y_test[[l]] <- as.matrix(y_test[[l]])
                } else {
                  y_test[[l]] <- matrix(y_test[[l]], nrow = 1)
                }
              }
            }
            nrow_y <- nrow(y_test[[l]])
            if ( nrow_y!=nrow_x ) stop(sprintf("The number of data points is inconsistent between 'x_test' and the element %i of 'y_test'.", l), call. = FALSE)
          }
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
        if ( identical(object$oos$x_test, if ( "metadata" %in% names(object$specs) ){x_test[[1]]} else {x_test} ) & identical(object$oos$y_test, y_test) & (method == 'mean_var')&("mean" %in% names(object$oos)) & (M == object$oos$M) ){
          if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
          if ( isTRUE(verb) ) message("Exporting lgp object without re-evaluation ...", appendLF = FALSE)
          if ( isTRUE(verb) ) Sys.sleep(0.5)
          if ( isTRUE(verb) ) message(" done")
          return(object)
        } else if ( identical(object$oos$x_test, if ( "metadata" %in% names(object$specs) ){x_test[[1]]} else {x_test}) & identical(object$oos$y_test, y_test) & (method == 'sampling')&("median" %in% names(object$oos)) & (M == object$oos$M) ){
          if ( sample_size == object$oos$sample_size ) {
            if ( isTRUE(verb) ) message(" OOS re-evaluation not needed.")
            if ( isTRUE(verb) ) message("Exporting lgp object without re-evaluation ...", appendLF = FALSE)
            if ( isTRUE(verb) ) Sys.sleep(0.5)
            if ( isTRUE(verb) ) message(" done")
            return(object)
          } else {
            if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
            if ( isTRUE(verb) ) message("Start re-evaluation: ")
          }
        } else {
          if ( isTRUE(verb) ) message(" OOS re-evaluation needed.")
          if ( isTRUE(verb) ) message("Start re-evaluation: ")
        }
      }
    }

    if ( isTRUE(verb) ) message("Initializing the OOS ...", appendLF = FALSE)
    if ( "metadata" %in% names(object$specs) ){
      dat <- list('x_test' = x_test[[1]],'y_test' = y_test)
    } else {
      dat <- list('x_test' = x_test,'y_test' = y_test)
    }
    if ( !is.list(y_test) ) y_test <- list(y_test)
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Calculating the OOS ...", appendLF = FALSE)
    if ( identical(cores,as.integer(1)) ){
      res <- object$emulator_obj$predict(x = x_test, method = method, sample_size = sample_size, m = M)
    } else {
      res <- object$emulator_obj$ppredict(x = x_test, method = method, sample_size = sample_size, m = M, core_num = cores)
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
    dat[["M"]] <- M
    if (method == "sampling"){
      dat[["sample_size"]] <- sample_size
    }
    object$oos <- dat
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    return(object)
    #For other cases
  } else {
    stop("Both 'x_test' and 'y_test' must be provided for validation of a linked (D)GP emulator.", call. = FALSE)
  }
}
