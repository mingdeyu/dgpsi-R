#' @title Prediction from GP, DGP, or linked (D)GP emulators
#'
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function implements prediction from GP, DGP, or linked (D)GP emulators.
#'
#' @param object an instance of the `gp`, `dgp`, or `lgp` class.
#' @param x the testing input data:
#' * if `object` is an instance of the `gp` or `dgp` class, `x` is a matrix where each row is an input testing data point and each column is an input dimension.
#' * `r lifecycle::badge("deprecated")` if `object` is an instance of the `lgp` class created by [lgp()] without specifying argument `struc` in data frame form, `x` can be either a matrix or a list:
#'    - if `x` is a matrix, its rows are treated as instances of the `Global` inputs. In this case, it is assumed that the only global input to the system is the input to the
#'      emulators in the first layer and there is no global input to emulators in other layers.
#'    - if `x` is a list, it should have *L* (the number of layers in an emulator system) elements. The first element
#'      is a matrix that represents the global testing input data that feed into the emulators in the first layer of the system. The
#'      remaining *L-1* elements are *L-1* sub-lists, each of which contains a number (the same number of emulators in
#'      the corresponding layer) of matrices (rows being testing input data points and columns being input dimensions) that represent the
#'      global testing input data to the emulators in the corresponding layer. The matrices must be placed in the sub-lists based on how
#'      their corresponding emulators are placed in `struc` argument of [lgp()]. If there is no global input data to a certain emulator,
#'      set `NULL` in the corresponding sub-list of `x`.
#'
#'    **This option for linked (D)GP emulators is deprecated and will be removed in the next release.**
#' * `r new_badge("new")` If `object` is an instance of the `lgp` class created by [lgp()] with argument `struc` in data frame form,
#'   `x` must be a matrix representing the global input, where each row corresponds to a test data point and each column represents a global input dimension.
#'   The column indices in `x` must align with the indices specified in the `From_Output` column of the `struc` data frame (used in [lgp()]),
#'   corresponding to rows where the `From_Emulator` column is `"Global"`.
#' @param method `r new_badge("updated")` the prediction approach to use: either the mean-variance approach (`"mean_var"`) or the sampling approach (`"sampling"`).
#'     The mean-variance approach returns the means and variances of the predictive distributions, while the sampling approach generates samples from predictive distributions
#'     using the derived means and variances. For DGP emulators with a categorical likelihood (`likelihood = "Categorical"` in [dgp()]), `method` is only applicable
#'     when `full_layer = TRUE`. In this case, the sampling approach generates samples from the GP nodes in all hidden layers using the derived means and variances,
#'     and subsequently propagates these samples through the categorical likelihood. By default, the method is set to `"sampling"` for DGP emulators with Poisson, Negative Binomial, and
#'     Categorical likelihoods, and to `"mean_var"` otherwise.
#' @param mode `r new_badge("new")` whether to predict the classes (`"label"`) or probabilities (`"proba"`) of different classes when `object` is a DGP emulator with a categorical likelihood.
#'      Defaults to `"label"`.
#' @param full_layer a bool indicating whether to output the predictions of all layers. Defaults to `FALSE`. Only used when `object` is a DGP or a linked (D)GP emulator.
#' @param sample_size the number of samples to draw for each given imputation if `method = "sampling"`. Defaults to `50`.
#' @param M `r new_badge("new")` the size of the conditioning set for the Vecchia approximation in the emulator prediction. Defaults to `50`. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation.
#' @param cores the number of processes to be used for prediction. If set to `NULL`, the number of processes is set to `max physical cores available %/% 2`. Defaults to `1`.
#' @param chunks the number of chunks that the testing input matrix `x` will be divided into for multi-cores to work on.
#'     Only used when `cores` is not `1`. If not specified (i.e., `chunks = NULL`), the number of chunks is set to the value of `cores`.
#'     Defaults to `NULL`.
#' @param ... N/A.
#' @return
#' * If `object` is an instance of the `gp` class:
#'   1. if `method = "mean_var"`: an updated `object` is returned with an additional slot called `results` that contains two matrices named `mean`
#'      for the predictive means and `var` for the predictive variances. Each matrix has only one column with its rows
#'      corresponding to testing positions (i.e., rows of `x`).
#'   2. if `method = "sampling"`: an updated `object` is returned with an additional slot called `results` that contains a matrix whose rows correspond
#'      to testing positions and columns correspond to `sample_size` number of samples drawn from the predictive distribution of GP.
#' * If `object` is an instance of the `dgp` class:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains two
#'      matrices named `mean` for the predictive means and `var` for the predictive variances respectively. Each matrix has its rows corresponding to testing
#'      positions and columns corresponding to DGP global output dimensions (i.e., the number of GP/likelihood nodes in the final layer).
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains two
#'      sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L* (i.e., the number of layers)
#'      matrices named `layer1, layer2,..., layerL`. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      output dimensions (i.e., the number of GP/likelihood nodes from the associated layer).
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains *D* (i.e., the number
#'      of GP/likelihood nodes in the final layer) matrices named `output1, output2,..., outputD`. Each matrix has its rows corresponding to testing positions and
#'      columns corresponding to samples of size: `B * sample_size`, where `B` is the number of imputations specified in [dgp()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains *L* (i.e., the number
#'      of layers) sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents samples drawn from the GP/likelihood nodes in the corresponding layer,
#'      and contains *D* (i.e., the number of GP/likelihood nodes in the corresponding layer) matrices named `output1, output2,..., outputD`. Each matrix gives samples
#'      of the output from one of *D* GP/likelihood nodes, and has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `B * sample_size`, where `B` is the number of imputations specified in [dgp()].
#' * `r new_badge("new")` If `object` is an instance of the `dgp` class with a categorical likelihood:
#'   1. if `full_layer = FALSE` and `mode = "label"`: an updated `object` is returned with an additional slot called `results` that contains one matrix named `label`.
#'      The matrix has rows corresponding to testing positions and columns corresponding to sample labels of size: `B * sample_size`, where `B` is the number
#'      of imputations specified in [dgp()].
#'   2. if `full_layer = FALSE` and `mode = "proba"`, an updated `object` is returned with an additional slot called `results`. This slot contains *D* matrices (where
#'      *D* is the number of classes in the training output), where each matrix gives probability samples for the corresponding class with its rows corresponding to testing
#'      positions and columns containing probabilities. The number of columns of each matrix is `B * sample_size`, where `B` is the number of imputations
#'      specified in the [dgp()] function.
#'   3. if `method = "mean_var"` and `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains *L* (i.e., the number
#'      of layers) sub-lists named `layer1, layer2,..., layerL`. Each of first `L-1` sub-lists contains two matrices named `mean` for the predictive means and `var`
#'      for the predictive variances of the GP nodes in the associated layer. Rows of each matrix correspond to testing positions.
#'      - when `mode = "label"`, the sub-list `LayerL` contains one matrix named `label`. The matrix has its rows corresponding to testing positions and columns
#'        corresponding to label samples of size: `B * sample_size`. `B` is the number of imputations specified in [dgp()].
#'      - when `mode = "proba"`, the sub-list `LayerL` contains *D* matrices (where *D* is the number of classes in the training output), where each matrix gives probability
#'        samples for the corresponding class with its rows corresponding to testing positions and columns containing probabilities. The number of columns of each matrix
#'        is `B * sample_size`. `B` is the number of imputations specified in [dgp()].
#'   4. if `method = "sampling"` and `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains *L* (i.e., the number
#'      of layers) sub-lists named `layer1, layer2,..., layerL`. Each of first `L-1` sub-lists represents samples drawn from the GP nodes in the
#'      corresponding layer, and contains *D* (i.e., the number of GP nodes in the corresponding layer) matrices named `output1, output2,..., outputD`. Each matrix
#'      gives samples of the output from one of *D* GP nodes, and has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `B * sample_size`.
#'      - when `mode = "label"`, the sub-list `LayerL` contains one matrix named `label`. The matrix has its rows corresponding to testing positions and columns
#'        corresponding to label samples of size: `B * sample_size`.
#'      - when `mode = "proba"`, the sub-list `LayerL` contains *D* matrices (where *D* is the number of classes in the training output), where each matrix gives probability
#'        samples for the corresponding class with its rows corresponding to testing positions and columns containing probabilities. The number of columns of each matrix
#'        is `B * sample_size`.
#'
#'      `B` is the number of imputations specified in [dgp()].
#' * `r new_badge("updated")` If `object` is an instance of the `lgp` class:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that
#'      contains two sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list
#'      contains *K* (same number of emulators in the final layer of the system) matrices named using the `ID`s of the corresponding emulators in the final layer.
#'      Each matrix has rows corresponding to global testing positions and columns corresponding to output dimensions of the associated emulator
#'      in the final layer.
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains
#'      two sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L*
#'      (i.e., the number of layers in the emulated system) components named `layer1, layer2,..., layerL`. Each component represents a layer
#'      and contains *K* (same number of emulators in the corresponding layer of the system) matrices named using the `ID`s of the corresponding emulators in that layer.
#'      Each matrix has its rows corresponding to global testing positions and columns corresponding to output dimensions of the associated
#'      GP/DGP emulator in the corresponding layer.
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains
#'      *K* (same number of emulators in the final layer of the system) sub-lists named using the `ID`s of the corresponding emulators in the final layer. Each
#'      sub-list contains *D* matrices, named `output1, output2,..., outputD`, that correspond to the output
#'      dimensions of the GP/DGP emulator. Each matrix has rows corresponding to testing positions and columns corresponding to samples
#'      of size: `B * sample_size`, where `B` is the number of imputations specified in [lgp()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains
#'      *L* (i.e., the number of layers of the emulated system) sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents a layer
#'      and contains *K* (same number of emulators in the corresponding layer of the system) components named using the `ID`s of the corresponding emulators in that layer.
#'      Each component contains *D* matrices, named `output1, output2,..., outputD`, that correspond to
#'      the output dimensions of the GP/DGP emulator. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      samples of size: `B * sample_size`, where `B` is the number of imputations specified in [lgp()].
#'
#'   If `object` is an instance of the `lgp` class created by [lgp()] without specifying the `struc` argument in data frame form, the `ID`s, that are used as names of sub-lists or
#'   matrices within `results`, will be replaced by `emulator1`, `emulator2`, and so on.
#'
#' The `results` slot will also include:
#' * `r new_badge("new")` the value of `M`, which represents the size of the conditioning set for the Vecchia approximation, if used, in the emulator prediction.
#' * the value of `sample_size` if `method = "sampling"`.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name predict
NULL


#' @rdname predict
#' @method predict dgp
#' @export
predict.dgp <- function(object, x, method = NULL, mode = 'label', full_layer = FALSE, sample_size = 50, M = 50, cores = 1, chunks = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) {
    if ( ncol(object$data$X)!=1 ){
      x <- matrix(x, nrow = 1)
    } else {
      x <- as.matrix(x)
    }
  }
  if ( ncol(x) != ncol(object$data$X) ) stop("'x' must have the same number of dimensions as the training input.", call. = FALSE)

  L = object$constructor_obj$n_layer
  final_node <- object$specs[[paste('layer', L, sep="")]][['node1']]
  if ("type" %in% names(final_node) && final_node$type == "Categorical") {
    is.categorical <- TRUE
    is.Poisson <- FALSE
    is.NegBin <- FALSE
    encoder <- object$constructor_obj$all_layer[[L]][[1]]$class_encoder
    index_to_label <- as.character(encoder$classes_)
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
    if ( method=='mean_var' && is.categorical && !full_layer){
      method = 'sampling'
    }
  }

  sample_size <- as.integer(sample_size)
  M <- as.integer(M)
  if( !is.null(chunks) ) {
    chunks <- as.integer(chunks)
    if ( chunks < 1 ) stop("'chunks' must be >= 1.", call. = FALSE)
  }
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  rep_x <- pkg.env$np$unique(x, return_inverse=TRUE, axis=0L)
  x_unique <- rep_x[[1]]
  rep <- rep_x[[2]] + 1

  if ( identical(cores,as.integer(1)) ){
    if (is.categorical) {
      res <- object$emulator_obj$classify(reticulate::np_array(x_unique), mode, method, full_layer, sample_size, M)
    } else {
      res <- object$emulator_obj$predict(x_unique, method, full_layer, sample_size, M)
    }
  } else {
    if (is.categorical) {
      res <- object$emulator_obj$pclassify(reticulate::np_array(x_unique), mode, method, full_layer, sample_size, M, chunks, cores)
    } else {
      res <- object$emulator_obj$ppredict(x_unique, method, full_layer, sample_size, M, chunks, cores)
    }
  }

  if (method == 'mean_var'){
    if (full_layer) {
      if (is.categorical){
        named_res <- vector('list', L)
        for (l in 1:(L-1)) {
          names(named_res)[l] <- paste('layer', l, sep="")
          named_res[[l]][[1]] <- res[[1]][[l]][rep,,drop=F]
          named_res[[l]][[2]] <- res[[2]][[l]][rep,,drop=F]
          names(named_res[[l]])[1] <- 'mean'
          names(named_res[[l]])[2] <- 'var'
        }
        names(named_res)[L] <- paste('layer', L, sep="")
        for (k in 1:length(res[[3]])) {
          named_res[[L]][[k]] <- res[[3]][[k]][rep,,drop=F]
          if ( mode == 'label' ){
            names(named_res[[L]])[k] <- 'label'
          } else {
            names(named_res[[L]])[k] <- index_to_label[k]
          }
        }
      } else {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (l in 1:length(named_res$mean)) {
          named_res$mean[[l]] <- named_res$mean[[l]][rep,,drop=F]
          named_res$var[[l]] <- named_res$var[[l]][rep,,drop=F]
          names(named_res$mean)[l] <- paste('layer', l, sep="")
          names(named_res$var)[l] <- paste('layer', l, sep="")
        }
      }
    } else {
      named_res <- list("mean" = res[[1]][rep,,drop=F], "var" = res[[2]][rep,,drop=F])
    }
  } else if (method == 'sampling') {
    if (full_layer) {
      named_res <- res
      for (l in 1:length(named_res)) {
        names(named_res)[l] <- paste('layer', l, sep="")
        for (k in 1:length(named_res[[l]])) {
          named_res[[l]][[k]] <- named_res[[l]][[k]][rep,,drop=F]
          if ( l == length(named_res) && is.categorical ) {
            if ( mode == 'label' ){
              names(named_res[[l]])[k] <- 'label'
            } else {
              names(named_res[[l]])[k] <- index_to_label[k]
            }
          } else {
            names(named_res[[l]])[k] <- paste('output', k, sep="")
          }
        }
      }
    } else {
      named_res <- res
      for (l in 1:length(named_res)) {
        named_res[[l]] <- named_res[[l]][rep,,drop=F]
        if ( is.categorical ) {
          if ( mode == 'label' ){
            names(named_res)[l] <- 'label'
          } else {
            names(named_res)[l] <- index_to_label[l]
          }
        } else {
          names(named_res)[l] <- paste('output', l, sep="")
        }
      }
    }
  }

  object$results <- named_res
  object$results[["M"]] <- M
  if (method == "sampling" | is.categorical){
    object$results[["sample_size"]] <- sample_size
  }
  return(object)
}


#' @rdname predict
#' @method predict lgp
#' @export
predict.lgp <- function(object, x, method = NULL, full_layer = FALSE, sample_size = 50, M = 50, cores = 1, chunks = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"lgp") ) stop("'object' must be an instance of the 'lgp' class.", call. = FALSE)
  if ( "metadata" %in% names(object$specs) ){
    if ( !("emulator_obj" %in% names(object)) ){
      stop("'object' is not activated for predictions. Please set `activate = TRUE` in `lgp()` to activate the emulator.", call. = FALSE)
    }
    if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
    x <- unname(x)
    global_dim <- unique(subset(object$specs$struc, object$specs$struc[["From_Emulator"]] == "Global")$From_Output)
    if ( is.vector(x) ) {
      if ( global_dim!=1 ){
        x <- matrix(x, nrow = 1)
        if ( ncol(x)<global_dim ) stop(sprintf("'x' has missing dimensions. Expected %d, found %d in 'x'.",
                                                global_dim, ncol(x)), call. = FALSE)
      } else {
        x <- as.matrix(x)
      }
    }
    num_layers <- length(object$constructor)
    x_list <- vector("list", num_layers)
    x_list[[1]] <- x
    for (l in 2:num_layers) {
      layer_metadata <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == l)
      layer_metadata <- layer_metadata[order(layer_metadata$Pos_in_Layer), ]
      layer_matrices <- vector("list", nrow(layer_metadata))
      for (k in seq_len(nrow(layer_metadata))) {
        emulator_id <- layer_metadata$Emulator[k]
        input_connections <- subset(object$specs$struc, object$specs$struc[["To_Emulator"]] == emulator_id)
        global_outputs <- input_connections$From_Output[input_connections$From_Emulator == "Global"]
        if ( length(global_outputs)>0 ) {
          layer_matrices[[k]] <- x[,global_outputs,drop=F]
        }
      }
      x_list[[l]] <- layer_matrices
    }
    x <- x_list
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
    if ( !is.list(x) ) {
      if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
      x <- unname(x)
      if ( is.vector(x) ) {
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
          x <- as.matrix(x)
        } else {
          x <- matrix(x, nrow = 1)
        }
      }
    } else {
      for ( l in 1:length(x) ){
        if ( l==1 ){
          if ( !is.matrix(x[[l]])&!is.vector(x[[l]]) ) {
            stop("The first element of 'x' must be a vector or a matrix.", call. = FALSE)
          } else {
            x[[l]] <- unname(x[[l]])
            if ( is.vector(x[[l]]) ) {
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
                x[[l]] <- as.matrix(x[[l]])
              } else {
                x[[l]] <- matrix(x[[l]], nrow = 1)
              }
            }
            nrow_x <- nrow(x[[l]])
          }
        } else {
          for ( k in 1:length(x[[l]]) ){
            if ( !is.matrix(x[[l]][[k]])&!is.null(x[[l]][[k]])&!is.vector(x[[l]][[k]]) ) stop(sprintf("The element %i in the sublist %i of 'x' must be a vector, a matrix, or 'NULL'.", k, l), call. = FALSE)
            if ( is.matrix(x[[l]][[k]])|is.vector(x[[l]][[k]]) ){
              x[[l]][[k]] <- unname(x[[l]][[k]])
              if (is.vector(x[[l]][[k]])) {
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
                  x[[l]][[k]] <- as.matrix(x[[l]][[k]])
                } else {
                  x[[l]][[k]] <- matrix(x[[l]][[k]], nrow = 1)
                }
              }
              if ( nrow(x[[l]][[k]])!=nrow_x ) {
                stop(sprintf("The element %i in the sublist %i of 'x' has inconsistent number of data points with the first element of 'x'.", k, l), call. = FALSE)
              }
            }
          }
        }
      }
    }
  }


  if ( is.null(method) ){
    method = 'mean_var'
  } else {
    if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)
  }

  sample_size <- as.integer(sample_size)
  M <- as.integer(M)
  if( !is.null(chunks) ) {
    chunks <- as.integer(chunks)
    if ( chunks < 1 ) stop("'chunks' must be >= 1.", call. = FALSE)
  }
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  if ( identical(cores,as.integer(1)) ){
    res <- object$emulator_obj$predict(x, method, full_layer, sample_size, M)
  } else {
    res <- object$emulator_obj$ppredict(x, method, full_layer, sample_size, M, chunks, cores)
  }

  if (method == 'mean_var'){
    if (full_layer) {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      for (l in 1:length(named_res$mean)) {
        names(named_res$mean)[l] <- paste('layer', l, sep="")
        names(named_res$var)[l] <- paste('layer', l, sep="")
        for (k in 1:length(named_res$mean[[l]])) {
          if ( "metadata" %in% names(object$specs) ){
            emulator_id <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == l & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
            names(named_res$mean[[l]])[k] <- emulator_id
            names(named_res$var[[l]])[k] <- emulator_id
          } else {
            names(named_res$mean[[l]])[k] <- paste('emulator', k, sep="")
            names(named_res$var[[l]])[k] <- paste('emulator', k, sep="")
          }
        }
      }
    } else {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      if ( "metadata" %in% names(object$specs) ) final_layer <- max(object$specs$metadata$Layer)
      for (k in 1:length(named_res$mean)) {
        if ( "metadata" %in% names(object$specs) ){
          emulator_id <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == final_layer & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
          names(named_res$mean)[k] <- emulator_id
          names(named_res$var)[k] <- emulator_id
        } else {
          names(named_res$mean)[k] <- paste('emulator', k, sep="")
          names(named_res$var)[k] <- paste('emulator', k, sep="")
        }
      }
    }
  } else if (method == 'sampling') {
    if (full_layer) {
      named_res <- list()
      for (l in 1:length(res)) {
        name1 <- paste('layer', l, sep="")
        for (k in 1:length(res[[l]])){
          if ( "metadata" %in% names(object$specs) ){
            name2 <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == l & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
          } else {
            name2 <- paste('emulator', k, sep="")
          }
          for (n in 1:dim(res[[l]][[k]])[1]) {
            name3 <- paste('output', n, sep="")
            named_res[[name1]][[name2]][[name3]] <- res[[l]][[k]][n,,]
          }
        }
      }
    } else {
      named_res <- list()
      if ( "metadata" %in% names(object$specs) ) final_layer <- max(object$specs$metadata$Layer)
      for (k in 1:length(res)) {
        if ( "metadata" %in% names(object$specs) ) {
          name1 <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == final_layer & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
        } else {
          name1 <- paste('emulator', k, sep="")
        }
        for (n in 1:dim(res[[k]])[1]) {
          name2 <- paste('output', n, sep="")
          named_res[[name1]][[name2]] <- res[[k]][n,,]
        }
      }
    }
  }

  object$results <- named_res
  object$results[["M"]] <- M
  if (method == "sampling"){
    object$results[["sample_size"]] <- sample_size
  }
  return(object)
}


#' @rdname predict
#' @method predict gp
#' @export
predict.gp <- function(object, x, method = NULL, sample_size = 50, M = 50, cores = 1, chunks = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) {
    if ( ncol(object$data$X)!=1 ){
      x <- matrix(x, nrow = 1)
    } else {
      x <- as.matrix(x)
    }
  }
  if ( ncol(x) != ncol(object$data$X) ) stop("'x' must have the same number of dimensions as the training input.", call. = FALSE)

  if ( is.null(method) ){
    method = 'mean_var'
  } else {
    if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)
  }

  sample_size <- as.integer(sample_size)
  M <- as.integer(M)
  if( !is.null(chunks) ) {
    chunks <- as.integer(chunks)
    if ( chunks < 1 ) stop("'chunks' must be >= 1.", call. = FALSE)
  }
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  rep_x <- pkg.env$np$unique(x, return_inverse=TRUE, axis=0L)
  x_unique <- rep_x[[1]]
  rep <- rep_x[[2]] + 1

  if ( identical(cores,as.integer(1)) ){
    res <- object$emulator_obj$predict(x_unique, method, sample_size, M)
  } else {
    res <- object$emulator_obj$ppredict(x_unique, method, sample_size, M, chunks, cores)
  }

  if (method == 'mean_var'){
    object$results <- list("mean" = res[[1]][rep,,drop=FALSE], "var" = res[[2]][rep,,drop=FALSE])
  } else if (method == 'sampling'){
    object$results <- res[rep,,drop=FALSE]
  }
  object$results[["M"]] <- M
  if (method == "sampling"){
    object$results[["sample_size"]] <- sample_size
  }
  return(object)
}

