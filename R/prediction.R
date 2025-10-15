#' @title Prediction from GP, DGP, or linked (D)GP emulators
#'
#' @description This function implements prediction from GP, DGP, or linked (D)GP emulators.
#'
#' @param object an instance of the `gp`, `dgp`, or `lgp` class.
#' @param x the testing input data:
#' * if `object` is an instance of the `gp` or `dgp` class, `x` is a matrix where each row is an input testing data point and each column is an input dimension.
#' * if `object` is an instance of the `lgp` class, `x` must be a matrix representing the global input, where each row corresponds to a test data point and each column represents a global input dimension.
#'   The column indices in `x` must align with the indices specified in the `From_Output` column of the `struc` data frame (used in [lgp()]),
#'   corresponding to rows where the `From_Emulator` column is `"Global"`.
#' @param method `r new_badge("updated")` the prediction approach to use: either the mean-variance approach (`"mean_var"`) or the sampling approach (`"sampling"`).
#'     The mean-variance approach returns the means and variances of the predictive distributions, while the sampling approach generates samples from predictive distributions
#'     using the derived means and variances. Defaults to `"mean_var"`.
#' @param full_layer a bool indicating whether to output the predictions of all layers. Defaults to `FALSE`. Only used when `object` is a DGP or a linked (D)GP emulator.
#' @param sample_size the number of samples to draw for each given imputation if `method = "sampling"`. Defaults to `50`.
#' @param M the size of the conditioning set for the Vecchia approximation in the emulator prediction. Defaults to `50`. This argument is only used if the emulator `object`
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
#' * `r new_badge("updated")` If `object` is an instance of the `dgp` class:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains two
#'      matrices named `mean` for the predictive means and `var` for the predictive variances respectively. Each matrix has its rows corresponding to testing
#'      positions and columns corresponding to DGP global output dimensions (i.e., the number of GP/likelihood nodes in the final layer). If the likelihood node
#'      is categorical, the matrices contain the predictive means and variances of the class probabilities, with columns corresponding to different classes.
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains two
#'      sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L* (i.e., the number of layers)
#'      matrices named `layer1, layer2,..., layerL`. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      output dimensions (i.e., the number of GP/likelihood nodes from the associated layer). If the likelihood node is categorical, the matrices named `layerL`
#'      in both `mean` and `var` contain the predictive means and variances of the class probabilities, respectively, with columns corresponding to different classes.
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains *D* (i.e., the number
#'      of GP/likelihood nodes in the final layer) matrices named `output1, output2,..., outputD`. If the likelihood node in the final layer is categorical, `results`
#'      contains *D* matrices (where *D* is the number of classes) of sampled class probabilities, each named according to its corresponding class label. Each matrix in `results`
#'      has its rows corresponding to testing positions and columns corresponding to samples of size: `B * sample_size`, where `B` is the number of imputations specified in [dgp()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains *L* (i.e., the number
#'      of layers) sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents samples drawn from the GP/likelihood nodes in the corresponding layer,
#'      and contains *D* (i.e., the number of GP/likelihood nodes in the corresponding layer) matrices named `output1, output2,..., outputD`. If the likelihood node in the final
#'      layer is categorical, `layerL` contains *D* matrices (where *D* is the number of classes) of sampled class probabilities, each named according to its corresponding class
#'      label. Each matrix has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `B * sample_size`, where `B` is the number of imputations specified in [dgp()].
#' * If `object` is an instance of the `lgp` class:
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
#' The `results` slot will also include:
#' * the value of `M`, which represents the size of the conditioning set for the Vecchia approximation, if used, in the emulator prediction.
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
predict.dgp <- function(object, x, method = "mean_var", full_layer = FALSE, sample_size = 50, M = 50, cores = 1, chunks = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
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
    encoder <- object$constructor_obj$all_layer[[L]][[1]]$class_encoder
    index_to_label <- as.character(encoder$classes_)
  } else {
    is.categorical <- FALSE
  }


  if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)

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
    res <- object$emulator_obj$predict(x_unique, method, full_layer, sample_size, M)
  } else {
    res <- object$emulator_obj$ppredict(x_unique, method, full_layer, sample_size, M, chunks, cores)
  }

  if (method == 'mean_var'){
    if (full_layer) {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      LL <- length(named_res$mean)
      for (l in 1:LL) {
        named_res$mean[[l]] <- named_res$mean[[l]][rep,,drop=F]
        named_res$var[[l]] <- named_res$var[[l]][rep,,drop=F]
        names(named_res$mean)[l] <- paste('layer', l, sep="")
        names(named_res$var)[l] <- paste('layer', l, sep="")
      }
      if (is.categorical){
        if (length(index_to_label)==2){
          named_res$mean[[LL]] <- cbind(1-named_res$mean[[LL]], named_res$mean[[LL]])
          named_res$var[[LL]] <- cbind(named_res$var[[LL]], named_res$var[[LL]])
        }
        colnames(named_res$mean[[LL]]) <- index_to_label
        colnames(named_res$var[[LL]]) <- index_to_label
      }
    } else {
      named_res <- list("mean" = res[[1]][rep,,drop=F], "var" = res[[2]][rep,,drop=F])
      if (is.categorical){
        if (length(index_to_label)==2){
          named_res$mean <- cbind(1-named_res$mean, named_res$mean)
          named_res$var <- cbind(named_res$var, named_res$var)
        }
        colnames(named_res$mean) <- index_to_label
        colnames(named_res$var) <- index_to_label
      }
    }
  } else if (method == 'sampling') {
    if (full_layer) {
      named_res <- res
      LL <- length(named_res)
      for (l in 1:LL) {
        names(named_res)[l] <- paste('layer', l, sep="")
        for (k in 1:length(named_res[[l]])) {
          named_res[[l]][[k]] <- named_res[[l]][[k]][rep,,drop=F]
          if ( l == LL && is.categorical ) {
            if (length(index_to_label)==2){
              names(named_res[[l]])[k] <- index_to_label[k+1]
            } else {
              names(named_res[[l]])[k] <- index_to_label[k]
            }
          } else {
            names(named_res[[l]])[k] <- paste('output', k, sep="")
          }
        }
      }
      if (is.categorical && length(index_to_label)==2){
        named_res[[LL]][[2]] <- 1 - named_res[[LL]][[1]]
        names(named_res[[LL]])[2] <- index_to_label[1]
      }
    } else {
      named_res <- res
      for (l in 1:length(named_res)) {
        named_res[[l]] <- named_res[[l]][rep,,drop=F]
        if ( is.categorical ) {
          if (length(index_to_label)==2){
            names(named_res)[l] <- index_to_label[l+1]
          } else {
            names(named_res)[l] <- index_to_label[l]
          }
        } else {
          names(named_res)[l] <- paste('output', l, sep="")
        }
      }
      if (is.categorical && length(index_to_label)==2){
        named_res[[2]] <- 1 - named_res[[1]]
        names(named_res)[2] <- index_to_label[1]
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
#' @method predict lgp
#' @export
predict.lgp <- function(object, x, method = "mean_var", full_layer = FALSE, sample_size = 50, M = 50, cores = 1, chunks = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"lgp") ) stop("'object' must be an instance of the 'lgp' class.", call. = FALSE)

    if ( !("emulator_obj" %in% names(object)) ){
      stop("'object' is not activated for predictions. Please set `activate = TRUE` in `lgp()` to activate the emulator.", call. = FALSE)
    }
    if ( reticulate::py_is_null_xptr(object$emulator_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
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


  if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)

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
            emulator_id <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == l & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
            names(named_res$mean[[l]])[k] <- emulator_id
            names(named_res$var[[l]])[k] <- emulator_id
        }
      }
    } else {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      final_layer <- max(object$specs$metadata$Layer)
      for (k in 1:length(named_res$mean)) {
          emulator_id <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == final_layer & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
          names(named_res$mean)[k] <- emulator_id
          names(named_res$var)[k] <- emulator_id
      }
    }
  } else if (method == 'sampling') {
    if (full_layer) {
      named_res <- list()
      for (l in 1:length(res)) {
        name1 <- paste('layer', l, sep="")
        for (k in 1:length(res[[l]])){
            name2 <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == l & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
          for (n in 1:dim(res[[l]][[k]])[1]) {
            name3 <- paste('output', n, sep="")
            named_res[[name1]][[name2]][[name3]] <- res[[l]][[k]][n,,]
          }
        }
      }
    } else {
      named_res <- list()
      final_layer <- max(object$specs$metadata$Layer)
      for (k in 1:length(res)) {
          name1 <- subset(object$specs$metadata, object$specs$metadata[["Layer"]] == final_layer & object$specs$metadata[["Pos_in_Layer"]] == k)$Emulator
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
predict.gp <- function(object, x, method = "mean_var", sample_size = 50, M = 50, cores = 1, chunks = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) {
    if ( ncol(object$data$X)!=1 ){
      x <- matrix(x, nrow = 1)
    } else {
      x <- as.matrix(x)
    }
  }
  if ( ncol(x) != ncol(object$data$X) ) stop("'x' must have the same number of dimensions as the training input.", call. = FALSE)

  if ( method!='mean_var' & method!='sampling' ) stop("'method' can only be either 'mean_var' or 'sampling'.", call. = FALSE)

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

