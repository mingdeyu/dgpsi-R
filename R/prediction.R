#' @title Predictions from GP, DGP, or linked (D)GP emulators
#'
#' @description This function implements prediction from GP, DGP, or linked (D)GP emulators.
#'
#' @param obj a GP, DGP, or linked (D)GP object produced by [gp()] (after applying [train()]), [emulator()], or [lgp()].
#' @param x a matrix where each row is an input testing data point and each column is an input dimension.
#' @param method the prediction approach: mean-variance (`"mean_var"`) or sampling (`"sampling"`) approach. Defaults to `"mean_var"`.
#' @param full_layer whether to output the predictions of all layers. Defaults to `FALSE`. Only used when `obj` is for DGP and linked (D)GP.
#' @param sample_size the number of samples to draw for each given imputation if `method = "sampling"`. Defaults to `50`.
#'
#' @return
#' * **(GP case)** If `obj` is produced by [train()] on a object built with [gp()]:
#'   1. if `method = "mean_var"`: a named list is returned that contains two matrices named `mean` for the predictive means
#'      and `var` for the predictive variances. Each matrix has only one column with its rows
#'     corresponding to testing positions (i.e., row of `x`).
#'   2. if `method = "sampling"`: a matrix is returned whose rows corresponding to testing positions and columns corresponding
#'      to `sample_size` number of samples drawn from the predictive distribution of GP.
#' * **(DGP case)** If `obj` is produced by [emulator()]:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: a named list is returned. The list contains two matrices named `mean`
#'      for the predictive means and `var` for the predictive variances respectively. Each matrix has its rows corresponding to testing
#'      positions and columns corresponding to DGP global output dimensions (i.e., the number of GP/likelihood nodes in the final layer).
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: a named list is returned. The list contains two sub-lists named `mean`
#'      for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L* (i.e., the number of layers)
#'      matrices named `layer1, layer2,..., layerL`. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      output dimensions (i.e., the number of GP nodes from the associated layer and in case of the final layer, it may be the number of the
#'      likelihood nodes).
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: a named list is returned. The list contains *D* (i.e., the number of GP/likelihood
#'      nodes in the final layer) matrices named `output1, output2,..., outputD`. Each matrix has its rows corresponding to testing positions and
#'      columns corresponding to samples of size: `N * sample_size`, where `N` is the number of imputations specified in [emulator()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: a named list is returned. The list contains *L* (i.e., the number of layers)
#'      sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents samples drawn from the GPs/likelihoods in the corresponding layers,
#'      and contains *D* (i.e., the number of GP nodes in the corresponding layer or likelihood nodes in the final layer)
#'      matrices named `output1, output2,..., outputD`. Each matrix gives samples of the output from one of *D* GPs/likelihoods at the
#'      testing positions, and has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `N * sample_size`, where `N` is the number of imputations specified in [emulator()].
#' * **(Linked (D)GP case)** If `obj` is produced by [lgp()]:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: a named list is returned. The list contains two sub-lists named `mean`
#'      for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *M* number (same number
#'      of computer models in the final layer of the system) of matrices named `emulator1, emulator2,..., emulatorM`. Each matrix has
#'      its rows corresponding to global testing positions and columns corresponding to output dimensions of the associated emulator
#'      in the final layer.
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: a named list is returned. The list contains two sub-lists named `mean`
#'      for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L* (i.e., the number of layers
#'      of the emulated system) components named `layer1, layer2,..., layerL`. Each component represents a layer and contains *M* number
#'      (same number of emulators in the corresponding layer of the system) of matrices named `emulator1, emulator2,..., emulatorM`. Each
#'      matrix has its rows corresponding to global testing positions and columns corresponding to output dimensions of the associated
#'      GP/DGP emulator in the corresponding layer.
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: a named list is returned. The list contains *M* number (same number of emulators
#'      in the final layer of the system) of sub-lists named `emulator1, emulator2,..., emulatorM`. Each sub-list corresponds to an emulator
#'      in the final layer, and contains *D* matrices, named `output1, output2,..., outputD`, that correspond to the output
#'      dimensions of the GP/DGP emulator. Each matrix has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `N * sample_size`, where `N` is the number of imputations specified in [emulator()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: a named list is returned. The list contains *L* (i.e., the number of layers of
#'      the emulated system) sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents a layer and contains *M* number (same
#'      number of emulators in the corresponding layer of the system) of components named `emulator1, emulator2,..., emulatorM`. Each component
#'      corresponds to an emulator in the associated layer, and contains *D* matrices, named `output1, output2,..., outputD`, that correspond to
#'      the output dimensions of the GP/DGP emulator. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      samples of size: `N * sample_size`, where `N` is the number of imputations specified in [emulator()].
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
predict <- function(obj, x, method = 'mean_var',full_layer = FALSE, sample_size = 50) {
  sample_size <- as.integer(sample_size)
  if (pkg.env$py_buildin$type(obj)$'__name__' == 'gp'){
     res <- obj$predict(x, method, sample_size)
     if (method == 'mean_var'){
       named_res <- list("mean" = res[[1]], "var" = res[[2]])
     } else if (method == 'sampling'){
       named_res <- res
       }
  } else if (pkg.env$py_buildin$type(obj)$'__name__' == 'emulator') {
    res <- obj$predict(x, method, full_layer, sample_size)
    if (method == 'mean_var'){
      if (full_layer) {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (l in 1:length(named_res$mean)) {
          names(named_res$mean)[l] <- paste('layer', l, sep="")
          names(named_res$var)[l] <- paste('layer', l, sep="")
        }
      } else {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
      }
    } else if (method == 'sampling') {
      if (full_layer) {
        named_res <- res
        for (l in 1:length(named_res)) {
          names(named_res)[l] <- paste('layer', l, sep="")
          for (k in 1:length(named_res[[l]])) {
            names(named_res[[l]])[k] <- paste('output', k, sep="")
            }
        }
      } else {
        named_res <- res
        for (l in 1:length(named_res)) {
          names(named_res)[l] <- paste('output', l, sep="")
        }
      }
    }
  } else if (pkg.env$py_buildin$type(obj)$'__name__' == 'lgp') {
    res <- obj$predict(x, method, full_layer, sample_size)
    if (method == 'mean_var'){
      if (full_layer) {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (l in 1:length(named_res$mean)) {
          names(named_res$mean)[l] <- paste('layer', l, sep="")
          names(named_res$var)[l] <- paste('layer', l, sep="")
          for (k in 1:length(named_res$mean[[l]])) {
              names(named_res$mean[[l]])[k] <- paste('emulator', k, sep="")
              names(named_res$var[[l]])[k] <- paste('emulator', k, sep="")
          }
        }
      } else {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (k in 1:length(named_res$mean)) {
          names(named_res$mean)[k] <- paste('emulator', k, sep="")
          names(named_res$var)[k] <- paste('emulator', k, sep="")
        }
      }
    } else if (method == 'sampling') {
      if (full_layer) {
        named_res <- list()
        for (l in 1:length(res)) {
          name1 <- paste('layer', l, sep="")
          for (k in 1:length(res[[l]])){
            name2 <- paste('emulator', k, sep="")
            for (n in 1:dim(res[[l]][[k]])[1]) {
              name3 <- paste('output', n, sep="")
              named_res[[name1]][[name2]][[name3]] <- res[[l]][[k]][n,,]
            }
          }
        }
      } else {
        named_res <- list()
        for (k in 1:length(res)) {
          name1 <- paste('emulator', k, sep="")
          for (n in 1:dim(res[[k]])[1]) {
              name2 <- paste('output', n, sep="")
              named_res[[name1]][[name2]] <- res[[k]][n,,]
          }
        }
      }
    }
  }
  return(named_res)
}


#' @title Multi-core predictions from GP, DGP, or linked (D)GP emulators
#'
#' @description This function implements multi-core prediction from GP, DGP, or linked (D)GP emulators.
#'
#' @param obj same as that in [predict()].
#' @param x same as that in [predict()].
#' @param method same as that in [predict()].
#' @param full_layer same as that in [predict()].
#' @param sample_size same as that in [predict()].
#' @param chunk_num the number of chunks that the testing input matrix `x` will be divided into.
#'     Defaults to `NULL`. If not specified, the number of chunks is set to `core_num`.
#' @param core_num the number of cores/workers to be used. Defaults to `NULL`. If not specified,
#'     the number of cores is set to `(max physical cores available - 1)`.
#'
#' @return Same as that in [predict()].
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
ppredict <- function(obj, x, method = 'mean_var', full_layer = FALSE, sample_size = 50, chunk_num = NULL, core_num = NULL) {
  sample_size <- as.integer(sample_size)
  if(!is.null(chunk_num)){
    chunk_num <- as.integer(chunk_num)
  }
  if(!is.null(core_num)){
    core_num <- as.integer(core_num)
  }

  if (pkg.env$py_buildin$type(obj)$'__name__' == 'gp'){
    res <- obj$ppredict(x, method, sample_size, chunk_num, core_num)
    if (method == 'mean_var'){
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
    } else if (method == 'sampling'){
      named_res <- res
    }
  } else if (pkg.env$py_buildin$type(obj)$'__name__' == 'emulator') {
    res <- obj$ppredict(x, method, full_layer, sample_size, chunk_num, core_num)
    if (method == 'mean_var'){
      if (full_layer) {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (l in 1:length(named_res$mean)) {
          names(named_res$mean)[l] <- paste('layer', l, sep="")
          names(named_res$var)[l] <- paste('layer', l, sep="")
        }
      } else {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
      }
    } else if (method == 'sampling') {
      if (full_layer) {
        named_res <- res
        for (l in 1:length(named_res)) {
          names(named_res)[l] <- paste('layer', l, sep="")
          for (k in 1:length(named_res[[l]])) {
            names(named_res[[l]])[k] <- paste('output', k, sep="")
          }
        }
      } else {
        named_res <- res
        for (l in 1:length(named_res)) {
          names(named_res)[l] <- paste('output', l, sep="")
        }
      }
    }
  } else if (pkg.env$py_buildin$type(obj)$'__name__' == 'lgp') {
    res <- obj$ppredict(x, method, full_layer, sample_size, chunk_num, core_num)
    if (method == 'mean_var'){
      if (full_layer) {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (l in 1:length(named_res$mean)) {
          names(named_res$mean)[l] <- paste('layer', l, sep="")
          names(named_res$var)[l] <- paste('layer', l, sep="")
          for (k in 1:length(named_res$mean[[l]])) {
            names(named_res$mean[[l]])[k] <- paste('emulator', k, sep="")
            names(named_res$var[[l]])[k] <- paste('emulator', k, sep="")
          }
        }
      } else {
        named_res <- list("mean" = res[[1]], "var" = res[[2]])
        for (k in 1:length(named_res$mean)) {
          names(named_res$mean)[k] <- paste('emulator', k, sep="")
          names(named_res$var)[k] <- paste('emulator', k, sep="")
        }
      }
    } else if (method == 'sampling') {
      if (full_layer) {
        named_res <- list()
        for (l in 1:length(res)) {
          name1 <- paste('layer', l, sep="")
          for (k in 1:length(res[[l]])){
            name2 <- paste('emulator', k, sep="")
            for (n in 1:dim(res[[l]][[k]])[1]) {
              name3 <- paste('output', n, sep="")
              named_res[[name1]][[name2]][[name3]] <- res[[l]][[k]][n,,]
            }
          }
        }
      } else {
        named_res <- list()
        for (k in 1:length(res)) {
          name1 <- paste('emulator', k, sep="")
          for (n in 1:dim(res[[k]])[1]) {
            name2 <- paste('output', n, sep="")
            named_res[[name1]][[name2]] <- res[[k]][n,,]
          }
        }
      }
    }
  }
  return(named_res)
}
