pack_bundle <- function(object, first_time, N, eval, new_wave, rmse_records, freq, N_acq_ind, target, type, X, Y, y_cand, n_wave = NULL, design_info = NULL, x_test = NULL, y_test = NULL, istarget = NULL, target_points = NULL){

  if ( isTRUE(first_time) ){
    object$design <- list()
    object[["design"]][["wave1"]][["N"]] <- N
    if ( is.null(eval) ){
      object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
    } else {
      object[["design"]][["wave1"]][["metric"]] <- unname(rmse_records)
    }
    object[["design"]][["wave1"]][["freq"]] <- freq[2]
    object[["design"]][["wave1"]][["enrichment"]] <- unname(N_acq_ind)
    if( !is.null(target) ) {
      object[["design"]][["wave1"]][["target"]] <- target
      object[["design"]][["wave1"]][["reached"]] <- if(length(target)==1) { all(istarget) } else { istarget }
    }
    object[["design"]][["type"]] <- type
    if ( identical(type, 'oos') ){
      object[["design"]][["x_test"]] <- x_test
      object[["design"]][["y_test"]] <- y_test
    }
    object[['data']][['X']] <- X
    object[['data']][['Y']] <- Y
  } else {
    object$design <- design_info
    if (new_wave){
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- unname(N_acq_ind)
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- if(length(target)==1) { all(istarget) } else { istarget }
      }
    } else {
      object[["design"]][[paste('wave', n_wave, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]], unname(rmse_records))
      } else {
        object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]] <- rbind( object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]], unname(rmse_records))
      }
      object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]], unname(N_acq_ind))
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave, sep="")]][["reached"]] <- if(length(target)==1) { all(istarget) } else { istarget }
      }
    }
    object[['data']][['X']] <- X
    object[['data']][['Y']] <- Y
  }

  if ( is.null(y_cand) ) {
    if (!all(sapply(target_points, is.null))) object[["design"]][['exclusion']] <- target_points
  }
  return(object)
}


pack_dgp <- function(object, first_time, N, eval, new_wave, rmse_records, freq, N_acq, n_dim_Y, target, type, y_cand, n_wave = NULL, design_info = NULL, x_test = NULL, y_test = NULL, istarget = NULL, target_points = NULL){
  if ( isTRUE(first_time) ){
    object$design <- list()
    object[["design"]][["wave1"]][["N"]] <- N
    if ( is.null(eval) ){
      object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
    } else {
      object[["design"]][["wave1"]][["metric"]] <- unname(rmse_records)
    }
    object[["design"]][["wave1"]][["freq"]] <- freq[2]
    if (inherits(object,"bundle")){
      object[["design"]][["wave1"]][["enrichment"]] <- matrix(rep(N_acq, n_dim_Y), nrow = length(N_acq), byrow = F)
    } else {
      object[["design"]][["wave1"]][["enrichment"]] <- N_acq
    }
    if( !is.null(target) ) {
      object[["design"]][["wave1"]][["target"]] <- target
      if (inherits(object,"bundle")){
        object[["design"]][["wave1"]][["reached"]] <- if(length(target)==1) {all(istarget)} else {istarget}
      } else {
        object[["design"]][["wave1"]][["reached"]] <- istarget
      }
    }
    object[["design"]][["type"]] <- type
    if ( identical(type, 'oos') ){
      object[["design"]][["x_test"]] <- x_test
      object[["design"]][["y_test"]] <- y_test
    }
  } else {
    object$design <- design_info
    if (new_wave){
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
      if (inherits(object,"bundle")){
        for (w in 1:n_wave){
          current_enrichments <- object[["design"]][[paste('wave', w, sep="")]][["enrichment"]]
          object[["design"]][[paste('wave', w, sep="")]][["enrichment"]] <- matrix(rep(current_enrichments, n_dim_Y), nrow = length(current_enrichments), byrow = F)
        }
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- matrix(rep(N_acq, n_dim_Y), nrow = length(N_acq), byrow = F)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- N_acq
      }
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["target"]] <- target
        if (inherits(object,"bundle")){
          for (w in 1:n_wave){
            object[["design"]][[paste('wave', w, sep="")]][["reached"]] <- rep(object[["design"]][[paste('wave', w, sep="")]][["reached"]], length(target))
          }
          object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- if(length(target)==1) {all(istarget)} else {istarget}
        } else {
          object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- istarget
        }
      }
    } else {
      object[["design"]][[paste('wave', n_wave, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]], unname(rmse_records))
      } else {
        object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]], unname(rmse_records))
      }
      if (inherits(object,"bundle")){
        for (w in 1:n_wave){
          current_enrichments <- object[["design"]][[paste('wave', w, sep="")]][["enrichment"]]
          object[["design"]][[paste('wave', w, sep="")]][["enrichment"]] <- matrix(rep(current_enrichments, n_dim_Y), nrow = length(current_enrichments), byrow = F)
        }
        object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]], matrix(rep(N_acq, n_dim_Y), nrow = length(N_acq), byrow = F))
      } else {
        object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]] <- c(object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]], N_acq)
      }
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave, sep="")]][["target"]] <- target
        if (inherits(object,"bundle")){
          for (w in 1:(n_wave-1)){
            object[["design"]][[paste('wave', w, sep="")]][["reached"]] <- rep(object[["design"]][[paste('wave', w, sep="")]][["reached"]], length(target))
          }
          object[["design"]][[paste('wave', n_wave, sep="")]][["reached"]] <- if(length(target)==1) {all(istarget)} else {istarget}
        } else {
          object[["design"]][[paste('wave', n_wave, sep="")]][["reached"]] <- istarget
        }
      }
    }
  }

  if ( is.null(y_cand) ) {
    if (!is.null(target_points)) object[["design"]][['exclusion']] <- target_points
  }
  return(object)
}

pack_gp <- function(object, first_time, N, eval, new_wave, rmse_records, freq, N_acq, target, type, y_cand, n_wave = NULL, design_info = NULL, x_test = NULL, y_test = NULL, istarget = NULL, target_points = NULL){

  if ( isTRUE(first_time) ){
    object$design <- list()
    object[["design"]][["wave1"]][["N"]] <- N
    if ( is.null(eval) ){
      object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
    } else {
      object[["design"]][["wave1"]][["metric"]] <- unname(rmse_records)
    }
    object[["design"]][["wave1"]][["freq"]] <- freq[2]
    object[["design"]][["wave1"]][["enrichment"]] <- N_acq
    if( !is.null(target) ) {
      object[["design"]][["wave1"]][["target"]] <- target
      object[["design"]][["wave1"]][["reached"]] <- istarget
    }
    object[["design"]][["type"]] <- type
    if ( identical(type, 'oos') ){
      object[["design"]][["x_test"]] <- x_test
      object[["design"]][["y_test"]] <- y_test
    }
  } else {
    object$design <- design_info
    if (new_wave){
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <-  N_acq
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- istarget
      }
    } else {
      object[["design"]][[paste('wave', n_wave, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]], unname(rmse_records))
      } else {
        object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]] <- rbind(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]], unname(rmse_records))
      }
      object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]] <- c(object[["design"]][[paste('wave', n_wave, sep="")]][["enrichment"]], N_acq)
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave, sep="")]][["reached"]] <- istarget
      }
    }
  }

  if ( is.null(y_cand) ) {
    if (!is.null(target_points)) object[["design"]][['exclusion']] <- target_points
  }

  return(object)
}

