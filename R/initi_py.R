pkg.env <- new.env(parent = emptyenv())
pkg.env$dgpsi <- NULL
pkg.env$py_buildin <- NULL
pkg.env$np <- NULL

#' @title 'python' environment initialization
#'
#' @description This function initializes the 'python' environment for the package.
#'
#' @param py_ver a string that gives the 'python' version to be installed. If `py_ver = NULL`, the default 'python'
#'     version '3.9.13' will be installed.
#' @param dgpsi_ver a string that gives the 'python' version of 'dgpsi' to be used.
#'     If `dgpsi_ver = NULL`, the latest 'python' version of 'dgpsi' will be used.
#'
#' @return No return value, called to install required 'python' environment.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#'
#' @md
#' @export
init_py <- function(py_ver = NULL, dgpsi_ver = NULL) {
  if ( is.null(py_ver) ) py_ver <- '3.9.13'
  if ( is.null(dgpsi_ver) ) {
    dgpsi_ver <- 'dgpsi==2.1.5'
    env_name <- 'dgp_si_R_2_1_5'
  } else {
    env_name <- paste('dgp_si_R_', gsub(".", "_", dgpsi_ver,fixed=TRUE), sep = "")
    dgpsi_ver <- paste('dgpsi==', dgpsi_ver, sep = "")
  }
  #Check if there is any conda binary installed, if not, request to install it.
  if (is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))){
    ans <- readline(prompt="I am unable to find a conda binary. Do you want me to install it for you? (Y/N) ")
    #If the user would like to have the conda binary to be installed
    if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
      message("Installing the Conda binary...")
      reticulate::install_miniconda()
      conda_path <- reticulate::conda_binary()
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
    } else{
      stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-initialize the Python environment.", call. = FALSE)
    }
  } else {
    conda_path <- reticulate::conda_binary()
    no_dgpsi <- inherits(tryCatch(reticulate::conda_python(envname = env_name, conda = conda_path), error = identity), "error")
    if (no_dgpsi){
      if (any(grepl('^dgp_si_R', reticulate::conda_list(conda = conda_path)$name))){
        conda_list <- reticulate::conda_list(conda = conda_path)$name
        dgpsi_list <- conda_list[grepl('^dgp_si_R', conda_list)]
        ans <- readline(prompt="I found Python environment(s) for other versions of the package. Do you want me to remove them? (Y/N) ")
        if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
          message(sprintf("Removing Python environment(s): %s.\n", paste(dgpsi_list, collapse = ', ')))
          for (item in dgpsi_list) {
            reticulate::conda_remove(envname = item, conda = conda_path)
          }
          message("Done.")
        }
        install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
      } else {
        ans <- readline(prompt="Is this your first time installing the package? (Y/N) ")
        if ( tolower(ans)=='n'|tolower(ans)=='no' ){
              message("I am unable to find the required Python environment. It may be because your conda binary has changed. I am re-setting it for you now ...")
              install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
              } else {
                install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
              }
      }
    }
  }

  message("Connecting to Python ...", appendLF = FALSE)
  warning_error_handler(with_warning_handler(reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE)))
  message(" done")

  message("Importing required Python modules ...", appendLF = FALSE)
  assign('dgpsi', reticulate::import("dgpsi"), pkg.env)
  assign('py_buildin', reticulate::import_builtins(), pkg.env)
  assign('np', reticulate::import("numpy"), pkg.env)
  message(" done")
  Sys.sleep(0.5)
  message("The Python environment for 'dgpsi' is successfully loaded.")
}

install_dgpsi <- function(env_name, py_ver, conda_path, dgpsi_ver) {
  message(sprintf("Setting up the Python environment for %s ...\n", dgpsi_ver))
  reticulate::conda_create(envname = env_name, python_version = py_ver, conda = conda_path)
  message("Installing the required Python packages ...")
  if (Sys.info()[["sysname"]] == "Darwin" & Sys.info()[["machine"]] == "arm64"){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*accelerate"') , conda = conda_path)
  } else if (grepl("Intel",benchmarkme::get_cpu()$model_name)){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*mkl"') , conda = conda_path)
  } else {
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver) , conda = conda_path)
  }
}

with_warning_handler <- function(...)
{
  Sys.unsetenv("RETICULATE_PYTHON")
  withCallingHandlers(..., warning = function(w)
  { condition <- conditionMessage(w)
  reg1 <- "Previous request to"
  reg2 <- "will be ignored. It is superseded by request to"
  if(grepl(reg1, condition) & grepl(reg2, condition)) invokeRestart("muffleWarning")
  })
}

warning_error_handler <- function(...){
  tryCatch(...,
           error = function(r)
           { condition <- conditionMessage(r)
           reg <- "failed to initialize requested version of Python"
           if(grepl(reg, condition)) {
             cat("NOTE: please clear your R workspace and delete the workspace image file '.RData' before restarting the R session.")
           } else {
             message(paste("ERROR in", condition))
           }
           })
}
