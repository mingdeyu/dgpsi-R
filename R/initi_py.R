pkg.env <- new.env(parent = emptyenv())
pkg.env$dgpsi <- NULL
pkg.env$py_buildin <- NULL

#' Initialize the Python environment
#' @export
init_py <- function() {
  py_ver <- '3.9.13'
  dgpsi_ver <- 'dgpsi==2.1.1'
  env_name <- 'dgp_si_R_2_1_1'
  #Check if there is any conda binary installed, if not, request to install it.
  if (is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))){
    ans <- readline(prompt="I am unable to find a conda binary. Do you want me to install it for you? (Y/N) ")
    #If the user would like to have the conda binary to be installed
    if (ans == 'Y'){
      message("Installing the Conda binary...")
      reticulate::install_miniconda()
      conda_path <- reticulate::conda_binary()
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
    } else{
      stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-initialize the Python environment", call. = FALSE)
    }
  } else {
    conda_path <- reticulate::conda_binary()
    no_dgpsi <- inherits(tryCatch(reticulate::conda_python(envname = env_name, conda = conda_path), error = identity), "error")
    if (no_dgpsi){
      message("I am unable to find the required Python environment. It may be because it is your first time using the package or your conda binary has changed. I am (re)setting it for you now...")
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
    }
  }
  reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE)
  #py_buildin <<- reticulate::import_builtins()
  #dgpsi <<- reticulate::import("dgpsi")
  assign('dgpsi', reticulate::import("dgpsi"), pkg.env)
  assign('py_buildin', reticulate::import_builtins(), pkg.env)
  message("The Python environment is loaded.")
}

install_dgpsi <- function(env_name, py_ver, conda_path, dgpsi_ver) {
  message("Setting up the Python environment...")
  reticulate::conda_create(envname = env_name, python_version = py_ver, conda = conda_path)
  message("Installing the required Python packages...")
  if (Sys.info()[["sysname"]] == "Darwin" & Sys.info()[["machine"]] == "arm64"){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*accelerate"') , conda = conda_path)
  } else if (grepl("Intel",benchmarkme::get_cpu()$model_name)){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*mkl"') , conda = conda_path)
  } else {
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver) , conda = conda_path)
  }
}
