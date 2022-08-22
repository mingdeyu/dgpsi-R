pkg.env <- new.env(parent = emptyenv())
pkg.env$dgpsi <- NULL
pkg.env$py_buildin <- NULL
pkg.env$sys <- NULL

#' Initialize the Python environment
#' @note See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#'
#' @md
#' @export
init_py <- function() {
  py_ver <- '3.9.13'
  dgpsi_ver <- 'dgpsi==2.1.2'
  env_name <- 'dgp_si_R_2_1_2'
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
      stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-initialize the Python environment.", call. = FALSE)
    }
  } else {
    conda_path <- reticulate::conda_binary()
    no_dgpsi <- inherits(tryCatch(reticulate::conda_python(envname = env_name, conda = conda_path), error = identity), "error")
    if (no_dgpsi){
      if (any(grepl('^dgp_si_R', reticulate::conda_list(conda = conda_path)$name))){
        conda_list <- reticulate::conda_list(conda = conda_path)$name
        dgpsi_list <- conda_list[grepl('^dgp_si_R', conda_list)]
        ans <- readline(prompt="I found Python environment(s) for old package versions. Do you want me to remove them for you? (Y/N) ")
        if (ans == 'Y'){
          message(sprintf("Removing Python environment(s): %s.\n", paste(dgpsi_list, collapse = ', ')))
          for (item in dgpsi_list) {
            reticulate::conda_remove(envname = item, conda = conda_path)
          }
          message("Done.")
        }
        install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
      } else {
        ans <- readline(prompt="Is this your first time installing the package? (Y/N) ")
        if (ans == 'N'){
              message("I am unable to find the required Python environment. It may be because your conda binary has changed. I am re-setting it for you now ...")
              install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
              } else {
                install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
              }
      }
    }
  }

  Sys.unsetenv("RETICULATE_PYTHON")
  with_warning_handler(reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE))

  assign('dgpsi', reticulate::import("dgpsi"), pkg.env)
  assign('py_buildin', reticulate::import_builtins(), pkg.env)
  assign('sys', reticulate::import("sys"), pkg.env)
  message("The Python environment for dgpsi is successfully loaded.")
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
  withCallingHandlers(..., warning = function(w)
  { condition <- conditionMessage(w)
  reg1 <- "Previous request to"
  reg2 <- "will be ignored. It is superseded by request to"
  if(grepl(reg1, condition) & grepl(reg2, condition)) invokeRestart("muffleWarning")
  })
}
