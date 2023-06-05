pkg.env <- new.env(parent = emptyenv())
pkg.env$dgpsi <- NULL
pkg.env$py_buildin <- NULL
pkg.env$np <- NULL
pkg.env$copy <- NULL
pkg.env$py_gc <- NULL

#' @title 'python' environment initialization
#'
#' @description This function initializes the 'python' environment for the package.
#'
#' @param py_ver a string that gives the 'python' version to be installed. If `py_ver = NULL`, the default 'python'
#'     version '3.9.13' will be installed.
#' @param dgpsi_ver a string that gives the 'python' version of 'dgpsi' to be used. If `dgpsi_ver = NULL`,
#' * the latest 'python' version of 'dgpsi' will be used, if the package is installed from CRAN;
#' * the development 'python' version of 'dgpsi' will be used, if the package is installed from GitHub.
#' @param reinstall a bool that indicates whether to reinstall the 'python' version of 'dgpsi' specified
#'    in `dgpsi_ver` if it has already been installed. This argument is useful when the development version
#'    of the R package is installed and one may want to regularly update the development 'python' version
#'    of 'dgpsi'. Defaults to `FALSE`.
#' @param uninstall a bool that indicates whether to uninstall the 'python' version of 'dgpsi' specified
#'    in `dgpsi_ver` if it has already been installed. This argument is useful when the 'python' environment
#'    is corrupted and one wants to completely uninstall and reinstall it. Defaults to `FALSE`.
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
init_py <- function(py_ver = NULL, dgpsi_ver = NULL, reinstall = FALSE, uninstall = FALSE) {
  if ( is.null(py_ver) ) py_ver <- '3.9.13'
  if ( is.null(dgpsi_ver) ) {
    ##For devel version
    #dgpsi_ver <- c('cython>=0.29.30', 'dill>=0.3.2', 'jupyter>=1.0.0', 'matplotlib-base>=3.2.1', 'numba >=0.51.2',
    #               'numpy >=1.18.2', 'pathos >=0.2.9', 'psutil >=5.8.0', 'pybind11 >=2.10.0', 'pythran >=0.11.0',
    #               'scikit-build >=0.15.0', 'scikit-learn >=0.22.0', 'scipy >=1.4.1', 'tqdm >=4.50.2', 'tabulate >=0.8.7')
    #env_name <- 'dgp_si_R_2_1_6_9000'
    ##For release version
    dgpsi_ver <- 'dgpsi==2.2.0'
    env_name <- 'dgp_si_R_2_2_0'
  } else {
    env_name <- paste('dgp_si_R_', gsub(".", "_", dgpsi_ver,fixed=TRUE), sep = "")
    dgpsi_ver <- paste('dgpsi==', dgpsi_ver, sep = "")
  }
  #Check if there is any conda binary installed, if not, request to install it.
  restart <- FALSE
  if (is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))){
    ans <- readline(prompt="I am unable to find a conda binary. Do you want me to install it for you? (Y/N) ")
    #If the user would like to have the conda binary to be installed
    if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
      message("Installing the Conda binary...")
      reticulate::install_miniconda()
      conda_path <- reticulate::conda_binary()
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
      restart <- TRUE
    } else{
      stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-initialize the Python environment.", call. = FALSE)
    }
  } else {
    conda_path <- reticulate::conda_binary()
    no_dgpsi <- inherits(tryCatch(reticulate::conda_python(envname = env_name, conda = conda_path), error = identity), "error")
    if (no_dgpsi){
      if (any(grepl('^dgp_si_R', reticulate::conda_list(conda = conda_path)$name))){
        #conda_list <- reticulate::conda_list(conda = conda_path)$name
        #dgpsi_list <- conda_list[grepl('^dgp_si_R', conda_list)]
        #cat("I found Python environment(s) for other versions of the package.")
        #ans <- readline(prompt="Do you want me to remove them? (Y/N) ")
        #if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
        #  message(sprintf("Removing Python environment(s): %s.\n", paste(dgpsi_list, collapse = ', ')))
        #  for (item in dgpsi_list) {
        #    reticulate::conda_remove(envname = item, conda = conda_path)
        #  }
        #  message("Done.")
        #}
        install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
        restart <- TRUE
      } else {
        ans <- readline(prompt="Is this your first time using the package? (Y/N) ")
        if ( tolower(ans)=='n'|tolower(ans)=='no' ){
              message("I am unable to find the required Python environment. It may be because your conda binary has changed.")
              cat("I am re-setting it for you now ...")
              install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
              restart <- TRUE
              } else {
                install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
                restart <- TRUE
              }
      }
    } else {
      if ( uninstall ){
        reticulate::conda_remove(envname = env_name, conda = conda_path)
        restart <- TRUE
        message("Uninstallation finished. Please restart R and run 'init_py()' to reinstall the Python environment.")
      } else {
        if (isTRUE(reinstall)) {
          install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver, reinsatll = TRUE)
          #if (grepl('9000',env_name)) {
          #  reticulate::conda_install(envname = env_name, packages = c("git+https://github.com/mingdeyu/DGP.git") , conda = conda_path, pip = TRUE, pip_options = c('--no-deps', '--force-reinstall'))
          #} else {
          #  reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver) , conda = conda_path)
          #}
          #if (Sys.info()[["sysname"]] == 'Linux'){
          #  libstdc_path <- paste(gsub("bin.*$", "", conda_path), 'envs/', env_name, '/lib/libstdc++.so.6.0.30', sep='')
          #  libstdc_sys_path <- "/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
          #  system(paste("sudo rm",libstdc_sys_path))
          #  system(paste("sudo ln -s", libstdc_path, libstdc_sys_path))
          #}
          #message("Installation finished. Please restart R.")
          restart <- TRUE
        }
      }
    }
  }

  if ( isFALSE(restart) ){
    message("Connecting to Python ...", appendLF = FALSE)
    warning_error_handler(with_warning_handler(reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE)))
    message(" done")

    message("Importing required Python modules ...", appendLF = FALSE)
    assign('dgpsi', reticulate::import("dgpsi"), pkg.env)
    assign('py_buildin', reticulate::import_builtins(), pkg.env)
    assign('np', reticulate::import("numpy"), pkg.env)
    assign('copy', reticulate::import("copy"), pkg.env)
    assign('py_gc', reticulate::import("gc"), pkg.env)
    message(" done")
    Sys.sleep(0.5)
    message("The Python environment for 'dgpsi' is successfully loaded.")
  }
}

install_dgpsi <- function(env_name, py_ver, conda_path, dgpsi_ver, reinsatll = FALSE) {
  if (!reinsatll) message(sprintf("Setting up the Python environment for %s ...\n", dgpsi_ver))
  if (!reinsatll) reticulate::conda_create(envname = env_name, python_version = py_ver, conda = conda_path)
  if (reinsatll) {
    message("Re-installing the required Python packages ...")
  } else {
    message("Installing the required Python packages ...")
  }
  if (Sys.info()[["sysname"]] == "Darwin" & Sys.info()[["machine"]] == "arm64"){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*accelerate"') , conda = conda_path)
  } else if (grepl("Intel",benchmarkme::get_cpu()$model_name)){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*mkl"') , conda = conda_path)
    reticulate::conda_install(envname = env_name, packages = c('icc_rt') , channel = c("numba"), conda = conda_path)
  } else {
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver) , conda = conda_path)
  }
  if (grepl('9000',env_name)) {
    if (reinsatll) {
      reticulate::conda_install(envname = env_name, packages = c("git+https://github.com/mingdeyu/DGP.git") , conda = conda_path, pip = TRUE, pip_options = c('--no-deps', '--force-reinstall'))
    } else {
      reticulate::conda_install(envname = env_name, packages = c("git+https://github.com/mingdeyu/DGP.git") , conda = conda_path, pip = TRUE)
    }
  }
  if (Sys.info()[["sysname"]] == 'Linux'){
    libstdc_path <- paste(gsub("bin.*$", "", conda_path), 'envs/', env_name, '/lib/libstdc++.so.6.0.30', sep='')
    libstdc_sys_path <- "/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
    system(paste("sudo rm",libstdc_sys_path))
    system(paste("sudo ln -s", libstdc_path, libstdc_sys_path))
  }
  message("Installation finished. Please restart R.")
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
