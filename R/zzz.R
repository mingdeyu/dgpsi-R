.onAttach <- function(...){
  # Retrieve Current Version
  this.version = utils::packageVersion("dgpsi")

  ## Print on Screen
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

  # echo output to screen
  packageStartupMessage("============================================================================")
  packageStartupMessage("                 Welcome to the package 'dgpsi'", " v", this.version,"!")
  packageStartupMessage("                       Copyright (C) 2022-", this.year, sep="")
  packageStartupMessage("             Maintainer: Deyu Ming (deyu.ming.16@ucl.ac.uk)")
  packageStartupMessage("============================================================================")
  packageStartupMessage("----------------------------------------------------------------------------")
  packageStartupMessage(" Run a function from the package to install/activate the Python environment ")
  packageStartupMessage("----------------------------------------------------------------------------")
}
