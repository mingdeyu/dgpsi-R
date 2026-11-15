.onAttach <- function(...){
  # Retrieve Current Version
  this.version = utils::packageVersion("dgpsi")

  # Wesbite
  website <- if (grepl("[.-]9000$", this.version)) {
    "https://mingdeyu.github.io/dgpsi-R/dev/"
  } else {
    "https://mingdeyu.github.io/dgpsi-R/"
  }

  ## Print on Screen
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

  # echo output to screen
  packageStartupMessage("===========================================================================")
  packageStartupMessage("                 Welcome to the package 'dgpsi'", " v", this.version,"!")
  packageStartupMessage("                         Copyright (C) 2022-", this.year, sep="")
  packageStartupMessage("                  Maintainer: Deyu Ming (deyu.ming@ucl.ac.uk)")
  packageStartupMessage("             Package website: ", website)
  packageStartupMessage("===========================================================================")
  packageStartupMessage("---------------------------------------------------------------------------")
  packageStartupMessage(" Run init_py() from the package to install/activate the Python environment ")
  packageStartupMessage("---------------------------------------------------------------------------")
}
