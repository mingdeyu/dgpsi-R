% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{unpack}
\alias{unpack}
\title{Unpack a bundle of (D)GP emulators}
\usage{
unpack(object)
}
\arguments{
\item{object}{an instance of the class \code{bundle}.}
}
\value{
A named list that contains individual emulators (named \verb{emulator1,...,emulatorS}) packed in \code{object},
where \code{S} is the number of emulators in \code{object}.
}
\description{
This function unpacks a bundle of (D)GP emulators safely so any further manipulations of unpacked individual emulators
will not impact the ones in the bundle.
}
\details{
See further examples and tutorials at \url{https://mingdeyu.github.io/dgpsi-R/}.
}
\examples{
\dontrun{

# See pack() for an example.
}
}
