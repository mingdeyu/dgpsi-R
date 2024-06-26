% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{set_vecchia}
\alias{set_vecchia}
\title{Add or remove the Vecchia approximation}
\usage{
set_vecchia(object, vecchia = TRUE, M = 25, ord = NULL)
}
\arguments{
\item{object}{an instance of the S3 class \code{gp}, \code{dgp}, or \code{lgp}.}

\item{vecchia}{a boolean or a list of booleans to indicate the addition or removal of the Vecchia approximation:
\itemize{
\item if \code{object} is an instance of the \code{gp} or \code{dgp} class, \code{vecchia} is a boolean that indicates
either addition (\code{vecchia = TRUE}) or removal (\code{vecchia = FALSE}) of the Vecchia approximation from \code{object}.
\item if \code{object} is an instance of the \code{lgp} class, \code{x} can be a boolean or a list of booleans:
\itemize{
\item if \code{vecchia} is a boolean, it indicates either addition (\code{vecchia = TRUE}) or removal (\code{vecchia = FALSE}) of
the Vecchia approximation from all individual (D)GP emulators contained in \code{object}.
\item if \code{vecchia} is a list of booleans, it should have same shape as \code{struc} that was supplied to \code{\link[=lgp]{lgp()}}. Each boolean
in the list indicates if the corresponding (D)GP emulator contained in \code{object} shall have the Vecchia approximation
added or removed.
}
}}

\item{M}{the size of the conditioning set for the Vecchia approximation in the (D)GP emulator training. Defaults to \code{25}.}

\item{ord}{an R function that returns the ordering of the input to the (D)GP emulator for the Vecchia approximation. The
function must satisfy the following basic rules:
\itemize{
\item the first argument represents the lengthscale-scaled input to the GP emulator or the lengthscale-scaled input to a GP node
of the DGP emulator.
\item the output of the function is a vector of indices that gives the ordering of the input to the GP emulator or the input to
the GP nodes of the DGP emulator.
}

If \code{ord = NULL}, the default random ordering is used. Defaults to \code{NULL}.}
}
\value{
An updated \code{object} with the Vecchia approximation either added or removed.
}
\description{
This function adds or removes the Vecchia approximation from a GP, DGP or linked (D)GP emulator
constructed by \code{\link[=gp]{gp()}}, \code{\link[=dgp]{dgp()}} or \code{\link[=lgp]{lgp()}}.
}
\details{
See further examples and tutorials at \url{https://mingdeyu.github.io/dgpsi-R/}.
}
\note{
This function is useful for quickly switching between Vecchia and non-Vecchia approximations for an existing emulator
without the need to reconstruct the emulator. If the emulator was built without the Vecchia approximation, the function
can add it, and if the emulator was built with the Vecchia approximation, the function can remove it. If the current
state already matches the requested state, the emulator remains unchanged.
}
