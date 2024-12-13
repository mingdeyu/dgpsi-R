% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{set_id}
\alias{set_id}
\title{Set Emulator ID}
\usage{
set_id(object, id)
}
\arguments{
\item{object}{an emulator object to which the ID will be assigned.}

\item{id}{a unique identifier for the emulator as either a numeric or character
string. Ensure this ID does not conflict with other emulator IDs, especially
when used in linked emulations.}
}
\value{
The updated \code{object}, with the assigned ID stored in its \code{id} slot.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#new}{\figure{lifecycle-new.svg}{options: alt='[New]'}}}{\strong{[New]}}

This function assigns a unique identifier to an emulator.
}
\details{
See further examples and tutorials at \url{https://mingdeyu.github.io/dgpsi-R/dev/}.
}
\examples{
\dontrun{

# See lgp() for an example.
}
}
