% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alm.R
\name{alm}
\alias{alm}
\alias{alm.gp}
\alias{alm.dgp}
\alias{alm.bundle}
\title{Locate the next design point for a (D)GP emulator or a bundle of (D)GP emulators using ALM}
\usage{
alm(object, x_cand, ...)

\method{alm}{gp}(object, x_cand, batch_size = 1, M = 50, workers = 1, ...)

\method{alm}{dgp}(object, x_cand, batch_size = 1, M = 50, workers = 1, aggregate = NULL, ...)

\method{alm}{bundle}(object, x_cand, batch_size = 1, M = 50, workers = 1, aggregate = NULL, ...)
}
\arguments{
\item{object}{can be one of the following:
\itemize{
\item the S3 class \code{gp}.
\item the S3 class \code{dgp}.
\item the S3 class \code{bundle}.
}}

\item{x_cand}{a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
from which the next design point(s) are determined. If \code{object} is an instance of the \code{bundle} class, \code{x_cand} could also
be a list with the length equal to the number of emulators contained in the \code{object}. Each slot in \code{x_cand} is a matrix
that gives a candidate set for each emulator included in the bundle. See \emph{Note} section below for further information.}

\item{...}{any arguments (with names different from those of arguments used in \code{\link[=alm]{alm()}}) that are used by \code{aggregate}
can be passed here.}

\item{batch_size}{an integer that gives the number of design points to be chosen.
Defaults to \code{1}.}

\item{M}{the size of the conditioning set for the Vecchia approximation in the criterion calculation. This argument is only used if the emulator \code{object}
was constructed under the Vecchia approximation. Defaults to \code{50}.}

\item{workers}{the number of processes to be used for the criterion calculation. If set to \code{NULL},
the number of processes is set to \verb{max physical cores available \%/\% 2}. Defaults to \code{1}.}

\item{aggregate}{an R function that aggregates scores of the ALM across different output dimensions (if \code{object} is an instance
of the \code{dgp} class) or across different emulators (if \code{object} is an instance of the \code{bundle} class). The function should be specified in the
following basic form:
\itemize{
\item the first argument is a matrix representing scores. The rows of the matrix correspond to different design points. The number of columns
of the matrix equals to:
\itemize{
\item the emulator output dimension if \code{object} is an instance of the \code{dgp} class; or
\item the number of emulators contained in \code{object} if \code{object} is an instance of the \code{bundle} class.
}
\item the output should be a vector that gives aggregations of scores at different design points.
}

Set to \code{NULL} to disable the aggregation. Defaults to \code{NULL}.}
}
\value{
\itemize{
\item If \code{object} is an instance of the \code{gp} class, a vector is returned with the length equal to \code{batch_size}, giving the positions (i.e., row numbers)
of next design points from \code{x_cand}.
\item If \code{object} is an instance of the \code{dgp} class, a matrix is returned with row number equal to \code{batch_size} and column number equal to one (if \code{aggregate}
is not \code{NULL}) or the output dimension (if \code{aggregate} is \code{NULL}), giving positions (i.e., row numbers) of next design points from \code{x_cand} to be added
to the DGP emulator across different outputs. If \code{object} is a DGP emulator with either \code{Hetero} or \code{NegBin} likelihood layer, the returned matrix has
two columns with the first column giving positions of next design points from \code{x_cand} that correspond to the mean parameter of the normal or negative Binomial
distribution, and the second column giving positions of next design points from \code{x_cand} that correspond to the variance parameter of the normal distribution or
the dispersion parameter of the negative Binomial distribution.
\item If \code{object} is an instance of the \code{bundle} class, a matrix is returned with row number equal to \code{batch_size} and column number equal to the number of
emulators in the bundle, giving positions (i.e., row numbers) of next design points from \code{x_cand} to be added to individual emulators.
}
}
\description{
This function searches from a candidate set to locate the next design point(s) to be added to a (D)GP emulator
or a bundle of (D)GP emulators using the Active Learning MacKay (ALM), see the reference below.
}
\details{
See further examples and tutorials at \url{https://mingdeyu.github.io/dgpsi-R/}.
}
\note{
\itemize{
\item The column order of the first argument of \code{aggregate} must be consistent with the order of emulator output dimensions (if \code{object} is an instance of the
\code{dgp} class), or the order of emulators placed in \code{object} if \code{object} is an instance of the \code{bundle} class;
\item If \code{x_cand} is supplied as a list when \code{object} is an instance of \code{bundle} class and a \code{aggregate} function is provided, the matrices in \code{x_cand} must have
common rows (i.e., the candidate sets of emulators in the bundle have common input locations) so the \code{aggregate} function can be applied.
\item Any R vector detected in \code{x_cand} will be treated as a column vector and automatically converted into a single-column
R matrix.
}
}
\examples{
\dontrun{

# load packages and the Python env
library(lhs)
library(dgpsi)

# construct a 1D non-stationary function
f <- function(x) {
 sin(30*((2*x-1)/2-0.4)^5)*cos(20*((2*x-1)/2-0.4))
}

# generate the initial design
X <- maximinLHS(10,1)
Y <- f(X)

# training a 2-layered DGP emulator with the global connection off
m <- dgp(X, Y, connect = F)

# generate a candidate set
x_cand <- maximinLHS(200,1)

# locate the next design point using ALM
next_point <- alm(m, x_cand = x_cand)
X_new <- x_cand[next_point,,drop = F]

# obtain the corresponding output at the located design point
Y_new <- f(X_new)

# combine the new input-output pair to the existing data
X <- rbind(X, X_new)
Y <- rbind(Y, Y_new)

# update the DGP emulator with the new input and output data and refit
m <- update(m, X, Y, refit = TRUE)

# plot the LOO validation
plot(m)
}
}
\references{
MacKay, D. J. (1992). Information-based objective functions for active data selection. \emph{Neural Computation}, \strong{4(4)}, 590-604.
}
