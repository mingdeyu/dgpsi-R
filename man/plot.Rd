% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot}
\alias{plot}
\alias{plot.dgp}
\alias{plot.lgp}
\alias{plot.gp}
\title{Validation plots of a constructed GP, DGP, or linked (D)GP emulator}
\usage{
\method{plot}{dgp}(
  x,
  x_test = NULL,
  y_test = NULL,
  dim = NULL,
  method = "mean_var",
  style = 1,
  min_max = TRUE,
  color = "turbo",
  type = "points",
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)

\method{plot}{lgp}(
  x,
  x_test = NULL,
  y_test = NULL,
  dim = NULL,
  method = "mean_var",
  style = 1,
  min_max = TRUE,
  color = "turbo",
  type = "points",
  M = 50,
  verb = TRUE,
  force = FALSE,
  cores = 1,
  ...
)

\method{plot}{gp}(
  x,
  x_test = NULL,
  y_test = NULL,
  dim = NULL,
  method = "mean_var",
  style = 1,
  min_max = TRUE,
  color = "turbo",
  type = "points",
  verb = TRUE,
  M = 50,
  force = FALSE,
  cores = 1,
  ...
)
}
\arguments{
\item{x}{can be one of the following emulator classes:
\itemize{
\item the S3 class \code{gp}.
\item the S3 class \code{dgp}.
\item the S3 class \code{lgp}.
}}

\item{x_test}{same as that of \code{\link[=validate]{validate()}}.}

\item{y_test}{same as that of \code{\link[=validate]{validate()}}.}

\item{dim}{if \code{dim = NULL}, the index of an emulator's input will be shown on the x-axis in validation plots. Otherwise, \code{dim} indicates
which dimension of an emulator's input will be shown on the x-axis in validation plots:
\itemize{
\item If \code{x} is an instance of the \code{gp} of \code{dgp} class, \code{dim} is an integer.
\item If \code{x} is an instance of the \code{lgp} class, \code{dim} can be
\enumerate{
\item an integer referring to the dimension of the global input to emulators in the first layer of a linked emulator system; or
\item a vector of three integers referring to the dimension (specified by the third integer) of the global input to an emulator
(specified by the second integer) in a layer (specified by the first integer) that is not the first layer of a linked emulator
system.
}
}

This argument is only used when \code{style = 1} and the emulator input is at least two-dimensional. Defaults to \code{NULL}.}

\item{method}{same as that of \code{\link[=validate]{validate()}}.}

\item{style}{either \code{1} or \code{2}, indicating two different types of validation plots.}

\item{min_max}{a bool indicating if min-max normalization will be used to scale the testing output, RMSE, predictive mean and std from the
emulator. Defaults to \code{TRUE}.}

\item{color}{a character string indicating the color map to use when \code{style = 2}:
\itemize{
\item \code{'magma'} (or \code{'A'})
\item \code{'inferno'} (or \code{'B'})
\item \code{'plasma'} (or '\code{C}')
\item \code{'viridis'} (or \code{'D'})
\item \code{'cividis'} (or \code{'E'})
\item \code{'rocket'} (or \code{'F'})
\item \code{'mako'} (or \code{'G'})
\item \code{'turbo'} (or \code{'H'})
}

Defaults to \code{'turbo'} (or \code{'H'}).}

\item{type}{either \code{'line'} or \verb{'points}, indicating whether to draw testing data in the OOS validation plot as a line or
individual points when the input of the emulator is one-dimensional and \code{style = 1}. Defaults to \code{'points'}}

\item{verb}{a bool indicating if the trace information on plotting will be printed during the function execution.
Defaults to \code{TRUE}.}

\item{M}{same as that of \code{\link[=validate]{validate()}}.}

\item{force}{same as that of \code{\link[=validate]{validate()}}.}

\item{cores}{same as that of \code{\link[=validate]{validate()}}.}

\item{...}{N/A.}
}
\value{
A \code{patchwork} object.
}
\description{
This function draws validation plots of a GP, DGP, or linked (D)GP emulator.
}
\details{
See further examples and tutorials at \url{https://mingdeyu.github.io/dgpsi-R/}.
}
\note{
\itemize{
\item \code{\link[=plot]{plot()}} calls \code{\link[=validate]{validate()}} internally to obtain validation results for plotting. However, \code{\link[=plot]{plot()}} will not export the
emulator object with validation results. Instead, it only returns the plotting object. For small-scale validations (i.e., small
training or testing data points), direct execution of \code{\link[=plot]{plot()}} is fine. However, for moderate- to large-scale validations,
it is recommended to first run \code{\link[=validate]{validate()}} to obtain and store validation results in the emulator object, and then supply the
object to \code{\link[=plot]{plot()}}. This is because if an emulator object has the validation results stored, each time when \code{\link[=plot]{plot()}}
is invoked, unnecessary evaluations of repetitive LOO or OOS validation will not be implemented.
\item \code{\link[=plot]{plot()}} uses information provided in \code{x_test} and \code{y_test} to produce the OOS validation plots. Therefore, if validation results
are already stored in \code{x}, unless \code{x_test} and \code{y_test} are identical to those used by \code{\link[=validate]{validate()}}, \code{\link[=plot]{plot()}} will re-evaluate OOS
validations before plotting.
\item Any R vector detected in \code{x_test} and \code{y_test} will be treated as a column vector and automatically converted into a single-column
R matrix. Thus, if \code{x_test} or \code{y_test} is a single testing data point with multiple dimensions, it must be given as a matrix.
\item The returned \code{patchwork} object contains the \code{ggplot2} objects. One can modify the included individual ggplots
by accessing them with double-bracket indexing. See \url{https://patchwork.data-imaginist.com/} for further information.
}
}
\examples{
\dontrun{

# See gp(), dgp(), or lgp() for an example.
}
}
