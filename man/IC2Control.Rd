% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IC2Control.R
\name{IC2Control}
\alias{IC2Control}
\title{Control IC2Pro object}
\usage{
IC2Control(Maxiter = 2000, convergence = 1e-06, Idigits = 4, Pdigits = 4)
}
\arguments{
\item{Maxiter}{The maximum number of iterations is specified by default as 2000.}

\item{convergence}{Specify the convergence criterion, the default is 1e-6.}

\item{Idigits}{The number of decimal places for the survival interval values.}

\item{Pdigits}{The number of decimal places for the survival probability values.}
}
\value{
list of Maxiter, convergence, Idigits, Pdigits.
}
\description{
Control IC2Pro object
}
\examples{
IC2Control()
}
