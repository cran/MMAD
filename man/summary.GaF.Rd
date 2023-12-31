% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.GaF.R
\name{summary.GaF}
\alias{summary.GaF}
\title{Summary of parameter estimates of a gamma frailty model}
\usage{
\method{summary}{GaF}(object, digits = 4, ...)
}
\arguments{
\item{object}{Output from a call to GaF.}

\item{digits}{The desired number of digits after the decimal point. Default of 4 digits is used.}

\item{...}{Additional arguments}
}
\value{
Summary for \code{GaFrailtyMM} objects.
}
\description{
This function returns the result of the \code{GaFrailtyMM} function
}
\examples{


library(survival)
result <- GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), data=kidney)

summary(result,digits=4)



}
\seealso{
\code{\link{GaFrailtyMM}}
}
\keyword{methods}
