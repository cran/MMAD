% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataDesc.R
\docType{data}
\name{lung}
\alias{lung}
\title{NCCTG Lung Cancer Data}
\format{
An object of class \code{data.frame} with 228 rows and 10 columns.
}
\usage{
lung
}
\value{
Kidney infection data contains the following fields:
\item{inst}{Institution code}
\item{time}{Survival time in days}
\item{status}{censoring status 1=censored, 2=dead}
\item{age}{Age in years}
\item{sex}{Male=1 Female=2}
\item{ph.ecog}{ECOG performance score as rated by the physician. 0=asymptomatic, 1= symptomatic but completely ambulatory}
\item{ph.karno}{Karnofsky performance score (bad=0-good=100) rated by physician}
\item{pat.karno}{Karnofsky performance score as rated by patient}
\item{meal.cal}{Calories consumed at meals}
\item{wt.loss}{Weight loss in last six months (pounds)}
}
\description{
Survival in patients with advanced lung cancer from the North Central Cancer Treatment Group. Performance scores rate how well the patient can perform usual daily activities.
}
\examples{
data = data(lung)
}
\references{
Finkelstein D.M. and Wolfe R.A.(1985). "A semiparametric model for regression analysis of interval-censored failure time data." \emph{Biometrics} \strong{41}, 933-945.
}
\keyword{datasets}
