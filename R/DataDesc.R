#' Kidney Infection Data
#'
#' @description The data consisted of the time to first and second infection relapse in 38 kidney disease patients using a portable dialysis machine.
#' Infection may occur where the catheter was inserted. Catheters are subsequently removed if infection develops and may be removed for other reasons,
#' in which case observations are censored.
#'
#' @return Kidney infection data contains the following fields:
#' \item{patient}{id}
#' \item{time}{time}
#' \item{status}{event status}
#' \item{age}{in years}
#' \item{sex}{1=male, 2=female}
#' \item{disease}{disease type (0=GN, 1=AN, 2=PKD, 3=Other)}
#' \item{frail}{frailty estimate from original paper}
#'
#'
#' @references McGilchrist C.A. and Aisbett C.W.(1991). "Regression with frailty in survival analysis." \emph{Biometrics} \strong{47}, 461-466.
#'
#' @examples
#' data = data(Kidney)

"kidney"

#' Breast Cosmesis Data
#'
#' @description The often used data set for interval censored data, described and given in full in Finkelstein and Wolfe (1985).
#'
#' @return Breast cosmesis data contains the following fields:
#' \item{left}{a numeric vector}
#' \item{right}{a numeric vector}
#' \item{treatment}{a factor with levels Rad and RadChem}
#'
#'
#' @references Finkelstein D.M. and Wolfe R.A.(1985). "A semiparametric model for regression analysis of interval-censored failure time data." \emph{Biometrics} \strong{41}, 933-945.
#'
#' @examples
#' data = data(bcos)

"bcos"

#' NCCTG Lung Cancer Data
#'
#' @description Survival in patients with advanced lung cancer from the North Central Cancer Treatment Group. Performance scores rate how well the patient can perform usual daily activities.
#'
#' @return Kidney infection data contains the following fields:
#' \item{inst}{Institution code}
#' \item{time}{Survival time in days}
#' \item{status}{censoring status 1=censored, 2=dead}
#' \item{age}{Age in years}
#' \item{sex}{Male=1 Female=2}
#' \item{ph.ecog}{ECOG performance score as rated by the physician. 0=asymptomatic, 1= symptomatic but completely ambulatory}
#' \item{ph.karno}{Karnofsky performance score (bad=0-good=100) rated by physician}
#' \item{pat.karno}{Karnofsky performance score as rated by patient}
#' \item{meal.cal}{Calories consumed at meals}
#' \item{wt.loss}{Weight loss in last six months (pounds)}
#'
#'
#' @references Finkelstein D.M. and Wolfe R.A.(1985). "A semiparametric model for regression analysis of interval-censored failure time data." \emph{Biometrics} \strong{41}, 933-945.
#'
#' @examples
#' data = data(lung)

"lung"

#' The children’s absenteeism data in Indonesia
#'
#' @description In a survey of Indonesian family life conducted by Strauss et al. the participants included 7000 households sampled from 321 communities randomly selected from 13 of the nation’s 26 Provinces,
#' in which 83\% of the Indonesian population lived. Among those households with one child per household, 437 household heads were asked questions about the health of their children.
#'
#' @return The children’s absenteeism data in Indonesia contains the following fields:
#' \item{y1}{The number of days the children missed their primary activities due to illness in the last four weeks}
#' \item{y2}{The number of days the children spent in bed due to illness in the last four weeks}
#'
#'
#' @references Huang X.F., Tian G.L., Zhang, C. and Jiang, X.(2017). "Type I multivariate zero-inflated generalized Poisson distribution with applications." \emph{Statistics and its Interface} \strong{10}(2), 291-311.
#' @references Strauss J., Beegle K., Sikoki B., Dawiyanto A., Herawati Y. and Witoelar Y.(2004). "The Third Wave of the Indonesia Family Life Survey (IFLS): Overview and Field Report, WR-144/1-NIA/NICHD, RAND Corporation, Santa Monica, CA."
#'
#' @examples
#' data = data(cadi)

"cadi"

#' Voluntary and involuntary job changes data
#'
#' @description Jung and Winkelmann(1993) provided data on both the numbers of voluntary and involuntary job changes of males during ten period 1974–1984.
#' The samples contain 2124 males who started their working career before or in 1974 and did not retire before 1984.
#'
#' @return Voluntary and involuntary job changes data contains the following fields:
#' \item{y1}{Job changes after experiencing an unemployment spell(assumed to be involuntary)}
#' \item{y2}{Direct job to job changes(which are assumed to be voluntary) }
#'
#'
#' @references Huang X.F., Tian G.L., Zhang, C. and Jiang, X.(2017). "Type I multivariate zero-inflated generalized Poisson distribution with applications." \emph{Statistics and its Interface} \strong{10}(2), 291-311.
#' @references Jung R.C. and Winkelmann R.(1993). "Two aspects of labor mobility: A bivariate Poisson regression approach." \emph{Empirical Economics} \strong{18}(3), 543–556.
#'
#' @examples
#' data = data(vijc)

"vijc"
