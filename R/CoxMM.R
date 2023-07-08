#' MM algorithm based on AD technology for Cox model
#'
#' @description Let \eqn{T_i, C_i} and \eqn{X_i = (x_{i1},\cdots, x_{iq})^T} denote the,
#' survival time, the censoring time and a \eqn{q} dimension vector of coefficients for the \eqn{i}-th individual, respectively. And assume the censoring time
#' \eqn{C_i} is independent of the survival time \eqn{T_i} are mutually independent, and \eqn{I_i = I(T_{i} \leqslant C_{i})} is the censoring indicator.
#' Then the instantaneous hazard rate function of \eqn{T_i} is
#'
#' \deqn{\lambda(t|X_i)=\lambda_{0}(t) \exp(X_{i}^{T} \beta)}
#'
#' where \eqn{\lambda_{0}(.)} is a baseline hazard rate and \eqn{\beta = (\beta_1, \cdots, \beta_q)^{T}} is a vector of regression parameters.
#' We denote \eqn{\Lambda} as the accumulative hazard rate. Then the observed data likelihood function is
#'
#' \deqn{ L(\alpha | Y_{obs}) = \prod_{i=1}^n (\lambda_{0}(t_i) \exp(X_{i}^{T} \beta))^{I_i} \exp(-\Lambda(t_i) \exp(X_{i}^{T} \beta)) }
#'
#' where \eqn{\alpha = (\beta, \Lambda)}. The \code{CoxMM} function is used to calculate the Cox model.
#'
#' @param formula A formula object, which contains on the left hand side an object of the type \code{Surv}
#' and on the right hand side is the terms, e.g. \code{formula=Surv(time, status) ~ x}.
#' @param data A \code{data.frame} in which to interpret the variables named in the formula.
#' @param beta A vector of unknown regression parameters, default is \code{NULL}.
#' If is \code{NULL}, then make all \code{beta=0} during calculation.
#' @param Maxiter The maximum number of iterations is specified by default as 2000.
#' @param convergence Specify the convergence criterion, the default is 1e-6.
#' @param ... Additional arguments

#'
#' @details The \code{CoxMM} function is used to calculate the Cox model using MM algorithms
#' based on AD technology. EM algorithms rely on the fact that, after profiling out the nonparametric component \eqn{\Lambda},
#' the resulting function is concave. However, when this assumption does not hold, maximizing the resulting function using Newton’s method becomes difficult,
#' especially when there are a large number of covariates. MM algorithms can avoid the
#' concavity requirement and bypass the need for Newton method and matrix inversion.
#'
#' @return An object of class \code{CoxMM} that contains the following fields: the Time, total amount of observations,
#' total number of failure events, the variable name, the \eqn{\beta}, the \eqn{\lambda}, the \eqn{\Lambda}, convergence result,
#' the log likelihood value, the standard deviation of the estimated \eqn{\beta}, the likelihood-based 95\% confidence interval for the \eqn{\beta}.
#'
#' @importFrom grDevices gray
#' @importFrom survival is.Surv Surv strata
#' @importFrom stats coef dnorm lm model.extract model.frame model.matrix pnorm terms
#'
#' @export
#'
#' @references D.R. Cox.(1972). 'Regression models and life tables.' \emph{Journal of the Royal Statistical Society(Series B)} \strong{34}(2), 187-220.
#' @references Zhang L.L. and Huang X.F.(2022). 'On MM algorithms for Cox model with right-censored data.' \emph{In International Conference on Cloud Computing, Internet of Things, and Computer Applications (CICA 2022)} \strong{12303}, 29-38.
#'
#' @examples
#' library(survival)
#' CoxMM(Surv(time, status) ~ age + sex, lung)
#'
CoxMM <- function(formula, data, beta = NULL, Maxiter = 2000, convergence = 1e-06,
                  ...) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' is missing or incorrect")
  }
  if (missing(data) || !inherits(data, "data.frame")) {
    stop("'data' is missing or not an object of type data.frame")
  }

  mf <- model.frame(formula, data)
  mx <- model.matrix(formula, data)

  X <- mx[, -c(1), drop = FALSE]
  namesX <- colnames(X)
  n <- dim(data)[1]
  q <- length(namesX)
  time <- matrix(mf[[1]][, 1], c(n, 1), byrow = TRUE)
  status <- matrix(mf[[1]][, 2], c(n, 1), byrow = TRUE)

  x <- matrix(X, n, q)
  event <- sum(status)
  vd <- as.vector(status)
  vy <- as.vector(time)
  if (is.null(beta) == TRUE) {
    be <- c(rep(0, q))
  } else {
    be <- beta
  }

  CoxLogLik <- function(lambda, x, d, vy, vd, beta) {
    La <- (cumsum(lambda[order(vy)]))[rank(vy)]
    ell <- sum(log(lambda[vd != 0])) + sum(d * (x %*% beta)) - sum(La * exp(x %*% beta))
    return(ell)
  }

  la <- rep(1/n, n)
  error <- 1
  log_ell <- CoxLogLik(la, x, status, vy, vd, be)
  ell <- c(log_ell)

  for (k in 1:Maxiter) {
    if (error > convergence) {
      re <- CoxProfile(x, status, time, be, la, n, q)

      be <- re$beta
      la <- re$lambda
      log_el <- CoxLogLik(la, x, status, vy, vd, be)
      ell <- append(ell, log_el)
      error <- abs(ell[k + 1] - ell[k])/(1 + abs(ell[k]))
    }

  }
  Lambda <- (cumsum(la[order(vy)]))[rank(vy)]
  if (error > convergence) {
    stop("Convergence result: Did not converge, try adjust 'Maxiter' or 'convergence'.")
  }

  ELL <- ell[length(ell)]
  Hessian_Matrix <- HessianC(n, q, x, status, vd, la, be, vy)
  std_err <- sqrt(diag(solve(Hessian_Matrix)))

  std_be <- std_err[1:q]

  # confidence intervals
  ci_be_lower <- be - 1.96 * std_be
  ci_be_upper <- be + 1.96 * std_be

  result <- list()
  result$call <- match.call()
  result$time <- time
  result$n <- n
  result$event <- event
  result$namesX <- namesX
  result$be <- be
  result$la <- la
  result$Lambda <- Lambda
  result$error <- error
  result$loglik <- ELL
  result$std_be <- std_be
  result$ci_be_lower <- ci_be_lower
  result$ci_be_upper <- ci_be_upper

  class(result) <- "Cox"
  return(invisible(result))
}
