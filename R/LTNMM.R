#' MM algorithm based on the AD method for left-truncated normal distribution
#'
#' @description The \code{LTNMM} function is used to calculate a left-truncated normal distribution model. A \eqn{ LTN(\mu, \sigma^2; a) } has the density function
#'
#' \deqn{f(y; \mu, \sigma^2; a) = \frac{1}{c \sqrt{2 \pi \sigma^{2}}} \exp{( -\frac{(y-\mu)^{2}}{2 \sigma^{2}} )} \centerdot I(y \geq a) }
#'
#' where \eqn{(\mu, \sigma^2)} are two unknown parameters, \eqn{ a } is a known constant, \eqn{ c = 1- \Phi(\frac{a-u}{ \sigma}) }, and  \eqn{ \Phi(\centerdot) } is the cdf of the standard normal distribution.
#'
#' @param formula A formula object which symbolically describes the model to calculated.
#' @param a A numeric scalar of the known left truncation value.
#' @param mu The mean of the normal distribution is set to NULL by default. If the distribution is truncated, we use estimates from OLS.
#' @param sigma The variance of the normal distribution is set to NULL by default. If the distribution is truncated, we use estimates from OLS.
#' @param data List that contains corresponding covariates. If none is provided then assumes objects are in user’s environment.
#' @param Maxiter The maximum number of iterations is specified by default as 2000.
#' @param convergence Specify the convergence criterion, the default is 1e-6.
#' @param ... Additional arguments

#'
#' @return An object of class \code{LTNMM} that contains the following fields: total amount of observations,
#' the number of iterations, convergence rate, the log likelihood value, estimated results for the unknown parameters,
#' the standard deviation of estimate for the unknown parameters, the likelihood-based 95\% confidence interval for the unknown parameters,
#' information criterion: AIC value and BIC value.

#' @details The \code{LTNMM} function is used to calculate a left-truncated normal distribution model using MM algorithms based on AD technology.
#' The \code{formula} parameter can be used to provide the data that needs to be calculated, such as \code{formula=y~1}. By default, the
#' \code{data} is provided by the user’s environment. The initial values of the mean and variance of the normal distribution are estimated using OLS.
#'
#'
#' @export
#' @references Tian G.L., Huang X.F., and Xu, J.(2019). 'An assembly and decomposition approach for constructing separable minorizing functions in a class of MM algorithms.' \emph{Statistica Sinica} \strong{29}(2), 961-982.
#'
#' @examples
#' y=c(8.7, 5.4, 8.9, 5.8, 6.2, 9.9, 7.5, 9.5, 6.5, 6.3); a=5
#' LTNMM(y~1, a=5)
#'
LTNMM <- function(formula, a, mu = NULL, sigma = NULL, data = sys.frame(sys.parent()),
                  Maxiter = 2000, convergence = 1e-06, ...) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' is missing or incorrect")
  }

  y <- model.frame(formula, data)[, 1]

  # check 'a'
  if (!is.numeric(a)) {
    stop("`a` must be numeric", call. = FALSE)
  }
  if (length(a) != 1) {
    stop("`a` must be scalars", call. = FALSE)
  }

  if (a == -Inf | a == Inf) {
    stop("`a` must be a finite value", call. = FALSE)
  }
  if (any(y < a)) {
    stop("observed values below specified truncation `a`", call. = FALSE)
  }

  # validate and process initial values of mu and sigma
  mu_null <- is.null(mu)
  sigma_null <- is.null(sigma)
  if (mu_null == TRUE & sigma_null == FALSE) {
    stop("Please give the value of mu", call. = FALSE)
  }
  if (mu_null == FALSE & sigma_null == TRUE) {
    stop("Please give the value of sigma", call. = FALSE)
  }

  if (mu_null == TRUE & sigma_null == TRUE) {
    lm_mod <- lm(y ~ 1)
    mu <- unname(coef(lm_mod))
    sigma <- unname((summary(lm_mod)$sigma)^2)
  }

  error <- 3
  n <- length(y)
  alpha0 <- c(mu, sigma)
  si <- sqrt(sigma)
  result <- list()

  # log-likelihood function
  log_ell <- -n * log(sigma)/2 - sum((y - mu)^2)/(2 * sigma) - n * log(1 -
                                                                         pnorm((a - mu)/si)) - n * log(2 * pi)/2
  el <- c(log_ell)

  for (k in 1:Maxiter) {
    if (error > convergence) {
      a1 <- (a - mu)/si
      w <- 1 - pnorm(a1)
      s1 <- (1 - w)/w
      tao <- exp(-(a - mu)^2/(2 * sigma))/sqrt(2 * pi * sigma)
      g <- tao/pnorm(a1)
      # deta = sigma - sigma*(a-mu)*g

      mu1 <- (mean(y) + s1 * (mu - sigma * g))/(1 + s1)
      deta <- sigma + (mu1 - mu)^2 - sigma * (a + mu - 2 * mu1) *
        g
      mu <- mu1
      sigma <- (sum((y - mu)^2)/n + s1 * deta)/(1 + s1)

      si <- sqrt(sigma)
      log_el <- -n * log(sigma)/2 - sum((y - mu)^2)/(2 * sigma) -
        n * log(1 - pnorm((a - mu)/si)) - n * log(2 * pi)/2
      el <- append(el, log_el)
      error <- abs(el[k + 1] - el[k])/(abs(el[k]) + 1)

      print_k <- k

    }
  }
  if (error > convergence) {
    stop("Convergence result: Did not converge, try adjust 'Maxiter' or 'convergence'.")
  }

  Fisher_Matrix <- -n * LTNFx(a, mu, sigma)
  std_err <- sqrt(diag(solve(Fisher_Matrix)))
  std_mu <- std_err[1]
  std_sigma <- std_err[2]

  # confidence intervals
  ci_mu_lower <- mu - 1.96 * std_mu
  ci_mu_upper <- mu + 1.96 * std_mu
  ci_sigma_lower <- sigma - 1.96 * std_sigma
  ci_sigma_upper <- sigma + 1.96 * std_sigma

  ELL <- el[length(el)]
  alpha <- c(mu, sigma)
  Rate <- LTN_CRate(a, mu, sigma)

  # add values of AIC and BIC
  aic <- (2 * length(alpha)) - (2 * ELL)
  bic <- log(length(y)) * length(alpha) - 2 * ELL
  info_criteria <- c(AIC = aic, BIC = bic)

  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$ELL <- ELL
  result$mu <- mu
  result$std_mu <- std_mu
  result$ci_mu_lower <- ci_mu_lower
  result$ci_mu_upper <- ci_mu_upper
  result$sigma <- sigma
  result$std_sigma <- std_sigma
  result$ci_sigma_lower <- ci_sigma_lower
  result$ci_sigma_upper <- ci_sigma_upper
  result$Rate <- Rate
  result$info_criteria <- info_criteria

  class(result) <- "LTN"
  return(invisible(result))
}
