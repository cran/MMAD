#' MM algorithm based on the AD method for type I multivariate zero-inflated generalized poisson distribution
#'
#' @description Let \eqn{Z \sim Bernoulli(1- \phi)},
#' \eqn{{x} = (X_1,\cdots, X_m)^T}, \eqn{X_i \sim GP(\lambda_i, \theta_i)}, for \eqn{i = 1,\cdots,m} , and \eqn{(Z,X_1,\cdots, X_m)}
#' are mutually independent. An \eqn{m} dimensional discrete random vector \eqn{{y}=(Y_1,\cdots, Y_m)^T} is said to have a Type I
#' multivariate zero-inflated generalized Poisson distribution(ZIGP) distribution if
#'
#' \deqn{ {y} \overset{\rm{d}}= Z{x}= \left\{ \begin{array}{ll}
#' {0} & \text{with probability} \ \lambda \\ {x} & \text{with probability} \ 1-\lambda \\ \end{array} \right. }
#'
#' where \eqn{\lambda \in [0,1)}, \eqn{{\lambda}=(\lambda_1,\cdots, \lambda_m)^T }, \eqn{{\theta}=(\theta_1,\cdots, \theta_m)^T },
#' \eqn{max(-1, -\lambda_i/q_i) < \theta_i \leqslant 1} and \eqn{q_i \geqslant 4} is the largest positive integer for each \eqn{\lambda_i + \theta_{i}q_{i}>0} when \eqn{\theta_{i}<0}.
#' The \code{ZIGPMM} function is used to calculate the Type I multivariate ZIGP model.
#'
#' @param data Data.frame or Matrix that contains corresponding covariates.
#' @param phi0 Probability value for the zero-inflated parameter for ZIGP model.
#' @param la The scale parameter for Generalized Poisson distribution model.
#' @param th The discrete parameter for Generalized Poisson distribution model.
#' @param Maxiter The maximum number of iterations is specified by default as 2000.
#' @param convergence Specify the convergence criterion, the default is 1e-6.
#' @param ... Additional arguments


#'@details The \code{ZIGPMM} function is used to calculate Type I multivariate zero-inflated generalized Poisson distribution model using MM algorithms
#'based on AD technology. \code{data} is provided by user by default, it can be a data frame or a matrix. In addition, the unknown parameters require users to give appropriate initial values,
#'where \code{0<=phi0<1}, each \code{th} should satisfy \code{0<=th<1}, and each \code{la} should be greater than 0.

#' @return An object of class \code{ZIGPMM} that contains the following fields: total amount of observations,
#' the number of iterations, convergence rate, the log likelihood value, estimated results for the unknown parameters,
#' the standard deviation of estimate for the unknown parameters, the likelihood-based 95\% confidence interval for the unknown parameters,
#' information criterion: AIC value and BIC value.
#'
#'
#' @export
#' @references Tian G.L., Huang X.F. and Xu, J.(2019). 'An assembly and decomposition approach for constructing separable minorizing functions in a class of MM algorithms.' \emph{Statistica Sinica} \strong{29}(2), 961-982.
#' @references Huang X.F., Tian G.L., Zhang, C. and Jiang, X.(2017). 'Type I multivariate zero-inflated generalized Poisson distribution with applications.' \emph{Statistics and its Interface} \strong{10}(2), 291-311.
#'
#' @examples
#' x1 <- c(0, 0, 0,38, 0,19,25, 0,25, 0)
#' x2 <- c(0, 0, 0,23, 0,51,24, 0,10, 0)
#' y <- cbind(x1, x2)
#' phi0 = 0.5; la = rep(1,2); th = rep(0.1,2)
#' ZIGPMM(y, phi0, la, th)
#'
ZIGPMM <- function(data, phi0, la, th, Maxiter = 2000, convergence = 1e-06,
                   ...) {
  if (missing(data) || (!is.data.frame(data) & !is.array(data))) {
    stop("'data' is missing or not a data.frame or a matrix or an array")
  }

  if (any(phi0 < 0) | any(phi0 > 1)) {
    stop("The probability `phi0` of Bernoulli distribution should be between 0 and 1",
         call. = FALSE)
  }

  if (length(phi0) == 0 | length(la) == 0 | length(th) == 0) {
    stop("The length of the variable `phi0` or `la` or `th` is not allowed to be 0",
         call. = FALSE)
  }

  if (length(la) != length(th)) {
    stop("The length of the variable `la` and `th` must be the same",
         call. = FALSE)
  }

  if (is.data.frame(data)) {
    y <- as.matrix(data)
  } else {
    y <- data
  }

  n <- nrow(y)
  m <- length(la)
  zero <- matrix(0, n, m)
  yz <- apply(1 * (y == zero), 1, prod)
  n0 <- sum(yz)
  Iy <- 1 - yz

  # log-likelihood function
  a <- n0 * log(phi0 + (1 - phi0) * prod(exp(-la))) + (n - n0) * log(1 -
                                                                       phi0)
  b1 <- sum(Iy * (y == 0) * la)
  b2 <- sum(Iy * (y != 0) * (log(la) + (y - 1) * log(la + th * y) - la -
                               th * y))
  b3 <- sum(Iy * (y != 0) * log(factorial(y)))
  log_ell <- a + b1 + b2 - b3
  el <- c(log_ell)

  error <- 3
  result <- list()

  for (k in 1:Maxiter) {
    if (error > convergence) {
      be0 <- phi0 + (1 - phi0) * prod(exp(-la))
      phi0 <- n0 * phi0/(n * be0)
      for (i in 1:m) {
        ala <- sum(Iy * (y[, i] != 0)) + sum(Iy * (y[, i] != 0) *
                                               la[i] * (y[, i] - 1)/(la[i] + th[i] * y[, i]))
        bla <- n * (1 - phi0)
        la[i] <- ala/bla
        at <- sum(Iy * (y[, i] != 0) * th[i] * y[, i] * (y[, i] -
                                                           1)/(la[i] + th[i] * y[, i]))
        bt <- sum(Iy * (y[, i] != 0) * y[, i])
        th[i] <- at/bt
      }

      a <- n0 * log(phi0 + (1 - phi0) * prod(exp(-la))) + (n - n0) *
        log(1 - phi0)
      b1 <- sum(Iy * (y == 0) * la)
      b2 <- sum(Iy * (y != 0) * (log(la) + (y - 1) * log(la + th *
                                                           y) - la - th * y))
      b3 <- sum(Iy * (y != 0) * log(factorial(y)))
      log_el <- a + b1 + b2 - b3
      el <- append(el, log_el)
      error <- abs(el[k + 1] - el[k])/(1 + abs(el[k]))

      print_k <- k
    }
  }
  if (error > convergence) {
    stop("Convergence result: Did not converge, try adjust 'Maxiter' or 'convergence'.")
  }

  Fisher_Matrix <- -n * ZIGPFx(phi0, la, th)
  std_err <- sqrt(diag(solve(Fisher_Matrix)))
  std_phi0 <- std_err[1]
  std_la <- std_err[2:(1 + m)]
  std_th <- std_err[(2 + m):(1 + 2 * m)]

  # confidence intervals
  ci_phi0_lower <- phi0 - 1.96 * std_phi0
  ci_phi0_upper <- phi0 + 1.96 * std_phi0

  ci_la_lower <- la - 1.96 * std_la
  ci_la_upper <- la + 1.96 * std_la

  ci_th_lower <- th - 1.96 * std_th
  ci_th_upper <- th + 1.96 * std_th


  ELL <- el[length(el)]
  alpha <- c(phi0, la, th)
  Rate <- ZIGP_Rate(phi0, la, th)

  # add values of AIC and BIC
  aic <- (2 * length(alpha)) - (2 * ELL)
  bic <- log(length(y)) * length(alpha) - 2 * ELL
  info_criteria <- c(AIC = aic, BIC = bic)

  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$ELL <- ELL
  result$phi0 <- phi0
  result$std_phi0 <- std_phi0
  result$ci_phi0_lower <- ci_phi0_lower
  result$ci_phi0_upper <- ci_phi0_upper

  result$la <- la
  result$std_la <- std_la
  result$ci_la_lower <- ci_la_lower
  result$ci_la_upper <- ci_la_upper

  result$th <- th
  result$std_th <- std_th
  result$ci_th_lower <- ci_th_lower
  result$ci_th_upper <- ci_th_upper

  result$Rate <- Rate
  result$info_criteria <- info_criteria

  class(result) <- "ZIGP"
  return(invisible(result))

}
