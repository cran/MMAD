#' Summary of parameter estimates of a LTN model
#'
#' @description This function returns the result of the \code{LTNMM} function
#'
#'
#' @aliases summary.LTN
#' @usage \method{summary}{LTN}(object, digits = 4, ...)
#' @param object Output from a call to LTN.
#' @param digits The desired number of digits after the decimal point. Default of 4 digits is used.
#' @param ... Additional arguments
#'
#' @return Summary for \code{LTNMM} objects.
#' @seealso \code{\link{LTNMM}}
#' @keywords methods
#' @method summary LTN
#' @export
#'
#' @examples
#'
#'
#' y=c(8.7, 5.4, 8.9, 5.8, 6.2, 9.9, 7.5, 9.5, 6.5, 6.3); a=5
#' result <- LTNMM(y~1, a=5)
#'
#' summary(result,digits=4)
#'
#'
#'
summary.LTN <- function(object, digits = 4, ...) {
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n"))
  cat("\n")
  cat("Model: ")
  cat("Left-Truncated Normal Distribution")
  cat("\n\n")
  LTN_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(LTN_cnt)
  cat("\n")
  converge_cnt <- unlist(object$print_k)
  cat("Number of iterations: ")
  cat(converge_cnt)
  cat("\n")

  print_ell <- round(object$ELL, digits)
  print_rate <- round(object$Rate, digits)
  cat("Convergence Rate: ")
  cat(print_rate)
  cat("\n\n")

  mu <- c(round(object$mu, digits), round(object$std_mu, digits), round(object$ci_mu_lower,
                                                                        digits), round(object$ci_mu_upper, digits))
  sigma <- c(round(object$sigma, digits), round(object$std_sigma, digits),
             round(object$ci_sigma_lower, digits), round(object$ci_sigma_upper,
                                                         digits))
  coef_ltn <- rbind(mu, sigma)
  colnames(coef_ltn) <- c("Estimate", "Std. Error", "Lower 95%-level",
                          "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef_ltn)
  cat("\n")

  cat("Log Likelihood: ")
  cat(as.character(print_ell))
  cat("\n")
  cat("Information Criterion: ")
  cat(paste0("AIC=", round(object$info_criteria[1], digits), " BIC=",
             round(object$info_criteria[2], digits)))
  cat("\n")
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm")
  cat("\n\n")
}
