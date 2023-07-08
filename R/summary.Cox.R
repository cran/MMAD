#' Summary of parameter estimates of a Cox model
#'
#' @description This function returns the result of the \code{CoxMM} function
#'
#'
#' @aliases summary.Cox
#' @usage \method{summary}{Cox}(object, digits = 4, ...)
#' @param object Output from a call to Cox.
#' @param digits The desired number of digits after the decimal point. Default of 4 digits is used.
#' @param ... Additional arguments
#'
#' @return Summary for \code{CoxMM} objects.
#' @seealso \code{\link{CoxMM}}
#' @keywords methods
#' @method summary Cox
#' @export
#'
#' @examples
#'
#' library(survival)
#' result <- CoxMM(Surv(time, status) ~ age + sex, lung)
#'
#' summary(result,digits=4)
#'
#'
#'
summary.Cox <- function(object, digits = 4, ...) {
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n"))
  cat("\n")
  cat("Model: ")
  cat("Cox Model")
  cat("\n\n")

  cat(paste0("n= ", object$n, ", number of events=", object$event))
  cat("\n\n")

  converge_value <- unlist(object$error)
  cat("Convergence result: ")
  cat(converge_value)
  cat("\n")
  loglik <- round(object$loglik, digits)
  cat("\n")

  m <- length(object$be)
  be <- c(round(object$be, digits), round(object$std_be, digits), round(object$ci_be_lower,
                                                                        digits), round(object$ci_be_upper, digits))
  be <- matrix(be, m, 4)
  rownames(be) <- object$namesX

  coef <- be
  colnames(coef) <- c("Estimate", "Std. Error", "Lower 95%-level", "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef)
  cat("\n")

  cat("Log Likelihood: ")
  cat(as.character(loglik))
  cat("\n")
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm")
  cat("\n\n")
}
