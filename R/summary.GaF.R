#' Summary of parameter estimates of a gamma frailty model
#'
#' @description This function returns the result of the \code{GaFrailtyMM} function
#'
#'
#' @aliases summary.GaF
#' @usage \method{summary}{GaF}(object, digits = 4, ...)
#' @param object Output from a call to GaF.
#' @param digits The desired number of digits after the decimal point. Default of 4 digits is used.
#' @param ... Additional arguments
#'
#' @return Summary for \code{GaFrailtyMM} objects.
#' @seealso \code{\link{GaFrailtyMM}}
#' @keywords methods
#' @method summary GaF
#' @export
#'
#' @examples
#'
#'
#' library(survival)
#' result <- GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), data=kidney)
#'
#' summary(result,digits=4)
#'
#'
#'
summary.GaF <- function(object, digits = 4, ...) {
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n"))
  cat("\n")
  cat("Model: ")
  cat("Gamma Frailty Survival Models")
  cat("\n\n")

  cat(paste0("n= ", object$N, ", nubetar of events=", object$TotalCen))
  cat("\n")

  converge_value <- unlist(object$print_err)
  cat("\n")
  cat("Convergence result: ")
  cat(converge_value)
  cat("\n")

  print_ell <- round(object$ELL, digits)
  cat("\n")

  m <- length(object$be)
  theta <- c(round(object$th, digits), round(object$std_th, digits),
             round(object$ci_th_lower, digits), round(object$ci_th_upper, digits))
  be <- c(round(object$be, digits), round(object$std_be, digits), round(object$ci_be_lower,
                                                                        digits), round(object$ci_be_upper, digits))

  be <- matrix(be, m, 4)

  rownames(be) <- object$namesX

  coef_Frailty <- rbind(theta, be)
  colnames(coef_Frailty) <- c("Estimate", "Std. Error", "Lower 95%-level",
                              "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef_Frailty)
  cat("\n")

  cat("Log Likelihood: ")
  cat(as.character(print_ell))
  cat("\n")
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm")
  cat("\n\n")
}
