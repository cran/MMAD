#' Summary of parameter estimates of a IC2 model
#'
#' @description This function returns the result of the \code{IC2MM} function
#'
#' @aliases summary.IC2
#' @usage \method{summary}{IC2}(object, ...)
#' @param object Output from a call to IC2.
#' @param ... Additional arguments
#'
#' @return Summary for \code{IC2MM} objects.
#' @seealso \code{\link{IC2MM}}
#' @keywords methods
#' @method summary IC2
#' @export
#'
#' @examples
#'
#'
#' library(survival)
#' result <- IC2MM(Surv(left, right, type = 'interval2') ~ treatment, bcos)
#'
#' summary(result)
#'
#'
#'
summary.IC2 <- function(object, ...) {
  tab <- object$df_tab
  k <- length(object$df_tab[[2]])
  # if more than one strata, then print each strata in order
  if (!is.null(object$strata) && length(object$strata) > 1) {
    cnt <- 1
    for (i in 1:length(object$strata)) {
      cat(paste(names(object$strata)[i], ":", sep = ""))
      endCnt <- cnt + object$strata[i] - 1
      cat("\n")
      tabi <- tab[cnt:endCnt, ]
      dimnames(tabi)[[1]] <- 1:object$strata[i]
      print(tabi)
      cnt <- endCnt + 1
    }
  } else {
    dimnames(tab)[[1]] <- 1:k
    print(tab)
  }
}
