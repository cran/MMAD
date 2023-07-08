#' Plot the Cox object
#'
#' @param x The Cox object, see \code{\link{CoxMM}}.
#' @param xlab x label, default is 'Time'.
#' @param ylab y label, default is 'Cumulative hazard'.
#' @param type type value, default is 's'.
#' @param lty lty value for line, default is 1.
#' @param lwd line width, default is 1.
#' @param col color parameter, default is gray(0).
#' @param digits The digits after the decimal point, default = 4.
#' @param ... Additional arguments
#'
#' @return the dataframe of 'Time' and accumulative hazard \eqn{\Lambda}.
#' @method plot Cox
#' @export
#'
#' @examples
#' library(survival)
#' result <- CoxMM(Surv(time, status) ~ age + sex, lung)
#'
#' plot(result)
#'
#'
plot.Cox <- function(x, xlab = "Time", ylab = "Cumulative hazard", type = "s",
                     lty = 1, lwd = 1, col = gray(0), digits = 4, ...) {
  LA <- sort(x$Lambda)
  time <- sort(x$time)

  tab <- data.frame(Time = round(time, digits), Lambda = round(LA, digits))
  tab <- unique(tab)
  rownames(tab) <- 1:nrow(tab)
  print(tab)

  # Add a figure frame The meaning of type='n' is not to add any
  # elements to the graph, but only to draw the coordinate axis
  XLIM <- range(time[time != Inf])
  YLIM <- range(LA[LA != Inf])
  plot(XLIM, YLIM, type = "n", xlab = xlab, ylab = ylab, ...)

  do.call("lines", c(list(x = time, y = LA), type = type, lty = lty,
                     lwd = lwd, col = col, ...))

}
