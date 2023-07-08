#' Plot the IC2 object
#'
#' @param x The IC2 object, see \code{\link{IC2MM}}.
#' @param xlab x label, default is 'Time'.
#' @param ylab y label, default is 'Survival'.
#' @param legend legend, default=NULL.
#' @param main figure title, default is 'Survival Function'
#' @param lty lty value for line, default is 1:9.
#' @param lwd line width, default is 1.
#' @param xleg positional parameters of the legend, default=0.
#' @param yleg positional parameters of the legend, default=0.15 .
#' @param col the color of the drawing, default=gray(0)
#' @param ... Additional arguments
#'
#' @return A list of arguments for the legend. Values are x, y, legend, fill, lty, bty, col.
#' @method plot IC2
#'
#' @export
#'
#' @examples
#'
#' library(survival)
#' result = IC2MM(Surv(left, right, type = 'interval2') ~ treatment, bcos)
#'
#' plot(result, col=c('red', 'blue'))
#'
#'
plot.IC2 <- function(x, xlab = "Time", ylab = "Survival", legend = NULL,
                     main = "Survival Function", lty = 1:9, lwd = 1, xleg = 0, yleg = 0.15,
                     col = gray(0), ...) {
  estpar <- list(lty = lty, lwd = lwd, col = col)

  # Add a figure frame
  YLIM <- c(0, 1)
  # The meaning of type='n' is not to add any elements to the
  # graph, but only to draw the coordinate axis
  time <- c(0, as.vector(x$s))
  XLIM <- range(time[time != Inf])
  plot(XLIM, YLIM, type = "n", xlab = xlab, ylab = ylab, main = main,
       ...)

  pickpari <- function(parlist, i) {
    picki <- function(x, i) {
      if (length(x) >= i) {
        out <- x[i]
      } else if (length(x) >= 1) {
        out <- x[1]
      } else {
        out <- 1
      }
      out
    }
    outlist <- parlist
    n <- length(parlist)
    for (j in 1:n) {
      outlist[[j]] <- picki(parlist[[j]], i)
    }
    outlist
  }

  lines.ic2 <- function(x, i, parlist = estpar) {
    parlist <- pickpari(parlist, i)

    s <- as.numeric(unlist(x$s))
    S <- as.numeric(unlist(x$S))
    locate0 <- which(s == 0)
    sart0 <- locate0[i]
    if (i == length(locate0)) {
      end0 <- length(s)
    } else {
      end0 <- (locate0[i + 1] - 1)
    }

    s <- s[sart0:end0]
    S <- S[sart0:end0]
    s[s == Inf] <- max(s[length(s) - 1])

    do.call("lines", c(list(x = s, y = S), parlist))
  }

  nstrata <- length(x$strata)
  if (nstrata == 0)
    nstrata <- 1

  if (nstrata > 1) {
    for (i in 1:nstrata) {
      lines.ic2(x, i)
    }
  } else {
    lines.ic2(x, 1)
  }

  xleg <- max(0, xleg)
  yleg <- max(0.15, yleg)

  legend.list <- list(x = xleg, y = yleg, legend = names(x$strata), lty = lty[1:nstrata],
                      bty = "n", col = col)

  if (is.null(legend)) {
    if (nstrata > 1)
      do.call("legend", legend.list)
  } else if (legend)
    do.call("legend", legend.list)

  invisible(legend.list)

}
