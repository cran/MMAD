#' Control IC2Pro object
#'
#' @param Maxiter The maximum number of iterations is specified by default as 2000.
#' @param convergence Specify the convergence criterion, the default is 1e-6.
#' @param Idigits The number of decimal places for the survival interval values.
#' @param Pdigits The number of decimal places for the survival probability values.
#'
#' @return list of Maxiter, convergence, Idigits, Pdigits.
#' @export
#'
#' @examples
#' IC2Control()
IC2Control <- function(Maxiter = 2000, convergence = 1e-06, Idigits = 4,
                       Pdigits = 4) {
  if (!is.numeric(convergence) || convergence <= 0)
    stop("value of 'convergence' must be > 0")
  if (!is.numeric(Maxiter) || Maxiter <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(Idigits) || Idigits <= 0)
    stop("decimal point of 'Interval' should be greater than 0")
  if (!is.numeric(Pdigits) || Pdigits <= 0)
    stop("decimal point of 'Probability' should be greater than 0")

  list(Maxiter = Maxiter, convergence = convergence, Idigits = Idigits,
       Pdigits = Pdigits)
}
