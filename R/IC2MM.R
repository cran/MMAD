#' MM algorithm based on the AD method for case II interval-censored data
#'
#'
#' @description The \code{IC2MM} function is used to calculate the case II interval-censored data model. A failure time study that consists of \eqn{n} independent subjects from a
#' homogeneous population with survival function \eqn{S_{(t)}}. Let \eqn{T_{i}} denote the survival time, and \eqn{i=1, \ldots, n}. Suppose that interval-censored data on the
#' \eqn{T_i} are observed and given by
#'
#' \deqn{Y_{obs} = \{ (L_{i}, R_{i}];  i=1, \ldots, n \} }
#'
#' where \eqn{T_i \in (L_{i}, R_{i}] }. Let \eqn{ \{s_i \}_{j=0}^{m} } denote the unique ordered elements of \eqn{ {0, L_{i}, R_{i}, i=1, \ldots, n } }.
#' Take \eqn{ \alpha_{ij} = I(s_{j} \in (L_{i}, R_{i}] ) } and \eqn{p_{j}= S(s_{j-1}) - S(s_{j}), j= 1, \ldots, m }. The log-likelihood function is
#'
#' \deqn{ \ell( {p} | Y_{obs}) = \sum_{i=1}^{n} \log (S(L_{i}) - S(R_{i}) ) = \sum_{i=1}^{n} \log \left( \sum_{j=1}^{m} \alpha_{ij} p_{j} \right)}
#'
#' where \eqn{{p} = (p_1, \ldots, p_m)^{T}  } and \eqn{ \sum_{j=1}^{m} p_{j} = 1 , p_{j} \geqslant 0}.
#'
#' @param formula A formula object, which contains on the left hand side an object of type = 'interval2' of the type \code{Surv}
#' e.g. \code{formula=Surv(L,R, type = 'interval2') ~ 1}
#' @param data A \code{data.frame} in which to interpret the variables named in the formula.
#' @param ... Additional arguments,  e.g. \code{control=IC2Control()}
#'
#' @details The \code{IC2MM} function allows the distributions for multiple strata of dataset to be stored as one \code{IC2} object, e.g. \code{data=bcos}.
#'
#' @return An object of class \code{IC2MM} that contains the following fields: \code{error}: convergence result; \code{strata}: dimensions of each \code{df_tab};
#' \code{s}: unique ordered elements of \eqn{ {0, L_{i}, R_{i}, Inf} }, if more than one strata, elements are concatenated; \code{S}: the survival function, if more than one strata, values are concatenated;
#' \code{df_tab}: the dataframe of survival intervals and survival probabilities for each interval, if more than one strata, dataframes are concatenated.
#' @seealso \code{\link{IC2Pro}}
#' @export
#'
#' @references Tian G.L., Huang X.F. and Xu, J.(2019). 'An assembly and decomposition approach for constructing separable minorizing functions in a class of MM algorithms.' \emph{Statistica Sinica} \strong{29}(2), 961-982.
#'
#' @examples
#' library(survival)
#' L <- c(1.5, 0.1, 1.5, 0.5, 0.4, 0.2, 0.9, 0.2, 0.08, 1.9)
#' R <- c(2.1, 2.9, 2.7, 1.9, 1.3, 1.4, 2.3, 0.5, 1.5, 4.6 )
#' data <- data.frame(L, R)
#' IC2MM(Surv(L,R, type = 'interval2') ~ 1, data )
#'
#' IC2MM(Surv(L,R, type = 'interval2') ~ 1, data, control=IC2Control(Pdigits=2) )
#'
IC2MM <- function(formula, data, ...) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' is missing or incorrect")
  }
  if (missing(data) || !inherits(data, "data.frame")) {
    stop("'data' is missing or not an object of type data.frame")
  }

  # Most of this function is copied or slightly modified from
  # survfit Copied starting from here:
  call <- match.call()
  if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Surv")) ||
      inherits(formula, "Surv")) {
    formula <- eval(parse(text = paste(deparse(call[[2]]), 1, sep = "~")))
    environment(formula) <- parent.frame()
  }


  m <- match.call(expand.dots = FALSE)
  m$... <- NULL
  Terms <- terms(formula, "strata")
  ord <- attr(Terms, "order")
  if (length(ord) & any(ord != 1))
    stop("Interaction terms are not valid for this function")
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  n <- nrow(m)
  ## left-hand-side of formula is 'response'
  Y <- model.extract(m, "response")
  ## change survfit code next few lines allow response to be
  ## numeric vector, treated as known time to event
  if (!is.Surv(Y)) {
    if (is.numeric(Y) & is.vector(Y))
      Y <- Surv(Y, rep(1, length(Y))) else stop("Response must be a survival object or numeric vector")
  }
  casewt <- model.extract(m, "weights")
  if (is.null(casewt))
    casewt <- rep(1, n)
  if (!is.null(attr(Terms, "offset")))
    warning("Offset term ignored")
  ll <- attr(Terms, "term.labels")
  if (length(ll) == 0)
    X <- factor(rep(1, n)) else X <- strata(m[ll])
  ### end of copied code from survfit.  do a separate fit for each
  ### level of the factor
  group <- levels(X)
  nstrata <- length(group)

  sbind <- function(x, y) {
    if (is.vector(x) & is.vector(y)) {
      if (is.list(x) & is.list(y)) {
        out <- list(x, y)
      } else out <- c(x, y)
    } else if (is.data.frame(x) & is.data.frame(y))
      out <- rbind(x, y)

    return(out)
  }

  SurvLR <- function(x) {
    type <- attr(x, "type")
    if (type != "interval") {
      stop("Surv object type='interval' only supported")
    }

    L <- R <- x[, 1]
    R[x[, 3] == 0] <- Inf
    L[x[, 3] == 2] <- 0
    R[x[, 3] == 3] <- x[x[, 3] == 3, 2]

    out <- data.frame(L = L, R = R)
    return(out)
  }

  ## change original left-hand-side of formula to list with L and R
  ## vectors representing left and right endpoints
  Y <- SurvLR(Y)
  for (i in 1:nstrata) {
    tempout <- IC2Pro(Y$L[X == group[i]], Y$R[X == group[i]], ...)
    pf <- tempout$df_tab
    tempout$strata <- length(pf[, 1])
    names(tempout$strata) <- group[i]
    if (i == 1) {
      icout <- tempout
    } else {
      icout <- mapply(sbind, icout, tempout)
    }
  }
  class(icout) <- "IC2"
  if (!is.null(attr(m, "na.action")))
    icout$na.action <- attr(m, "na.action")

  icout$call <- call
  return(invisible(icout))

}
