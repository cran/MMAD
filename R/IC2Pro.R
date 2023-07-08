#' Calculate non-parametric estimate for case II interval censored survival function
#'
#' @param L The numeric vector of left endpoints of censoring interval, the first element of Surv when type=’interval2’.
#' @param R The numeric vector of right endpoints of censoring interval, the second element of Surv function when type=’interval2’.
#' @param control An object as created by \code{IC2Control}
#' @param ... Additional arguments
#'
#' @return An object of class \code{IC2Pro} that contains the following fields: \code{error}: convergence result; \code{strata}: dimensions of \code{df_tab};
#' \code{s}: unique ordered elements of \eqn{ {0, L_{i}, R_{i}, Inf} }; \code{S}: the survival function;
#' \code{df_tab}: the data frame of survival intervals and survival probabilities for each interval.
#' @seealso \code{\link{IC2Control}}
#' @export
#'
#' @references Tian G.L., Huang X.F. and Xu, J.(2019). 'An assembly and decomposition approach for constructing separable minorizing functions in a class of MM algorithms.' \emph{Statistica Sinica} \strong{29}(2), 961-982.
#'
#' @examples
#' L <- c(1.4, 1.5, 1.3, 0.9, 0.4, 0.2, 0.5, 0.03, 1.7, 0.2)
#' R <- c(2.2, 3, 2.4, 1.2, 2.8, 0.3, 1.6, 2.5, 2.6, 3.4)
#' IC2Pro(L, R, control=IC2Control())
#'
IC2Pro <- function(L, R, control = IC2Control(), ...) {
  if (length(L) == 0) {
    stop("The length of interval-censored data equals 0", call. = FALSE)
  }
  if (length(L) != length(R)) {
    stop("The amount of interval-censored data (L, R) is not equal",
         call. = FALSE)
  }

  Maxiter <- control$Maxiter
  convergence <- control$convergence
  Idigits <- control$Idigits
  Pdigits <- control$Pdigits
  LR <- c(L, R)
  ss <- sort(LR)
  s <- unique(c(0, ss, Inf))
  m <- length(s)
  n <- length(L)
  alpha <- matrix(0, n, m)
  for (i in 1:n) {
    alpha[i, ] <- (s > L[i]) * (s <= R[i])
  }

  p <- rep(1/m, m)
  pp <- matrix(rep(p, each = n), n, m)
  log_ell <- sum(log(rowSums(alpha * pp)))
  el <- c(log_ell)

  error <- 3
  for (k in 1:Maxiter) {
    if (error > convergence) {
      A <- rowSums(alpha * pp)
      AA <- t(matrix(rep(A, each = m), m, n))
      B <- colSums(alpha * pp/AA)
      p <- B/sum(B)

      pp <- matrix(rep(p, each = n), n, m)
      log_el <- sum(log(rowSums(alpha * pp)))
      el <- append(el, log_el)
      error <- abs(el[k + 1] - el[k])/(1 + abs(el[k]))

    }
  }

  S <- rep(0, m)
  for (i in 1:m) {
    S[i] <- 1 - sum((s <= s[i]) * p)
  }

  s <- round(s, Idigits)
  k <- length(s) - 1
  Lbracket <- rep("(", k)
  Rbracket <- rep("]", k)
  intmapL <- s[1:k]
  intmapR <- s[1:k + 1]
  intname <- paste(Lbracket, intmapL, ",", intmapR, Rbracket, sep = "")

  tab <- data.frame(Interval = intname, Probability = round(p[1:k + 1],
                                                            Pdigits))
  df_tab <- tab[which(tab$Probability > 0), ]
  # rename
  row.names(df_tab) <- c(1:length(df_tab$Probability))
  strata <- length(df_tab)

  result <- list(error = error, strata = strata, s = s, S = S, df_tab = df_tab)
  class(result) <- c("IC2", "list")
  return(result)

}
