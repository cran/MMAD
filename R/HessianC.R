
HessianC <- function(N, q, x, d, vd, lambda, beta, vy) {
  # HessianC is used to calculate the Hessian matrix of the Cox
  # model

  La <- (cumsum(lambda[order(vy)]))[rank(vy)]
  E_0 <- as.vector(exp(x %*% beta))
  SUM_0 <- cumsum((E_0[order(vy)])[seq(N, 1, -1)])
  SUM_0 <- (SUM_0[seq(N, 1, -1)])[rank(vy)]

  df <- matrix(0, q, q)

  ber <- exp(x %*% beta)
  for (p in 1:q) {
    E_1 <- as.vector(x[, p] * exp(x %*% beta))
    SUM_1 <- cumsum((E_1[order(vy)])[seq(N, 1, -1)])
    SUM_1 <- (SUM_1[seq(N, 1, -1)])[rank(vy)]
    l_0r_p <- -vd * SUM_1/SUM_0^2
    L_0r_p <- (cumsum(l_0r_p[order(vy)]))[rank(vy)]

    for (j in 1:q) {
      E_1j <- as.vector(x[, j] * exp(x %*% beta))
      SUM_1j <- cumsum((E_1j[order(vy)])[seq(N, 1, -1)])
      SUM_1j <- (SUM_1j[seq(N, 1, -1)])[rank(vy)]
      l_0a_j <- -vd * SUM_1j/SUM_0^2
      L_0a_j <- (cumsum(l_0a_j[order(vy)]))[rank(vy)]

      E_2 <- x[, p] * x[, j] * exp(x %*% beta)
      SUM_2 <- cumsum((E_2[order(vy)])[seq(N, 1, -1)])
      SUM_2 <- (SUM_2[seq(N, 1, -1)])[rank(vy)]
      E_4 <- 2 * E_1 * E_1j
      SUM_3 <- cumsum((E_4[order(vy)])[seq(N, 1, -1)])
      SUM_3 <- (SUM_3[seq(N, 1, -1)])[rank(vy)]
      l_0ra_p <- -vd * SUM_2/SUM_0^2 + vd * SUM_3/SUM_0^3
      L_0ra_p <- (cumsum(l_0ra_p[order(vy)]))[rank(vy)]

      df[j, p] <- df[p, j] <- sum(rowSums(ber * (x[, p] * x[, j] *
                                                   La + x[, p] * L_0a_j + x[, j] * L_0r_p + L_0ra_p))) + sum(((l_0r_p *
                                                                                                                 l_0a_j - l_0ra_p * lambda)/lambda^2)[which(d > 0)])

    }

  }

  return(df)
}
