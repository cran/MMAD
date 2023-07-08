
HessianF <- function(N, q, x, d, vd, a, b, lambda, theta, beta, vy) {
  # HessianF is used to calculate the Hessian matrix of the Gamma
  # Frailty model

  La <- (cumsum(lambda[order(vy)]))[rank(vy)]
  La <- matrix(La, a, b)
  A <- 1/theta + rowSums(d)
  BE <- array(rep(beta, each = a * b), c(a, b, length(beta)))
  C <- 1/theta + rowSums(La * exp(apply(x * BE, c(1, 2), sum)))
  AC <- matrix(A/C, a, b)
  E_0 <- as.vector(AC * exp(apply(x * BE, c(1, 2), sum)))
  SUM_0 <- cumsum((E_0[order(vy)])[seq(N, 1, -1)])
  SUM_0 <- (SUM_0[seq(N, 1, -1)])[rank(vy)]

  df <- matrix(0, q + 1, q + 1)

  mc <- rowSums(La * exp(apply(x * BE, c(1, 2), sum)))
  ber <- exp(apply(x * BE, c(1, 2), sum))

  df[1, 1] <- sum(rowSums(d)/(theta^2) + 2 * digamma(1/theta)/(theta^3) +
                    trigamma(1/theta)/(theta^4) - 2 * digamma(A)/(theta^3) - trigamma(A)/(theta^4) -
                    2 * mc/((theta^2) * (1 + theta * mc)) - A * (mc^2)/(1 + theta *
                                                                          mc)^2 + (2/(theta^3)) * log(1 + theta * mc))


  for (p in 1:q) {
    E_1 <- as.vector(AC * x[, , p] * exp(apply(x * BE, c(1, 2), sum)))
    SUM_1 <- cumsum((E_1[order(vy)])[seq(N, 1, -1)])
    SUM_1 <- (SUM_1[seq(N, 1, -1)])[rank(vy)]
    l_0r_p <- -vd * SUM_1/SUM_0^2
    L_0r_p <- (cumsum(l_0r_p[order(vy)]))[rank(vy)]

    df[(p + 1), 1] <- df[1, (p + 1)] <- sum(rowSums(d) * rowSums(ber *
                                                                   (x[, , p] * La + L_0r_p))/(1 + theta * mc) - (1 + theta * rowSums(d)) *
                                              mc * rowSums(ber * (x[, , p] * La + L_0r_p))/(1 + theta * mc)^2)

    for (j in 1:q) {
      E_1j <- AC * x[, , j] * exp(apply(x * BE, c(1, 2), sum))
      SUM_1j <- cumsum((E_1j[order(vy)])[seq(N, 1, -1)])
      SUM_1j <- (SUM_1j[seq(N, 1, -1)])[rank(vy)]
      l_0a_j <- -vd * SUM_1j/SUM_0^2
      L_0a_j <- (cumsum(l_0a_j[order(vy)]))[rank(vy)]

      E_2 <- AC * x[, , p] * x[, , j] * exp(apply(x * BE, c(1, 2),
                                                  sum))
      SUM_2 <- cumsum((E_2[order(vy)])[seq(N, 1, -1)])
      SUM_2 <- (SUM_2[seq(N, 1, -1)])[rank(vy)]
      E_4 <- 2 * E_1 * E_1j
      SUM_3 <- cumsum((E_4[order(vy)])[seq(N, 1, -1)])
      SUM_3 <- (SUM_3[seq(N, 1, -1)])[rank(vy)]
      l_0ra_p <- -vd * SUM_2/SUM_0^2 + vd * SUM_3/SUM_0^3
      L_0ra_p <- (cumsum(l_0ra_p[order(vy)]))[rank(vy)]

      df[j + 1, p + 1] <- df[p + 1, j + 1] <- sum((1 + theta * rowSums(d)) *
                                                    rowSums(ber * (x[, , p] * x[, , j] * La + x[, , p] * L_0a_j +
                                                                     x[, , j] * L_0r_p + L_0ra_p))/(1 + theta * mc) - theta *
                                                    (1 + theta * rowSums(d)) * rowSums(ber * (x[, , p] * La +
                                                                                                L_0r_p))/(1 + theta * mc) * rowSums(ber * (x[, , j] * La +
                                                                                                                                             L_0a_j))/(1 + theta * mc)) + sum(((l_0r_p * l_0a_j - l_0ra_p *
                                                                                                                                                                                  lambda)/lambda^2)[which(d > 0)])

    }

  }

  return(df)
}
