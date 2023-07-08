
FrailtyPro <- function(N, q, x, d, a, b, lambda, theta, beta, vy, vd) {
  La <- (cumsum(lambda[order(vy)]))[rank(vy)]
  La <- matrix(La, a, b)
  A <- 1/theta + rowSums(d)
  BE <- array(rep(beta, each = a * b), c(a, b, length(beta)))
  C <- 1/theta + rowSums(La * exp(apply(x * BE, c(1, 2), sum)))
  AC <- matrix(A/C, a, b)
  E_0 <- as.vector(AC * exp(apply(x * BE, c(1, 2), sum)))
  SUM_0 <- cumsum((E_0[order(vy)])[seq(N, 1, -1)])
  SUM_0 <- (SUM_0[seq(N, 1, -1)])[rank(vy)]

  # beta #
  for (p in 1:q) {
    E_1 <- as.vector(AC * x[, , p] * exp(apply(x * BE, c(1, 2), sum)))
    AVE_X <- apply(abs(x), c(1, 2), sum)/abs(x[, , p])
    E_2 <- as.vector(AC * AVE_X * x[, , p]^2 * exp(apply(x * BE, c(1,
                                                                   2), sum)))

    SUM_1 <- cumsum((E_1[order(vy)])[seq(N, 1, -1)])
    SUM_1 <- (SUM_1[seq(N, 1, -1)])[rank(vy)]
    SUM_2 <- cumsum((E_2[order(vy)])[seq(N, 1, -1)])
    SUM_2 <- (SUM_2[seq(N, 1, -1)])[rank(vy)]

    DE_1 <- sum(d * x[, , p]) - sum(vd * SUM_1/SUM_0)
    DE_2 <- -sum(vd * SUM_2/SUM_0)

    beta[p] <- beta[p] - DE_1/DE_2
  }

  # th #
  Q01 <- a * (digamma(1/theta) + log(theta) - 1)/(theta^2) + sum(A/C -
                                                                   digamma(A) + log(C))/(theta^2)
  Q02 <- a * (3 - 2 * digamma(1/theta) - 2 * log(theta))/(theta^3) +
    2 * sum(digamma(A) - log(C) - A/C)/(theta^3) - a * trigamma(1/theta)/(theta^4)
  th <- theta - Q01/Q02
  if (th > 0) {
    theta <- th
  }

  # la #
  A <- 1/theta + rowSums(d)
  BE <- array(rep(beta, each = a * b), c(a, b, length(beta)))
  C <- 1/theta + rowSums(La * exp(apply(x * BE, c(1, 2), sum)))
  AC <- matrix(A/C, a, b)
  E_0 <- as.vector(AC * exp(apply(x * BE, c(1, 2), sum)))
  SUM_0 <- cumsum((E_0[order(vy)])[seq(N, 1, -1)])
  SUM_0 <- (SUM_0[seq(N, 1, -1)])[rank(vy)]
  lambda <- vd/SUM_0
  TB <- c(theta, beta)
  return(list(TB = TB, lambda = lambda))
}
