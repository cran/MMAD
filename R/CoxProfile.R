CoxProfile <- function(x, d, y, beta, lambda, N, q) {
  # compute
  vd <- as.vector(d)
  vy <- as.vector(y)
  E_0 <- as.vector(exp(x %*% beta))
  SUM_0 <- cumsum((E_0[order(vy)])[seq(N, 1, -1)])
  SUM_0 <- (SUM_0[seq(N, 1, -1)])[rank(vy)]

  # beta #
  for (p in 1:q) {
    E_1 <- as.vector(x[, p] * exp(x %*% beta))
    AVE_X <- apply(abs(x), 1, sum)/abs(x[, p])
    E_2 <- as.vector(AVE_X * x[, p]^2 * exp(x %*% beta))

    SUM_1 <- cumsum((E_1[order(vy)])[seq(N, 1, -1)])
    SUM_1 <- (SUM_1[seq(N, 1, -1)])[rank(vy)]
    SUM_2 <- cumsum((E_2[order(vy)])[seq(N, 1, -1)])
    SUM_2 <- (SUM_2[seq(N, 1, -1)])[rank(vy)]

    DE_1 <- sum(d * x[, p]) - sum(vd * SUM_1/SUM_0)
    DE_2 <- -sum(vd * SUM_2/SUM_0)

    beta[p] <- beta[p] - DE_1/DE_2
  }

  # la #
  E_0 <- as.vector(exp(x %*% beta))
  SUM_0 <- cumsum((E_0[order(vy)])[seq(N, 1, -1)])
  SUM_0 <- (SUM_0[seq(N, 1, -1)])[rank(vy)]
  lambda <- vd/SUM_0
  return(list(beta = beta, lambda = lambda))
}
