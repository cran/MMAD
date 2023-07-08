LogLik <- function(lambda, x, d, vy, vd, theta, beta, a, b) {
  La <- (cumsum(lambda[order(vy)]))[rank(vy)]
  La <- matrix(La, a, b)
  A <- 1/theta + rowSums(d)
  BE <- array(rep(beta, each = a * b), c(a, b, length(beta)))
  C <- 1/theta + rowSums(La * exp(apply(x * BE, c(1, 2), sum)))
  AC <- matrix(A/C, a, b)

  l1 <- sum(lgamma(A)) - a * (lgamma(1/theta) + log(theta)/theta) - sum(A *
                                                                          log(C))
  l2 <- sum(log(lambda[vd != 0]))
  l3 <- sum(d * (apply(x * BE, c(1, 2), sum)))

  ell <- l1 + l2 + l3

  return(ell)
}
