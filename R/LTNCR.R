
LTNFx <- function(a, mu, sigma) {
  a1 <- (a - mu)/(sqrt(sigma))
  dphi <- -a1 * exp(-0.5 * a1^2)/sqrt(2 * 3.14159265)
  c <- 1 - pnorm(a1)
  df <- matrix(0, 2, 2)
  df[1, 1] <- -1/sigma + dphi/(c * sigma) + (dnorm(a1))^2/(sigma * c^2)
  df[1, 2] <- df[2, 1] <- (a1 * dphi - dnorm(a1))/(2 * c * sqrt(sigma^3)) +
    a1 * (dnorm(a1))^2/(2 * c^2 * sqrt(sigma^3))
  df[2, 2] <- -1/(2 * sigma^2) + (a1^2 * dphi - a1 * dnorm(a1))/(4 *
                                                                   c * sigma^2) + a1^2 * (dnorm(a1))^2/(4 * c^2 * sigma^2)

  return(df)
}


LTNQx <- function(a, mu, sigma) {
  a1 <- (a - mu)/(sqrt(sigma))
  c <- 1 - pnorm(a1)
  s1 <- (1 - c)/c
  tao <- exp(-(a - mu)^2/(2 * sigma))/sqrt(2 * pi * sigma)
  g <- tao/pnorm(a1)

  dQ <- matrix(0, 2, 2)
  dQ[1, 1] <- -(1 + s1)/sigma
  dQ[1, 2] <- dQ[2, 1] <- s1 * g/sigma - dnorm(a1)/(c * sqrt(sigma^3))
  dQ[2, 2] <- -(1 + s1)/(2 * sigma^2) + s1 * (a - mu) * g/(sigma^2) -
    a1 * dnorm(a1)/(c * sigma^2)

  return(dQ)
}

LTN_CRate <- function(a, mu, sigma) {
  DF <- LTNFx(a, mu, sigma)
  DGI <- solve(LTNQx(a, mu, sigma))
  A <- DGI %*% DF
  rate <- 1 - min(eigen(A)$values)
  return(rate)
}
