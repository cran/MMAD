
ZIGPFx <- function(phi0, la, th) {
  m <- length(la)
  df <- matrix(0, 2 * m + 1, 2 * m + 1)
  a0 <- exp(-sum(la))

  r1 <- phi0 + (1 - phi0) * a0
  df[1, 1] <- -(1 - a0)^2/r1 - (1 - r1)/((1 - phi0)^2)
  df[2:(m + 1), 1] <- df[1, 2:(m + 1)] <- a0/r1
  df[(m + 2):(2 * m + 1), 1] <- df[1, (m + 2):(2 * m + 1)] <- 0

  for (i in 1:m) {
    for (j in 1:m) {
      df[j + 1, i + 1] <- df[i + 1, j + 1] <- phi0 * (1 - phi0) *
        a0/r1
      df[j + m + 1, i + 1] <- df[i + 1, j + m + 1] <- 0
      df[j + m + 1, i + m + 1] <- df[i + m + 1, j + m + 1] <- 0
    }

    df[i + 1, i + 1] <- -(1 - phi0)/la[i] + th[i] * (1 - phi0)/(la[i] +
                                                                  2 * th[i]) + phi0 * (1 - phi0) * a0/r1
    df[i + m + 1, i + m + 1] <- -la[i] * (1 - phi0)/(1 - th[i]) - 2 *
      la[i] * (1 - phi0)/(la[i] + 2 * th[i])
    df[i + m + 1, i + 1] <- df[i + 1, i + m + 1] <- -la[i] * (1 - phi0)/(la[i] +
                                                                           2 * th[i])
  }

  return(df)
}


ZIGPQx <- function(phi0, la, th) {
  m <- length(la)
  df <- matrix(0, 2 * m + 1, 2 * m + 1)
  dphi0 <- -1/phi0 - 1/(1 - phi0)
  dla <- -(1 - phi0)/la
  dpi <- -(1 - phi0) * la/(th * (1 - th))
  dth <- c(dphi0, dla, dpi)
  dQ <- diag(dth)
  return(dQ)
}


ZIGP_Rate <- function(phi0, la, th) {
  DF <- ZIGPFx(phi0, la, th)
  DGI <- solve(ZIGPQx(phi0, la, th))
  A <- DGI %*% DF
  rate <- 1 - min(eigen(A)$values)
  return(rate)
}
