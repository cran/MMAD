# CZIGPMM

ZIGP_sample <- function(phi0, phi, la, th, n) {
  GP <- function(x, la0, th0) {
    b <- gamma(x + 1)
    if (th0 < 1 && th0 >= 0) {
      p <- exp(-la0 - th0 * x) * la0 * (la0 + th0 * x)^(x - 1)/b
    }
    if (th0 < 0 && th0 > -1 && la0 + th0 * x <= 0) {
      p <- 0
    }
    if (th0 < 0 && th0 > -1 && la0 + th0 * x > 0) {
      p <- exp(-la0 - th0 * x) * la0 * (la0 + th0 * x)^(x - 1)/b
    }
    return(p)
  }

  m <- length(la)
  Y <- matrix(0, n, m)
  P <- matrix(0, 100, m)
  for (j in 1:m) {
    for (i in 0:99) {
      P[i + 1, j] <- GP(i, la[j], th[j])
    }
    a <- 0:99
    X <- sample(a, n, P[, j], replace = TRUE)
    Z <- sample(c(0, 1), n, c(phi[j], 1 - phi[j]), replace = TRUE)
    Y[, j] <- X * Z
  }
  z0 <- sample(c(0, 1), n, c(phi0, 1 - phi0), replace = TRUE)
  zz <- matrix(z0, n, m)
  y <- zz * Y
  return(y)
}

n <- 600
m <- 2
phi0 <- 0.1
phi <- c(0.2, 0.5)  #rep(0.2,m)
la <- c(2, 7)  #rep(2,m)
th <- c(0.2, 0.3)  #rep(0.3,m)
y <- ZIGP_sample(phi0, phi, la, th, n)

phi0 <- 0.1
phi <- rep(0.1, m)
la <- rep(1, m)
th <- rep(0.1, m)

result <- CZIGPMM(y, phi0, phi, la, th)
summary(result)

# 2
y1 <- c(rep(0, 2789), rep(1, 224), rep(0, 726), rep(0, 307), rep(0, 171),
        rep(0, 76), rep(0, 32), rep(0, 16), rep(0, 15), rep(0, 9), rep(1, 212),
        rep(1, 149), rep(1, 85), rep(1, 50), rep(1, 35), rep(1, 13), rep(1,
        5), rep(1, 9), rep(2, 49), rep(2, 34), rep(2, 38), rep(2, 11),
        rep(2, 23), rep(2, 7), rep(2, 5), rep(2, 3), rep(2, 4), rep(3, 8),
        rep(3, 10), rep(3, 6), rep(3, 2), rep(3, 1), rep(3, 1), rep(3, 2),
        rep(4, 8), rep(4, 8), rep(4, 2), rep(4, 2), rep(4, 3), rep(4, 1), rep(5,
        3), rep(5, 3), rep(5, 2), rep(5, 1), rep(6, 2), rep(6, 1), rep(6, 3),
        rep(6, 1), rep(6, 2), rep(6, 2), rep(6, 1), rep(7, 1), rep(7,3),
        rep(7, 2), rep(7, 1), rep(7, 2), rep(7, 1), rep(7, 2), rep(8,1),
        rep(8, 1), rep(8, 1), rep(8, 1), rep(8, 1), rep(9, 1))

y2 <- c(rep(0, 2789), rep(0, 224), rep(1, 726), rep(2, 307), rep(3, 171),
        rep(4, 76), rep(5, 32), rep(6, 16), rep(7, 15), rep(8, 9), rep(1, 212),
        rep(2, 149), rep(3, 85), rep(4, 50), rep(5, 35), rep(6, 13), rep(7,5),
        rep(8, 9), rep(0, 49), rep(1, 34), rep(2, 38), rep(3, 11),
        rep(4, 23), rep(5, 7), rep(6, 5), rep(7, 3), rep(8, 4), rep(0, 8),
        rep(1, 10), rep(2, 6), rep(3, 2), rep(4, 1), rep(5, 1), rep(6, 2),
        rep(0, 8), rep(1, 8), rep(2, 2), rep(3, 2), rep(4, 3), rep(5, 1), rep(0,3),
        rep(1, 3), rep(2, 2), rep(4, 1), rep(0, 2), rep(2, 1), rep(3, 3), rep(4, 1),
        rep(5, 2), rep(6, 2), rep(8, 1), rep(0, 1), rep(2,3), rep(3, 2), rep(4, 1),
        rep(5, 2), rep(6, 1), rep(8, 2), rep(0,1), rep(1, 1), rep(2, 1), rep(4, 1),
        rep(6, 1), rep(8, 1))

yy <- data.frame(y1, y2)
phi0 <- 0.5
phi <- rep(0.5, 2)
la <- rep(5, 2)
th <- rep(0.5, 2)

result <- CZIGPMM(data = yy, phi0, phi, la, th)
summary(result)

# 3
result <- CZIGPMM(data = vijc, phi0, phi, la, th)
summary(result)

# GaFrailtyMM
library(survival)
result <- GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), data = kidney)
summary(result)
plot(result)

# IC2MM

sample <- function(la, n) {
  u <- runif(n, 0, 1)
  w <- runif(n, 0, 1)
  U <- -log(u)/la
  W <- -log(w)/la
  L <- pmin(U, W)
  R <- pmax(U, W)
  return(list(L = L, R = R))
}
n <- 10
y <- sample(0.5, n)
L <- y$L
R <- y$R
data <- data.frame(L, R)

result <- IC2MM(Surv(L, R, type = "interval2") ~ 1, data, control = IC2Control(Pdigits = 3))
summary(result)
plot(result)

# 2
result <- IC2MM(Surv(left, right, type = "interval2") ~ treatment, bcos,
                IC2Control(Pdigits = 3))
summary(result)
plot(result, col = c("red", "blue"))

# LTNMM

LTN_sample <- function(a, mu, si2, n) {
  x <- rep(0, n)
  for (i in 1:n) {
    repeat {
      z <- rnorm(1, mu, sqrt(si2))
      if (z >= a) {
        break
      }
    }
    x[i] <- z
  }
  return(x)
}

a <- 5
n <- 1000
mu <- 7
si2 <- 4
y <- LTN_sample(a, mu, si2, n)

result <- LTNMM(y ~ 1, a = 5)
summary(result)

# ZIGPMM

ZIGP_sample <- function(phi0, la, th, n) {
  GP <- function(x, la0, th0) {
    b <- gamma(x + 1)
    if (th0 < 1 && th0 >= 0) {
      p <- exp(-la0 - th0 * x) * la0 * (la0 + th0 * x)^(x - 1)/b
    }
    if (th0 < 0 && th0 > -1 && la0 + th0 * x <= 0) {
      p <- 0
    }
    if (th0 < 0 && th0 > -1 && la0 + th0 * x > 0) {
      p <- exp(-la0 - th0 * x) * la0 * (la0 + th0 * x)^(x - 1)/b
    }
    return(p)
  }

  m <- length(la)
  Y <- matrix(0, n, m)
  P <- matrix(0, 100, m)
  for (j in 1:m) {
    for (i in 0:99) {
      P[i + 1, j] <- GP(i, la[j], th[j])
    }
    a <- 0:99
    X <- sample(a, n, P[, j], replace = TRUE)
    Y[, j] <- X
  }
  z0 <- sample(c(0, 1), n, c(phi0, 1 - phi0), replace = TRUE)
  zz <- matrix(z0, n, m)
  y <- zz * Y
  return(y)
}

n <- 600
m <- 2
phi0 <- 0.4
la <- rep(9, m)
th <- rep(0.7, m)

y <- ZIGP_sample(phi0, la, th, n)

phi0 <- 0.1
la <- rep(1, m)
th <- rep(0.1, m)

result <- ZIGPMM(y, phi0, la, th)
summary(result)

# 2

phi0 <- 0.5
la <- rep(5, 2)
th <- rep(0.5, 2)

result <- ZIGPMM(vijc, phi0, la, th)
summary(result)

# 3
result <- ZIGPMM(cadi, phi0, la, th)
summary(result)

# CoxMM

sam <- function(n, be, th, la) {
  q <- length(be)
  x <- matrix(rnorm(n * q, -1, 1), n, q)
  u <- runif(n)
  t <- -log(u)/(la * exp(x %*% be))
  cen <- 3.7
  I <- 1 * (t <= cen)
  t <- pmin(t, cen)
  return(list(x = x, I = I, t = t))
}

n <- 200
q <- 3
be.true <- matrix(c(-0.5, 1, 2), ncol = 1)
t.true <- c(0.1, 0.5)
LA.TRUE <- 2 * t.true

be.in <- be.true * 0.5

da <- sam(n, be = be.true, la = 2)
data <- da$x
time <- da$t
status <- da$I
be <- be.in
X1 <- data[, 1]
X2 <- data[, 2]
X3 <- data[, 3]
data <- data.frame(X1, X2, X3, time = time, status)

result <- CoxMM(Surv(time, status) ~ X1 + X2 + X3, data, beat = be)
summary(result)
plot(result)

# 2
result <- CoxMM(Surv(time, status) ~ age + sex, lung)
summary(result)
plot(result)

