source("gevNonlinear.R")

set.seed(12345)
aarray <- array(0, c(4, 4, 6))
barray <- array(0, c(4, 4, 6))
aarray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 4] <- 0
aarray[, , 5] <- 0
aarray[, , 6] <- 0
barray[, , 1] <- 0
barray[, , 2] <- 0
barray[, , 3] <- 0
barray[, , 4] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 5] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 6] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
bintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
p <- 6
n <- 4

dA <- function(theta, s) {
  return(aarray[, , s])
}

dB <- function(theta, s) {
  return(barray[, , s])
}

ddA <- function(theta, s, t) {
  return(matrix(0, n, n))
}

ddB <- function(theta, s, t) {
  return(matrix(0, n, n))
}

theA <- function(theta) {
  a <- aintercept
  for ( s in 1:p) {
    a <- a + theta[s] * aarray[, , s]
  }
  return(a)
}

theB <- function(theta) {
  p <- length(theta)
  b <- bintercept
  for (s in 1:p) {
    b <- b + theta[s] * barray[, , s]
  }
  return(b)
}
