source("gevdNonlinear.R")
set.seed(12345)
aintercept <- diag(1:4)
aperturber <- crossprod(matrix(rnorm(400), 100, 4)) / 100
p <- 1
n <- 4

theA <- function(theta) {
  return(aintercept + theta * aperturber)
}

theB <- function(theta) {
  return(diag(n))
}

dA <- function(theta, s) {
  return(aperturber)
}

dB <- function(theta, s) {
  return(matrix(0, n, n))
}

ddA <- function(theta, s, t) {
  return(matrix(0, n, n))
}

ddB <- function(theta, s, t) {
  return(matrix(0, n, n))
}