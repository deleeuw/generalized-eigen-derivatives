source("evdNonlinear.R")
set.seed(12345)
aintercept <- diag(1:4)
aperturber <- crossprod(matrix(rnorm(400), 100, 4)) / 100
p <- 1
n <- 4
hessianl <- TRUE
hessianx <- TRUE

theA <- function(theta) {
  return(aintercept + theta * aperturber)
}

dA <- function(theta, s) {
  return(aperturber)
}

ddA <- function(theta, s, t) {
  return(matrix(0, n, n))
}
