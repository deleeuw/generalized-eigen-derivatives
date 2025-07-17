source("gevNonlinear.R")

set.seed(12345)
aarray <- array(0, c(4, 4, 6))
barray <- array(0, c(4, 4, 6))
aarray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
bintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
n <- 4
p <- 1

dA <- function(theta, s = 1) {
  if (is.null(aarray)) {
    return(0)
  }
  das <- aarray[, , 1]
  for (r in 1:(d - 1)) {
    das <- das + (theta^r) * aarray[, , r + 1]
  }
  return(das)
}

dB <- function(theta, s = 1) {
  if (is.null(barray)) {
    return(0)
  }
  dbs <- barray[, , 1]
  for (r in 1:(d - 1)) {
    dbs <- dbs + (theta^r) * barray[, , r + 1]
  }
  return(dbs)
}

ddA <- function(theta, s = 1, t = 1) {
  daa <- aarray[, , 2]
  for (r in 1:(d - 2)) {
    daa <- daa + (r + 1) * (theta^r) * aarray[, , r + 2]
  }
  return(daa)
}

ddB <- function(theta, s = 1, t = 1) {
  dbb <- barray[, , 2]
  for (r in 1:(d - 2)) {
    dbb <- dbb + (r + 1) * (theta^r) * barray[, , r + 2]
  }
  return(dbb)
}

theA <- function(theta) {
  a <- aintercept
  for (r in 1:d) {
    a <- a + ((theta^r) / r) * aarray[, , r]
  }
  return(a)
}

theB <- function(theta) {
  b <- bintercept
  for (r in 1:d) {
    b <- b + ((theta^r) / r) * barray[, , r]
  }
  return(b)
}
