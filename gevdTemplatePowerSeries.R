source("gevdNonlinear.R")

set.seed(12345)
aarray <- array(0, c(4, 4, 3))
barray <- array(0, c(4, 4, 3))
aarray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
bintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
p <- 1
n <- 4
ma <- 3
mb <- 3
theta <- 1
hessianl <- TRUE
hessianx <- TRUE

theA <- function(theta) {
  a <- aintercept
  for (r in 1:ma) {
    a <- a + (theta^r) * aarray[, , r]
  }
  return(a)
}

theB <- function(theta) {
  b <- bintercept
  for (r in 1:mb) {
    b <- b + (theta^r) * barray[, , r]
  }
  return(b)
}

dA <- function(theta, s = 1) {
  if (ma == 0) {
    return(matrix(0, n, n))
  }
  das <- matrix(0, n, n)
  for (r in 1:ma) {
    das <- das + r * (theta^(r - 1)) * aarray[, , r]
  }
  return(das)
}

dB <- function(theta, s = 1) {
  if (mb == 0) {
    return(matrix(0, n, n))
  }
  dbs <- matrix(0, n, n)
  for (r in 1:mb) {
    dbs <- dbs + r * (theta^(r - 1)) * barray[, , r]
  }
  return(dbs)
}

ddA <- function(theta, s = 1, t = 1) {
  if (ma == 1) {
    return(matrix(0, n, n))
  }
  dda <- matrix(0, n, n)
  for (r in 2:ma) {
    dda <- dda + r * (r - 1) * (theta^(r - 2)) * aarray[, , r]
  }
  return(dda)
}

ddB <- function(theta, s = 1, t = 1) {
  if (mb == 1) {
    return(matrix(0, n, n))
  }
  ddb <- matrix(0, n, n)
  for (r in 2:ma) {
    ddb <- ddb + r * (r - 1) * (theta^(r - 2)) * barray[, , r]
  }
  return(ddb)
}
