source("gevNonlinear.R")

set.seed(12345)
aarray <- array(0, c(4, 4, 3))
barray <- array(0, c(4, 4, 3))
aarray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
bintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100


dA <- function(theta) {
  p <- length(theta)
  n <- dim(aarray[, , 1])[1]
  da <- array(0, c(n, n, p))
  for (s in 1:p) {
    da[, , s] <- aarray[, , 1]
    da[, , s] <- da[, , s] + theta * aarray[, , 2]
    da[, , s] <- da[, , s] + (theta^2) * aarray[, , 3]
  }
  return(da)
}

dB <- function(theta) {
  p <- length(theta)
  n <- dim(barray[, , 1])[1]
  db <- array(0, c(n, n, p))
  for (s in 1:p) {
    db[, , s] <- barray[, , 1]
    db[, , s] <- db[, , s] + theta * barray[, , 2]
    db[, , s] <- db[, , s] + (theta^2) * barray[, , 3]
  }
  return(db)
}

ddA <- function(theta) {
  daa <- aarray[, , 2]
  daa <- daa + 2 * theta * aarray[, , 3]
  return(daa)
}

ddB <- function(theta) {
  dbb <- barray[, , 2]
  dbb <- dbb + 2 * theta * barray[, , 3]
  return(dbb)
}

theA <- function(theta) {
  a <- aintercept
  a <- a + theta * aarray[, , 1]
  a <- a + ((theta^2) / 2) * aarray[, , 2]
  a <- a + ((theta^3) / 3) * aarray[, , 3]
  return(a)
}

theB <- function(theta) {
  b <- bintercept
  b <- b + theta * barray[, , 1]
  b <- b + ((theta^2) / 2) * barray[, , 2]
  b <- b + ((theta^3) / 3) * barray[, , 3]
  return(b)
}
