source("gevNonlinear.R")
library(datasets)
freq <- HairEyeColor


N <- sum(freq)
kcat <- dim(freq)
m <- length(kcat)
p <- prod(kcat)
n <- sum(kcat)

aarray <- array(0, c(n, n, p))
barray <- array(0, c(n, n, p))

ei <- function(i, n) {
  return(ifelse(i == 1:n, 1, 0))
}

theta <- NULL
l <- 1
for (i in 1:kcat[1]) {
  for (j in 1:kcat[2]) {
    for (k in 1:kcat[3]) {
      theta <- c(theta, HairEyeColor[i, j, k])
      indi <- c(ei(i, kcat[1]), ei(j, kcat[2]), ei(k, kcat[3]))
      aarray[, , l] <- outer(indi, indi)
      barray[, , l] <- m * diag(indi)
      l <- l + 1
    }
  }
}

theta <- theta / N

aintercept = 0
bintercept = 0

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
