source("gevdNonlinear.R")

data(glass, package = "anacor")
theta <- as.vector(as.matrix(glass))
nr <- 7
nc <- 7
n <- nc + nr
p <- nc * nr

indi <- function(s, n) {
  j <- ceiling(s / n)
  i <- s - n * (j - 1)
  return(list(i = i, j = j))
}

ei <- function(i, n) {
  return(ifelse(i == 1:n, 1, 0))
}

theF <- function(theta) {
  return(matrix(theta, nr, nc))
}

theG <- function(theta) {
  f <- matrix(theta, nr, nr)
  return(diag(rowSums(f)))
}

theH <- function(theta) {
  f <- matrix(theta, nc, nc)
  return(diag(colSums(f)))
}

dF <- function(theta, s) {
  k <- indi(s, nr)
  return(outer(ei(k$i, nr), ei(k$j, nc)))
}

dG <- function(theta, s) {
  return(diag(rowSums(dF(theta, s))))
}

dH <- function(theta, s) {
  return(diag(colSums(dF(theta, s))))
}

ddF <- function(theta, s, t) {
  return(matrix(0, nr, nc))
}

ddG <- function(theta, s, t) {
  return(matrix(0, nr, nr))
}

ddH <- function(theta, s, t) {
  return(matrix(0, nc, nc))
}

theA <- function(theta) {
  a <- matrix(0, n, n)
  a[1:nr, (nr + 1):n] <- theF(theta)
  return(a + t(a))
}

theB <- function(theta) {
  b <- matrix(0, n, n)
  b[1:nr, 1:nr] <- theG(theta)
  b[(nr + 1):n, (nr + 1):n] <- theH(theta)
  return(b)
}

dA <- function(theta, s) {
  a <- matrix(0, n, n)
  a[1:nr, (nr + 1):n] <- dF(theta, s)
  return(a + t(a))
}

dB <- function(theta, s) {
  b <- matrix(0, n, n)
  b[1:nr, 1:nr] <- dG(theta, s)
  b[(nr + 1):n, (nr + 1):n] <- dH(theta, s)
  return(b)
}

ddA <- function(theta, s, t) {
  a <- matrix(0, n, n)
  a[1:nr, (nr + 1):n] <- ddF(theta, s, t)
  return(a + t(a))
}

ddB <- function(theta, s, t) {
  b <- matrix(0, n, n)
  b[1:nr, 1:nr] <- ddG(theta, s, t)
  b[(nr + 1):n, (nr + 1):n] <- ddH(theta, s, t)
  return(b)
}