set.seed(12345)
amat <- crossprod(matrix(rnorm(400), 100, 4)) / 100
bmat <- crossprod(matrix(rnorm(400), 100, 4)) / 100
theta <- c(amat[outer(1:4, 1:4, ">=")], bmat[outer(1:4, 1:4, ">=")])
n <- 4
off <- 10
p <- 20
hessianl <- TRUE
hessianx <- TRUE

indi <- function(k, n) {
  nn <- 1:n
  y <- (nn * n) - choose(nn, 2) - (n - nn)
  j <- max(which(k >= y))
  i <- k + n - (((2 * n) + 1) * j  - j^2) / 2
  return(c(i, j))
}

theA <- function(theta) {
  a <- matrix(0, 4, 4)
  for (s in 1:10) {
    ij <- indi(s, 4)
    i <- ij[1]
    j <- ij[2]
    if (i == j) {
      a[i, i] <- theta[s]
    } else {
      a[i, j] <- a[j, i] <- theta[s]
    }
  }
  return(a)
}

theB <- function(theta) {
  b <- matrix(0, 4, 4)
  for (s in 1:10) {
    ij <- indi(s, 4)
    i <- ij[1]
    j <- ij[2]
    if (i == j) {
      b[i, i] <- theta[off + s]
    } else {
      b[i, j] <- b[j, i] <- theta[off + s]
    }
  }
  return(b)
}

dA <- function(theta, s) {
  da <- matrix(0, 4, 4)
  if (s > off) {
    return(da)
  }
  ij <- indi(s, 4)
  i <- ij[1]
  j <- ij[2]
  if (i == j) {
    da[i, i] <- 1
  } else {
    da[i, j] <- da[j, i] <- 1
  }
  return(da)
}

dB <- function(theta, s) {
  db <- matrix(0, 4, 4)
  if (s <= off) {
    return(db)
  }
  s <- s - off
  ij <- indi(s, 4)
  i <- ij[1]
  j <- ij[2]
  if (i == j) {
    db[i, i] <- 1
  } else {
    db[i, j] <- db[j, i] <- 1
  }
  return(db)
}

ddA <- function(theta, s, t) {
  return(matrix(0, 4, 4))
}

ddB <- function(theta, s, t) {
  return(matrix(0, 4, 4))
}
