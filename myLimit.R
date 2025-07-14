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
theta <- c(1, 1, 1, 1, 1, 1)

myGeigen <- function(a, b) {
  h <- eigen(solve(b, a))
  lbd <- h$values
  x <- h$vectors
  xbx <- apply(x, 2, function(z)
    sum(z * (b %*% z)))
  x <- x %*% diag(1 / sqrt(xbx))
  return(list(values = lbd, vectors = x))
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

ei <- function(i, n) {
  return(ifelse(i == (1:n), 1, 0))
}

myLimit <- function(theta, i, j, eps = 1e-5) {
  a <- theA(theta)
  b <- theB(theta)
  h0 <- myGeigen(a, b)
  ddx <- matrix(0, p, p)
  for (s in 1:p) {
    for (t in 1:p) {
      th2 <- theta + eps * (ei(s, p) + ei(t, p))
      ths <- theta + eps * ei(s, p)
      tht <- theta + eps * ei(t, p)
      a2 <- theA(th2)
      b2 <- theB(th2)
      as <- theA(ths)
      bs <- theB(ths)
      at <- theA(tht)
      bt <- theB(tht)
      h2 <- myGeigen(a2, b2)
      hs <- myGeigen(as, bs)
      ht <- myGeigen(at, bt)
      ddx[s, t] <- (h2$vectors[i, j] - hs$vectors[i, j] - ht$vectors[i, j] + h0$vectors[i, j]) / (eps^2)
    }
  }
  return(ddx)
}
