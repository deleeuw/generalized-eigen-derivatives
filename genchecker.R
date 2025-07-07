library(numDeriv)

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
pa <- 3
pb <- 3
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

checkFun <- function(theta,
                     pal,
                     pbl,
                     aarrayl,
                     barrayl,
                     ainterceptl,
                     binterceptl,
                     typel = 0) {
  n <- nrow(aarray[, , 1])
  a <- matrix(aintercept, n, n)
  for (s in 1:pa) {
    a <- a + theta[s] * aarray[, , s]
  }
  b <- matrix(bintercept, n, n)
  if (pb > 0) {
    for (s in (pa + 1):(pa + pb)) {
      b <- b + theta[s] * barray[, , s - pa]
    }
  }
  e <- myGeigen(a, b)
  l <- e$values
  x <- e$vectors
  if (typel == 0) {
    return(l)
  } else {
    return(x[, typel])
  }
}

h <- jacobian(
  checkFun,
  theta,
  pal = pa,
  pbl = pb,
  aarrayl = aarray,
  barrayl = barray,
  ainterceptl = 0,
  binterceptl = 0,
  typel = 2
)

print(h)