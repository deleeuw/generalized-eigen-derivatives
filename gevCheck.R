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
thetaa <- c(1, 1, 1)
thetab <- c(1, 1, 1)
theta <- c(thetaa, thetab)

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
                     aarrayl,
                     barrayl,
                     ainterceptl,
                     binterceptl,
                     typel = 0) {
  n <- nrow(aarrayl[, , 1])
  pal <- dim(aarrayl)[3]
  pbl <- dim(barrayl)[3]
  a <- matrix(ainterceptl, n, n)
  for (s in 1:pal) {
    a <- a + theta[s] * aarrayl[, , s]
  }
  b <- matrix(binterceptl, n, n)
  if (pbl > 0) {
    for (s in (pal + 1):(pal + pbl)) {
      b <- b + theta[s] * barrayl[, , s - pal]
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

for (type in  0:4) {
  h <- jacobian(
    checkFun,
    theta,
    aarrayl = aarray,
    barrayl = barray,
    ainterceptl = 0,
    binterceptl = 0,
    typel = type
  )
  print(h)
}

gevLim <- function(thetaa, thetab, aarray, barray, aintercept = 0, bintercept = 0, eps = 1e-6) {
  pa <- dim(aarray)[3]
  pb <- dim(barray)[3]  
  n <- nrow(aarray[, , 1]) 
  a <- matrix(aintercept, n, n)
  for (s in 1:pa) {
    a <- a + thetaa[s] * aarray[, , s]
  }
  b <- matrix(bintercept, n, n)
  if (pb > 0) {
    for (s in 1:pb) {
      b <- b + thetab[s] * barray[, , s]
    }
  }
  e <- myGeigen(a, b)
  x <- e$vectors
  l <- e$values
  dl <- matrix(0, n, pa + pb)
  dx <- array(0, c(n, pa + pb, n))
  for (s in 1:pa) {
    ae <- a + eps * aarray[, , s]
    ee <- myGeigen(ae, b)
    xe <- (ee$vectors - e$vectors) / eps
    le <- (ee$values - e$values) / eps
    dl[, s] <- le
    dx[, s, ] <- xe
  }
  for (s in 1:pb) {
    be <- b + eps * barray[, , s]
    ee <- myGeigen(a, be)
    xe <- (ee$vectors - e$vectors) / eps
    le <- (ee$values - e$values) / eps
    dl[, pa + s] <- le
    dx[, pa + s, ] <- xe
  }
  return(list(dl = dl, dx = dx))
}