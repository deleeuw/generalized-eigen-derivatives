library(numDeriv)

myGeigen <- function(a, b) {
  h <- eigen(solve(b, a))
  lbd <- h$values
  x <- h$vectors
  xbx <- apply(x, 2, function(z)
    sum(z * (b %*% z)))
  x <- x %*% diag(1 / sqrt(xbx))
  return(list(values = lbd, vectors = x))
}

dEigen <- function(a, b, da, db, eps = 1e-6) {
  n <- nrow(a)
  abeig <- myGeigen(a, b)
  abval <- abeig$values
  abvec <- abeig$vectors
  aeps <- a + eps * da
  beps <- b + eps * db
  abeigp <- myGeigen(aeps, beps)
  abvecp <- abeigp$vectors
  abvalp <- abeigp$values
  dvaln <- (abvalp - abval) / eps
  dvecn <- (abvecp - abvec) / eps
  dvalf <- rep(0, n)
  dvecf <- matrix(0, n, n)
  for (s in 1:n) {
    xs <- abvec[, s]
    ls <- abval[s]
    mm <- da - ls * db
    bb <- sum(xs * (db %*% xs))
    dvalf[s] <- sum(xs * (mm %*% xs))
    dg <- abval - ls
    dg[s] <- 1
    dg <- diag(1 / dg)
    dg[s, s] <- 0
    m1 <- abvec %*% dg %*% t(abvec)
    dvecf[, s] <- -(m1 %*% mm %*% xs + (1 / 2) * bb  * xs)
  }
  return(list(
    dvaln = dvaln,
    dvalf = dvalf,
    dvecn = dvecn,
    dvecf = dvecf
  ))
}

pEigenF <- function(a0, a1, a2, b0, b1, b2, theta, xi, eps = 1e-6) {
  n <- nrow(a0)
  a <- a0 + theta * a1 + (theta^2) * a2 / 2
  b <- b0 + xi * b1 + (xi^2) * b2 / 2
  abeig <- myGeigen(a, b)
  abval <- abeig$values
  abvec <- abeig$vectors
  dvalf <- rep(0, n)
  dvecf <- matrix(0, n, n)
  dtha <- a1 + theta * a2
  dxib <- b1 + xi * b2
  dvalf <- matrix(0, 2, n)
  for (s in 1:n) {
    xs <- abvec[, s]
    ls <- abval[s]
    dvalf[1, s] <- sum(xs * (dtha %*% xs))
    dvalf[2, s] <- -ls * sum(xs * (dxib %*% xs))
  }
  return(dvalf)
}

pEigenN <- function(a0, a1, a2, b0, b1, b2, theta, xi) {
  theFuncL <- function(x, aa0, aa1, aa2, bb0, bb1, bb2) {
    a <- a0 + x[1] * a1 + (x[1]^2) * a2 / 2
    b <- b0 + x[2] * b1 + (x[2]^2) * b2 / 2
    abeig <- myGeigen(a, b)
    return(abeig$values)
  }
  h <- jacobian(
    theFuncL,
    c(theta, xi),
    aa0 = a0,
    aa1 = a1,
    aa2 = a2,
    bb0 = b0,
    bb1 = b1,
    bb2 = b2
  )
  return(t(h))
}