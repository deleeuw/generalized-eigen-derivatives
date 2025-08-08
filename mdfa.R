set.seed(12345)
library(numDeriv)
data(cattell, package = "psych")
cmat <- cattell
r <- 3
m <- 12
n <- 15
p <- m * (r + 1)
nfix <- 0

theta <- (1:48) / 100

fix <- matrix(c(1, 2, 0, 1, 3, 0, 2, 3, 0), 3, 3, byrow = TRUE)

indi <- function(s) {
  if (s > (m * r)) {
    i <- s - (m * r)
    j <- r + i
  } else {
    j <- ceiling(s / m)
    i <- s - (j - 1) * m
  }
  return(c(i, j))
}

jndi <- function(i, j) {
  return(i + (j - 1) * m)
}

ei <- function(i, n) {
  return(ifelse(i == 1:n, 1, 0))
}


theT <- function(theta) {
  tmat <- matrix(0, m, m + r)
  for (s in 1:p) {
    if (s <= (m * r)) {
      ij <- indi(s)
      i <- ij[1]
      j <- ij[2]
    } else {
      i <- s - (m * r)
      j <- i + r
    }
    tmat[i, j] <- theta[s]
  }
  if (nfix > 0) {
    for (k in 1:nfix) {
      tmat[fix[k, 1], fix[k, 2]] <- fix[k, 3]
    }
  }
  return(tmat)
}

dT <- function(theta, s) {
  if (s <= (m * r)) {
    ij <- indi(s)
    i <- ij[1]
    j <- ij[2]
  } else {
    i <- s - (m * r)
    j <- i + r
  }
  return(outer(ei(i, m), ei(j, n)))
}

theA <- function(theta) {
  tmat <- theT(theta)
  return(crossprod(tmat, cmat %*% tmat))
}

dA <- function(theta, s) {
  tmat <- theT(theta)
  ts <- dT(theta, s)
  qs <- crossprod(ts, cmat %*% tmat)
  qs <- qs + t(qs)
  return(qs)
}

mdfaEigenDerivatives <- function(theta) {
  amat <- theA(theta)
  heig <- eigen(amat, symmetric = TRUE)
  hvec <- heig$vectors
  hval <- heig$values
  hinv <- 1 / sqrt(abs(hval))
  f <- sum(sqrt(hval[1:m]))
  g <- matrix(0, m, p)
  h <- array(0, c(m, p, p))
  halt <- array(0, c(m, p, p))
  for (nu in 1:m) {
    xnu <- hvec[, nu]
    lnu <- hval[nu]
    gil <- diag(ifelse((hval - lnu) == 0, 0, 1 / (hval - lnu)))
    ainv <- hvec %*% gil %*% t(hvec)
    for (s in 1:p) {
      ts <- dT(theta, s)
      qs <- dA(theta, s)
      g[nu, s] <- sum(xnu * qs %*% xnu)
      for (t in 1:p) {
        tt <- dT(theta, t)
        qt <- dA(theta, t)
        hsum <- 0
        for (eta in 1:n) {
          xeta <- hvec[, eta]
          leta <- hval[eta]
          if (eta == nu) {
            next
          }
          hsum <- hsum + sum(xnu * (qs %*% xeta)) * sum(xnu * (qt %*% xeta)) / (leta - lnu)
        }
        qst <- crossprod(ts, cmat %*% tt)
        qst <- qst + t(qst)
        h[nu, s, t] <- -2 * hsum + sum(xnu * (qst %*% xnu))
      }
    }
  }
  return(list(f = f, dl = g, ddl = h))
}

mdfaDerivatives <- function(theta) {
  amat <- theA(theta)
  heig <- eigen(amat, symmetric = TRUE)
  hvec <- heig$vectors[, 1:m]
  hval <- heig$values[1:m]
  hinv <- 1 / sqrt(abs(hval))
  h <- mdfaEigenDerivatives(theta)
  ds <- drop(hinv %*% h$dl) / 2
  dds <- matrix(0, p, p)
  for (nu in 1:m) {
    dds <- dds - ((hval[nu]^(-3 / 2)) * outer(h$dl[nu, ], h$dl[nu, ])) / 4
    dds <- dds + (hinv[nu] * h$ddl[nu, , ]) / 2
  }
  return(list(ds = ds, dds = dds))
}

mdfaEigenDerivativesNum <- function(theta) {
  theFunc <- function(theta, j) {
    amat <- theA(theta)
    hval <- eigen(amat, symmetric = TRUE)$values
    return(hval[j])
  }
  g <- matrix(0, m, p)
  h <- array(0, c(m, p, p))
  for (nu in 1:m) {
    g[nu, ] <- grad(theFunc, theta, j = nu, method.args = list(r = 4))
    h[nu, , ] <- hessian(theFunc, theta, j = nu, method.args = list(r = 4))
  }
  return(list(dl = g, ddl = h))
}

mdfaDerivativesNum <- function(theta) {
  theFunc <- function(theta) {
    amat <- theA(theta)
    hval <- eigen(amat, symmetric = TRUE)$values
    return(sum(sqrt(hval[1:m])))
  }
  ds <- grad(theFunc, theta)
  dds <- hessian(theFunc, theta)
  return(list(ds = ds, dds = dds))
}

