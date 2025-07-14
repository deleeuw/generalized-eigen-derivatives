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

gevHessian <- function(theta) {
  p <- length(theta)
  a <- theA(theta)
  b <- theB(theta)
  h <- myGeigen(a, b)
  l <- h$values
  x <- h$vectors
  n <- length(l)
  dl <- matrix(0, n, p)
  dx <- array(0, c(n, p, n))
  ddx <- array(0, c(p, p, n, n))
  for (i in 1:n) {
    xi <- x[, i]
    li <- l[i]
    mpl <- l - li
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wi <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dsas <- dsA(theta, s)
      dsbs <- dsB(theta, s)
      dfas <- dsas - li * dsbs
      dl[i, s] <- sum(xi * (dfas %*% xi))
      dx[, s, i] <- -wi %*% dfas %*% xi - 0.5 * sum(xi * (dsbs %*% xi)) * xi
    }
    for (s in 1:p) {
      dsas <- dsA(theta, s)
      dsbs <- dsB(theta, s)
      dfas <- dsas - li * dsbs
      for (t in 1:p) {
        dtwi <- 0
        for (j in 1:n) {
          xj <- x[, j]
          lj <- l[j]
          if (j == i) {
            next
          }
          dtxj <- dx[, t, j]
          dtwi <- dtwi + (outer(xj, dtxj) + outer(dtxj, xj)) / (lj - li)
          dtwi <- dtwi - ((dl[j, t] - dl[i, t]) / ((lj - li)^2)) * outer(xj, xj)
        }
        dsat <- dsA(theta, t)
        dsbt <- dsB(theta, t)
        dfat <- dsat - li * dsbt
        dxti <- dx[, t, i]
        accu <- 0
        accu <- accu - dtwi %*% dfas %*% xi
        accu <- accu - wi %*% (dstA(theta, s, t) - li * dstB(theta, s, t)) %*% xi
        accu <- accu + dl[i, t] * wi %*% dsbs %*% xi
        accu <- accu - wi %*% dfas %*% dxti
        accu <- accu - sum(dxti * (dsbs %*% xi)) * xi
        accu <- accu - 0.5 * sum(xi * (dsbs %*% xi)) * dxti
        accu <- accu - 0.5 * sum(xi * (dstB(theta, s, t) %*% xi)) * xi
        ddx[s, t, , i] <- accu
      }
    }
  }
  return(list(
    l = l,
    x = x,
    dl = dl,
    dx = dx,
    ddx = ddx
  ))
}

gevHessianNum <- function(theta) {
  theFunc <- function(theta, ii, jj) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    return(h$vectors[ii, jj])
  }
  ddx <- array(0, c(p, p, n, n))
  for (i in 1:n) {
    for (j in 1:n) {
      ddx[, , i, j] <- hessian(
        func = theFunc,
        x = theta,
        method.args = list(eps = 1e-6, r = 6),
        ii = i,
        jj = j
      )
    }
  }
  return(ddx)
}
