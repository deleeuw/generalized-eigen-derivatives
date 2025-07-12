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
  }
  ddx <- array(0, c(n, n, p, p))
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
          dtwi <- dtwi - ((dl[j, t] -dl[i, t]) / ((lj - li)^2)) * outer(xj, xj)
        }
        dsat <- dsA(theta, t)
        dsbt <- dsB(theta, t)
        dfat <- dsat - li * dsbt
        accu <- 0
        accu <- accu - dtwi %*% dfas %*% xi
        accu <- accu - wi %*% (dstA(theta, s, t) - li * dstB(theta, s, t)) %*% xi
        accu <- accu + dl[i, t] * wi %*% dsbs %*% xi
        accu <- accu - wi %*% dfas %*% dsbt %*% xi
        accu <- accu - sum(dx[, t, i] * (dsbs %*% xi)) * xi
        accu <- accu - 0.5 * sum(xi * (dsbs %*% xi)) * dx[, t, i]
        accu <- accu - 0.5 * sum(xi * (dstB(theta, s, t) %*% xi)) * xi
      }
    }
  }
  return(list(
    values = l,
    vectors = x,
    dvalues = dl,
    dvectors = dx,
    ddvectors = ddx
  ))
}

gevHessianNum <- function(theta) {
  theFunc <- function(theta, ss) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    return(h$values[ss])
  }
  dl <- jacobian(func = theFunc, x = theta, ss = 0)
  n <- nrow(theA(theta))
  dx <- array(0, c(n, p, n))
  ddl <- array(0, c(p, p, n))
  for (t in 1:n) {
    dx[, , t] <- jacobian(func = theFunc, x = theta, ss = t)
    ddl[, , t] <- hessian(func = theGunc, x = theta, ss = t)
  }
  return(list(dl = dl, dx = dx, ddl = ddl))
}
