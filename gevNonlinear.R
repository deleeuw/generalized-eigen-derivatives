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

gevNonlinear <- function(theta) {
  p <- length(theta)
  a <- theA(theta)
  b <- theB(theta)
  h <- myGeigen(a, b)
  l <- h$values
  x <- h$vectors
  n <- length(l)
  dl <- matrix(0, n, p)
  dx <- array(0, c(n, p, n))
  ddl <- array(0, c(p, p, n))
  for (i in 1:n) {
    xi <- x[, i]
    li <- l[i]
    mpl <- l - li
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wi <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dsas <- dsA(theta, s)
      dsbs <- dsB(theta, s)
      dfar <- dsas - li * dsbs
      dl[i, s] <- sum(xi * (dfar %*% xi))
      dx[, s, i] <- -wi %*% dfar %*% xi - 0.5 * sum(xi * (dsbs %*% xi)) * xi
      for (t in 1:p) {
        dsat <- dsA(theta, t)
        dsbt <- dsB(theta, t)
        dfat <- dsat - li * dsbt
        dfit <- dstA(theta, t) - li * dstB(theta, t)
        accu <- -2 * sum(xi * (dfar %*% wi %*% dfat %*% xi))
        accu <- accu + sum(xi * (dfit %*% xi))
        accu <- accu - sum(xi * (dsbt %*% xi)) * sum(xi * (dfar %*% xi))
        accu <- accu - sum(xi * (dsbs %*% xi)) * sum(xi * (dfat %*% xi))
        ddl[s, t, i] <- accu
      }
    }
  }
  return(list(
    values = l,
    vectors = x,
    dvalues = dl,
    dvectors = dx,
    ddvalues = ddl
  ))
}

gevNonlinearNum <- function(theta) {
  theFunc <- function(theta, ss) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    if (ss == 0) {
      return(h$values)
    } else {
      return(h$vectors[, ss])
    }
  }
  theGunc <- function(theta, ss) {
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
