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
  da <- dA(theta)
  db <- dB(theta)
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
    w <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dfar <- da[, , s] - li * db[, , s]
      dl[i, s] <- sum(xi * (dfar %*% xi))
      dx[, s, i] <- -w %*% dfar %*% xi - 0.5 * sum(xi * (db[, , s] %*% xi)) * xi
    }
  }
  return(list(
    values = l,
    vectors = x,
    dvalues = dl,
    dvectors = dx
  ))
}

gevNonlinearNum <- function(theta) {
  p <- length(theta)
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
  n <- length(theFunc(theta, 0))
  for (t in 0:n) {
    print(jacobian(func = theFunc, x = theta, ss = t))
  }
}
