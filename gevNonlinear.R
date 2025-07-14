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
  ddx <- array(0, c(p, p, n, n))
  for (i in 1:n) {
    xi <- x[, i]
    li <- l[i]
    mpl <- l - li
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wi <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dsa <- dA(theta, s)
      dsb <- dB(theta, s)
      dsab <- dsa - li * dsb
      dl[i, s] <- sum(xi * (dsab %*% xi))
      dx[, s, i] <- -wi %*% dsab %*% xi - 0.5 * sum(xi * (dsb %*% xi)) * xi
    }
    for (s in 1:p) {
      dsa <- dA(theta, s)
      dsb <- dB(theta, s)
      dsab <- dsa - li * dsb
      for (t in 1:p) {
        dta <- dA(theta, t)
        dtb <- dB(theta, t)
        dtab <- dta - li * dtb
        ddtp <- ddA(theta, s, t) - li * ddB(theta, s, t)
        dtxi <- dx[, t, i]
        accu <- -2 * sum(xi * (dsab %*% wi %*% dtab %*% xi))
        accu <- accu + sum(xi * (ddtp %*% xi))
        accu <- accu - sum(xi * (dtb %*% xi)) * sum(xi * (dsab %*% xi))
        accu <- accu - sum(xi * (dsb %*% xi)) * sum(xi * (dtab %*% xi))
        ddl[s, t, i] <- accu
        dtwi <- 0
        for (j in 1:n) {
          if (j != i) {
            xj <- x[, j]
            lj <- l[j]
            dtxj <- dx[, t, j]
            dtwi <- dtwi + (outer(xj, dtxj) + outer(dtxj, xj)) / (lj - li)
            dtwi <- dtwi - ((dl[j, t] - dl[i, t]) / ((lj - li)^2)) * outer(xj, xj)
          }
        }
        accu <- 0
        accu <- accu - dtwi %*% dsab %*% xi
        accu <- accu - wi %*% ddtp %*% xi
        accu <- accu + dl[i, t] * wi %*% dsb %*% xi
        accu <- accu - wi %*% dsab %*% dtxi
        accu <- accu - sum(dtxi * (dsb %*% xi)) * xi
        accu <- accu - 0.5 * sum(xi * (dsb %*% xi)) * dtxi
        accu <- accu - 0.5 * sum(xi * (ddB(theta, s, t) %*% xi)) * xi
        ddx[s, t, , i] <- accu
      }
    }
  }
  return(list(
    values = l,
    vectors = x,
    dvalues = dl,
    dvectors = dx,
    ddvalues = ddl,
    ddvectors = ddx
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
  theHunc <- function(theta, ss, ii) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    return(h$vectors[ii, ss])
  }
  h <- myGeigen(theA(theta), theB(theta))
  l <- h$values
  x <- h$vectors
  dl <- jacobian(func = theFunc, x = theta, ss = 0)
  n <- nrow(theA(theta))
  dx <- array(0, c(n, p, n))
  ddl <- array(0, c(p, p, n))
  ddx <- array(0, c(p, p, n, n))
  for (t in 1:n) {
    dx[, , t] <- jacobian(
      func = theFunc,
      x = theta,
      method.args = list(eps = 1e-6, r = 6),
      ss = t
    )
    ddl[, , t] <- hessian(
      func = theGunc,
      x = theta,
      method.args = list(eps = 1e-6, r = 6),
      ss = t
    )
  }
  for (i in 1:n) {
    for (t in 1:n) {
      ddx[, , i, t] <- hessian(
        func = theHunc,
        x = theta,
        method.args = list(eps = 1e-6, r = 6),
        ss = t,
        ii = i
      )
    }
  }
  return(list(
    values = l,
    vectors = x,
    dvalues = dl,
    dvectors = dx,
    ddvalues = ddl,
    ddvectors = ddx
  ))
}
