library(numDeriv)

matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

myGeigen <- function(a, b) {
  h <- eigen(solve(b, a))
  lbd <- h$values
  x <- h$vectors
  xbx <- apply(x, 2, function(z)
    sum(z * (b %*% z)))
  x <- x %*% diag(1 / sqrt(xbx))
  return(list(values = Re(lbd), vectors = Re(x)))
}

gevdNonlinear <- function(theta) {
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
  for (nu in 1:n) {
    xnu <- x[, nu]
    lnu <- l[nu]
    mpl <- l - lnu
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wnu <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dsa <- dA(theta, s)
      dsb <- dB(theta, s)
      dsab <- dsa - lnu * dsb
      dl[nu, s] <- sum(xnu * (dsab %*% xnu))
      dx[, s, nu] <- -wnu %*% dsab %*% xnu - 0.5 * sum(xnu * (dsb %*% xnu)) * xnu
    }
  }
  ddl <- gevdHessianValues(theta, x, l)
  ddx <- gevdHessianVectors(theta, x, l, dx, dl)
  return(list(
    a = a,
    b = b,
    l = l,
    x = x,
    dl = dl,
    dx = dx,
    ddl = ddl,
    ddx = ddx
  ))
}

gevdHessianValues <- function(theta, x, l) {
  ddl <- array(0, c(p, p, n))
  for (nu in 1:n) {
    xnu <- x[, nu]
    lnu <- l[nu]
    mpl <- l - lnu
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wnu <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dsa <- dA(theta, s)
      dsb <- dB(theta, s)
      dsab <- dsa - lnu * dsb
      for (t in 1:p) {
        dta <- dA(theta, t)
        dtb <- dB(theta, t)
        dtab <- dta - lnu * dtb
        ddtp <- ddA(theta, s, t) - lnu * ddB(theta, s, t)
        accu <- -2 * sum(xnu * (dsab %*% wnu %*% dtab %*% xnu))
        accu <- accu + sum(xnu * (ddtp %*% xnu))
        accu <- accu - sum(xnu * (dtb %*% xnu)) * sum(xnu * (dsab %*% xnu))
        accu <- accu - sum(xnu * (dsb %*% xnu)) * sum(xnu * (dtab %*% xnu))
        ddl[s, t, nu] <- accu
      }
    }
  }
  return(ddl)
}

gevdHessianVectors <- function(theta, x, l, dx, dl) {
  ddx <- array(0, c(p, p, n, n))
  for (nu in 1:n) {
    xnu <- x[, nu]
    lnu <- l[nu]
    for (s in 1:p) {
      dsa <- dA(theta, s)
      dsb <- dB(theta, s)
      dsab <- dsa - lnu * dsb
      for (t in 1:p) {
        dta <- dA(theta, t)
        dtb <- dB(theta, t)
        dtab <- dta - lnu * dtb
        dstb <- ddB(theta, s, t)
        ddtp <- ddA(theta, s, t) - lnu * dstb
        dtxi <- dx[, t, nu]
        dtli <- dl[nu, t]
        accu <- 0
        accu <- accu + sum(xnu * (dsb %*% xnu)) * dtxi / 2
        accu <- accu + sum(xnu * (dstb %*% xnu)) * xnu / 2
        accu <- accu + sum(xnu * (dsb %*% dtxi)) * xnu
        for (eta in 1:n) {
          if (eta == nu) {
            next
          }
          xeta <- x[, eta]
          dtlj <- dl[eta, t]
          leta <- l[eta]
          dlij <- leta - lnu
          dtxj <- dx[, t, eta]
          djab <- sum(xeta * (dsab %*% xnu))
          accj <- 0
          accj <- accj + sum(dtxj * (dsab %*% xnu))
          accj <- accj + sum(xeta * (dsab %*% dtxi))
          accj <- accj + sum(xeta * (ddtp %*% xnu))
          accj <- accj - dtli * sum(xeta * (dsb %*% xnu))
          accu <- accu + (accj / dlij) * xeta
          accu <- accu - ((dtlj - dtli) / (dlij^2)) * djab * xeta
          accu <- accu + (djab / dlij) * dtxj
        }
        ddx[s, t, , nu] <- -accu
      }
    }
  }
  return(ddx)
}

gevdNonlinearNum <- function(theta) {
  theEunc <- function(theta) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    return(h$values)
  }
  theFunc <- function(theta, ss) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    return(h$vectors[, ss])
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
  a <- theA(theta)
  b <- theB(theta)
  h <- myGeigen(a, b)
  l <- h$values
  x <- h$vectors
  dl <- jacobian(func = theEunc, x = theta, method.args = list(r = 6))
  n <- nrow(theA(theta))
  dx <- array(0, c(n, p, n))
  ddl <- array(0, c(p, p, n))
  ddx <- array(0, c(p, p, n, n))
  for (nu in 1:n) {
    dx[, , nu] <- jacobian(
      func = theFunc,
      x = theta,
      method.args = list(r = 6),
      ss = nu
    )
    ddl[, , nu] <- hessian(
      func = theGunc,
      x = theta,
      method.args = list(r = 6),
      ss = nu
    )
  }
  for (i in 1:n) {
    for (nu in 1:n) {
      ddx[, , i, nu] <- hessian(
        func = theHunc,
        x = theta,
        method.args = list(r = 6),
        ss = nu,
        ii = i
      )
    }
  }
  return(list(
    a = a,
    b = b,
    l = l,
    x = x,
    dl = dl,
    dx = dx,
    ddl = ddl,
    ddx = ddx
  ))
}
