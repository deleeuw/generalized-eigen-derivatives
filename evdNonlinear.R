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

evdNonlinear <- function(theta) {
  a <- theA(theta)
  h <- eigen(a)
  l <- h$values
  x <- h$vectors
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
      dl[nu, s] <- sum(xnu * (dsa %*% xnu))
      dx[, s, nu] <- -wnu %*% dsa %*% xnu
    }
    if (hessianl) {
      ddl <- evdHessianValues(theta, x, l)
    } else {
      ddl = NULL
    }
    if (hessianx) {
      ddx <- evdHessianVectors(theta, x, l, dx, dl)
    } else {
      ddx = NULL
    }
  }
  return(list(
    a = a,
    l = l,
    x = x,
    dl = dl,
    dx = dx,
    ddl = ddl,
    ddx = ddx
  ))
}

evdHessianValues <- function(theta, x, l) {
  ddl <- array(0, c(p, p, n))
  for (nu in 1:n) {
    xnu <- x[, nu]
    lnu <- l[nu]
    mpl <- l - lnu
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wnu <- x %*% mpl %*% t(x)
    for (s in 1:p) {
      dsa <- dA(theta, s)
      for (t in 1:p) {
        dta <- dA(theta, t)
        ddt <- ddA(theta, s, t)
        accu <- -2 * sum(xnu * (dsa %*% wnu %*% dta %*% xnu))
        accu <- accu + sum(xnu * (ddt %*% xnu))
        ddl[s, t, nu] <- accu
      }
    }
  }
  return(ddl)
}

evdHessianVectors <- function(theta, x, l, dx, dl) {
  ddx <- array(0, c(p, p, n, n))
  for (nu in 1:n) {
    xnu <- x[, nu]
    lnu <- l[nu]
    for (s in 1:p) {
      dsa <- dA(theta, s)
      for (t in 1:p) {
        dta <- dA(theta, t)
        ddtp <- ddA(theta, s, t)
        dtxi <- dx[, t, nu]
        dtli <- dl[nu, t]
        accu <- 0
        for (eta in 1:n) {
          if (eta == nu) {
            next
          }
          xeta <- x[, eta]
          dtlj <- dl[eta, t]
          leta <- l[eta]
          dlij <- leta - lnu
          dtxj <- dx[, t, eta]
          djab <- sum(xeta * (dsa %*% xnu))
          accj <- 0
          accj <- accj + sum(dtxj * (dsa %*% xnu))
          accj <- accj + sum(xeta * (dsa %*% dtxi))
          accj <- accj + sum(xeta * (ddtp %*% xnu))
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

evdNonlinearNum <- function(theta) {
  theEunc <- function(theta) {
    a <- theA(theta)
    h <- eigen(a)
    return(h$values)
  }
  theFunc <- function(theta, ss) {
    a <- theA(theta)
    h <- eigen(a)
    return(h$vectors[, ss])
  }
  theGunc <- function(theta, ss) {
    a <- theA(theta)
    h <- eigen(a)
    return(h$values[ss])
  }
  theHunc <- function(theta, ss, ii) {
    a <- theA(theta)
    h <- eigen(a)
    return(h$vectors[ii, ss])
  }
  a <- theA(theta)
  h <- eigen(a)
  l <- h$values
  x <- h$vectors
  dl <- jacobian(func = theEunc,
                 x = theta,
                 method.args = list(r = 6))
  n <- nrow(theA(theta))
  dx <- array(0, c(n, p, n))
  if (hessianl) {
    ddl <- array(0, c(p, p, n))
  } else {
    ddl = NULL
  }
  if (hessianx) {
    ddx <- array(0, c(p, p, n, n))
  } else {
    ddx <- NULL
  }
  for (nu in 1:n) {
    dx[, , nu] <- jacobian(
      func = theFunc,
      x = theta,
      method.args = list(r = 6),
      ss = nu
    )
    if (hessianl) {
      ddl[, , nu] <- hessian(
        func = theGunc,
        x = theta,
        method.args = list(r = 6),
        ss = nu
      )
    }
  }
  if (hessianx) {
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
  }
  return(list(
    a = a,
    l = l,
    x = x,
    dl = dl,
    dx = dx,
    ddl = ddl,
    ddx = ddx
  ))
}
