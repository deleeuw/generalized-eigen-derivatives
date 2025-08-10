source("gsvdTemplate.R")
library(numDeriv)
library(RSpectra)


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
  hval <- Re(h$values)
  ival <- order(hval, decreasing = TRUE)
  hval <- hval[ival]
  hvec <- Re(h$vectors[, ival])
  xbx <- apply(hvec, 2, function(z)
    sum(z * (b %*% z)))
  hvec <- hvec %*% diag(1 / sqrt(xbx))
  sgn <- sign(drop((1:n) %*% hvec))
  hvec <- hvec %*% diag(sgn)
  return(list(values = hval, vectors = hvec))
}



gsvdNonlinear <- function(theta) {
  a <- theA(theta)
  b <- theB(theta)
  h <- myGeigen(a, b)
  hval <- h$values
  hvec <- h$vectors
  xind <- 1:nr
  yind <- nr + 1:nc
  dl <- matrix(0, nc, p)
  dx <- array(0, c(nr, p, nc))
  dy <- array(0, c(nc, p, nc))
  for (s in 1:p) {
    ts <- crossprod(hvec[xind, ], dF(theta, s) %*% hvec[yind, ])
    us <- crossprod(hvec[xind, ], dG(theta, s) %*% hvec[xind, ])
    vs <- crossprod(hvec[yind, ], dH(theta, s) %*% hvec[yind, ])
    for (nu in 1:nc) {
      dl[nu, s] <- 2 * ts[nu, nu] - hval[nu] * (us[nu, nu] + vs[nu, nu])
      for (eta in 1:n) {
        if (eta == nu) {
          next
        }
        s1 <- 1 / (hval[eta] - hval[nu])
        s2 <- ts[eta, nu] + ts[nu, eta]
        s3 <- us[eta, nu] + vs[nu, eta]
        t1 <- s1 * (s2 - hval[nu] * s3) * hvec[xind, eta]
        dx[, s, nu] <- dx[, s, nu] - t1
        t1 <- s1 * (s2 - hval[nu] * s3) * hvec[yind, eta]
        t2 <- (us[nu, nu] + vs[nu, nu]) * hvec[yind, nu] / 2
        dy[, s, nu] <- dy[, s, nu] - t1
      }
      t2 <- (us[nu, nu] + vs[nu, nu]) * hvec[xind, nu] / 2
      dx[, s, nu] <- dx[, s, nu] - t2
      t2 <- (us[nu, nu] + vs[nu, nu]) * hvec[yind, nu] / 2
      dy[, s, nu] <- dy[, s, nu] - t2
    }
  }
  return(list(dl = dl, dx = dx, dy = dy))
}

gsvdNonlinearNum <- function(theta) {
  theFunc <- function(theta) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    return(h$values[1:nc])
  }
  theGunc <- function(theta) {
    a <- theA(theta)
    b <- theB(theta)
    h <- myGeigen(a, b)
    sgn <- sign(drop((1:n) %*% (h$vectors)))
    h$vectors <- h$vectors %*% diag(sgn)
    return(h$vectors[, nu])
  }
  dl <- jacobian(theFunc, theta)
  dx <- array(0, c(nr, p, nc))
  dy <- array(0, c(nc, p, nc))
  for (nu in 1:nc) {
    h <- jacobian(theGunc, theta)
    dx[, , nu] <- h[1:nr, ]
    dy[, , nu] <- h[nr + 1:nc, ]
  }
  return(list(dl = dl, dx = dx, dy = dy))
}
