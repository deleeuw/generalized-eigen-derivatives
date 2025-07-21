data(cattell, package = "psych")
cmat <- cattell
theta <- rep(.5, 12)
library(numDeriv)

matrixPower <- function(x, p) {
  h <- eigen(x)
  hvec <- h$vectors
  hval <- diag(h$values^p)
  return(tcrossprod(hvec %*% hval, hvec))
}

f1 <- function(x) {
  return((1 / x) + log(x) - 1)
}

df1 <- function(x) {
  return(-1 / (x^2) + (1 / x))
}

ddf1 <- function(x) {
  return(2 / (x^3) - (1 / (x^2)))
}

ff1 <- list(f1, df1, ddf1)

f2 <- function(x) {
  return(0.5 * (1 - x)^2)
}

df2 <- function(x) {
  return(-(1 - x))
}

ddf2 <- function(x) {
  return(x)
}

ff2 <- list(f2, df2, ddf2)

f3 <- function(x) {
  return(0.5 * (log(x))^2)
}

df3 <- function(x) {
  return(log(x) / x)
}

ddf3 <- function(x) {
  return((1 - log(x)) / (x^2))
}

ff3 <- list(f3, df3, ddf3)

f4 <- function(x) {
  return(0.5 * ((x - 1)^2) / (x^2))
}

df4 <- function(x) {
  return((x - 1) / (x^3))
}

ddf4 <- function(x) {
 return((3 - (2 * x)) / x^4)
}

ff4 <- list(f4, df4, ddf4)

swainDerivatives <- function(theta, cmat, ff, p) {
  n <- length(theta)
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  hvec <- h$vectors
  y <- cmat %*% hvec
  g <- drop((y[, -(1:p)]^2) %*% ff[[2]](hval[-(1:p)]))
  h <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (nu in (p + 1):n) {
        h[i, j] <- h[i, j] + ff[[3]](hval[nu]) * (y[i, nu]^2) * (y[j, nu]^2)
        for (eta in 1:n) {
          if (eta == nu) {
            next
          }
          h[i, j] <- h[i, j] - 2 * ff[[2]](hval[nu]) / (hval[eta] - hval[nu]) * y[i, nu] * y[i, eta] * y[j, nu] * y[j, eta]
        }
      }
    }
  }
  return(list(g = g, h = h))
}

swainFunction <- function(theta, cmat, ff, p) {
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  f <- sum(ff[[1]](hval[-(1:p)]))
  return(f)
}

swainGradient <- function(theta, cmat, ff, p) {
  n <- length(theta)
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  hvec <- h$vectors
  y <- cmat %*% hvec
  g <- drop((y[, -(1:p)]^2) %*% ff[[2]](hval[-(1:p)]))
  return(g)
}

swainHessian <- function(theta, cmat, ff, p) {
  n <- length(theta)
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  hvec <- h$vectors
  y <- cmat %*% hvec
  h <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (nu in (p + 1):n) {
        h[i, j] <- h[i, j] + ff[[3]](hval[nu]) * (y[i, nu]^2) * (y[j, nu]^2)
        for (eta in 1:n) {
          if (eta == nu) {
            next
          }
          h[i, j] <- h[i, j] - 2 * ff[[2]](hval[nu]) / (hval[eta] - hval[nu]) * y[i, nu] * y[i, eta] * y[j, nu] * y[j, eta]
        }
      }
    }
  }
  return(h)
}

swainNumDerivatives <- function(theta, cmat, ff, p) {
  cmat <- matrixPower(cmat, -1 / 2)
  theFunc <- function(theta) {
    return(sum(ff[[1]](eigen(
      cmat %*% diag(theta) %*% cmat
    )$values)[-(1:p)]))
  }
  g <- jacobian(theFunc, theta)
  h <- hessian(theFunc, theta)
  return(list(g = g, h = h))
}

swainNLM <- function(theta, cmat, ff, p) {
  res <- swainFunction(theta, cmat, ff, p)
  attr(res, "gradient") <- swainGradient(theta, cmat, ff, p)
  attr(res, "hessian") <- swainHessian(theta, cmat, ff, p)
  return(res)
}

#