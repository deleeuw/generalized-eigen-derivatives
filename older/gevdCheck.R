library(numDeriv)
source("gevdUtils.R")

partialCheck <- function(theta, a, b, s) {
  myFuncEval <- function(x, ain, bin, ss) {
    q <- length(ain)
    n <- nrow(ain[[1]])
    aa <- makeMatrix(x, ain)
    bb <- makeMatrix(x, bin)
    eig <- myGeigen(aa, bb)
    return(eig$val[ss])
  }
  myFuncEvec <- function(x, ain, bin, ss) {
    q <- length(ain)
    n <- nrow(ain[[1]])
    aa <- matrix(0, n, n)
    bb <- matrix(0, n, n)
    aa <- makeMatrix(x, ain)
    bb <- makeMatrix(x, bin)
    eig <- myGeigen(aa, bb)
    return(drop(eig$vec[, ss]))
  }
  grad <- grad(
    myFuncEval,
    x = theta,
    ain = a,
    bin = b,
    ss = s
  )
  jaco <- jacobian(
    myFuncEvec,
    x = theta,
    ain = a,
    bin = b,
    ss = s
  )
  return(list(dl = grad, dx = jaco))
}

hessianCheckEval <- function(theta, a, b, s) {
  myFuncEval <- function(x, a, b, s) {
    q <- length(a)
    n <- nrow(a[[1]])
    aa <- makeMatrix(x, a)
    bb <- makeMatrix(x, b)
    eig <- myGeigen(aa, bb)
    return(eig$val[s])
  }
  hess <- hessian(
    myFuncEval,
    x = theta,
    a = a,
    b = b,
    s = s
  )
  return(hess)
}

hessianCheckEvec <- function(theta, a, b, k, s) {
  myFuncEvec <- function(x, a, b, k, s) {
    q <- length(a)
    n <- nrow(a[[1]])
    aa <- matrix(0, n, n)
    bb <- matrix(0, n, n)
    aa <- makeMatrix(x, a)
    bb <- makeMatrix(x, b)
    eig <- myGeigen(aa, bb)
    return(eig$vec[k, s])
  }
  hess <- hessian(
    myFuncEvec,
    x = theta,
    a = a,
    b = b,
    k = k,
    s = s
  )
  return(hess)
}

