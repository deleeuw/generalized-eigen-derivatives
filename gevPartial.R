library(numDeriv)

partialGeigen <- function(theta, a, b, s) {
  q <- length(a)
  n <- nrow(a[[1]])
  aa <- makeMatrix(theta, a)
  bb <- makeMatrix(theta, b)
  eig <- myGeigen(aa, bb)
  xmat <- eig$vec
  labd <- eig$val
  xs <- xmat[, s]
  ls <- labd[s]
  dl <- rep(0, q)
  dx <- matrix(0, n, q)
  for (p in 1:q) {
    dl[p] <- sum(xs * ((a[[p]] - ls * b[[p]]) %*% xs))
    aa <- crossprod(xmat, (a[[p]] - ls * b[[p]]) %*% xs)
    bb <- 1 / ((labd - ls) + ei(s, n)) - ei(s, n)
    cc <- aa * bb
    dd <- drop(crossprod(xs, b[[p]] %*% xs))
    dx[, p] <- -xmat %*% cc - dd * xs / 2.0
  }
  return(list(dl = dl, dx = dx))
}

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

partialElementEval <- function(a, b, s) {
  n <- nrow(a)
  eig <- myGeigen(a, b)
  xmat <- eig$vec
  labd <- eig$val
  xs <- xmat[, s]
  ls <- labd[s]
  dla <- 2 * outer(xs, xs)
  diag(dla) <- diag(dla) / 2
  dlb <- -2 * ls * outer(xs, xs)
  diag(dlb) <- diag(dlb) / 2
  for (t in 1:n) {
    num <- outer(xs, xmat[, t])
    num <- num + t(num)
    den <- outer
  }
  return(list(dla = dla, dlb = dlb))
}

partialElementEvec <- function(a, b, k, s) {
  n <- nrow(a)
  eig <- myGeigen(a, b)
  xmat <- eig$vec
  labd <- eig$val
  xs <- xmat[, s]
  ls <- labd[s]
  dxa <- matrix(0, n, n)
  dxb <- matrix(0, n, n)
  for (t in 1:n) {
    if (t == s) {
      next
    }
    num <- outer(xs, xmat[, t])
    num <- num + t(num) 
    diag(num) <- xs * xmat[, t]
    num <- num / (labd[t] - ls)
    dxa <- dxa - num * xmat[k, t]
    dxb <- dxb + ls * num * xmat[k, t]
  }
  cer <- outer(xs, xs) * xmat[k, s]
  diag(cer) <- diag(cer) / 2
  dxb <- dxb - cer
  return(list(dxa = dxa, dxb = dxb))
}

hessianGeigenEval <- function(theta, a, b, s) {
  q <- length(a)
  n <- nrow(a[[1]])
  aa <- makeMatrix(theta, a)
  bb <- makeMatrix(theta, b)
  eig <- myGeigen(aa, bb)
  xmat <- eig$vec
  labd <- eig$val
  xs <- xmat[, s]
  ls <- labd[s]
  dl <- rep(0, q)
  dx <- matrix(0, n, q)
  for (p in 1:q) {
    dl[p] <- sum(xs * ((a[[p]] - ls * b[[p]]) %*% xs))
    aa <- crossprod(xmat, (a[[p]] - ls * b[[p]]) %*% xs)
    bb <- 1 / ((labd - ls) + ei(s, n)) - ei(s, n)
    cc <- aa * bb
    dd <- drop(crossprod(xs, b[[p]] %*% xs))
    dx[, p] <- -xmat %*% cc - dd * xs / 2.0
  }
  hess <- matrix(0, q, q)
  for (p in 1:q) {
    del <- (a[[p]] - ls * b[[p]]) %*% xs
    dal <- sum(xs * (b[[p]] %*% xs))
    for (r in 1:q) {
      hess[p, r] <- 2 * sum(dx[, r] * del) - dl[r] * dal
    }
  }
  return(hess)
}

hessianCheckEval <- function(theta, a, b, s) {
  myFuncEval <- function(x, ain, bin, ss) {
    q <- length(ain)
    n <- nrow(ain[[1]])
    aa <- makeMatrix(x, ain)
    bb <- makeMatrix(x, bin)
    eig <- myGeigen(aa, bb)
    return(eig$val[ss])
  }
  hess <- hessian(
    myFuncEval,
    x = theta,
    ain = a,
    bin = b,
    ss = s
  )
  return(hess)
}