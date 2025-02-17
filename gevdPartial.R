
gradGeigenEval <- function(theta, a, b, ind) {
  nthe <- length(theta)
  nind <- length(ind)
  aa <- makeMatrix(theta, a)
  bb <- makeMatrix(theta, b)
  eig <- myGeigen(aa, bb)
  xmat <- eig$vec
  labd <- eig$val
  dl <- matrix(0, nind, nthe)
  for (s in 1:nind) {
    si <- ind[s]
    xs <- xmat[, si]
    ls <- labd[si]
    for (p in 1:nthe) {
      dl[s, p] <- sum(xs * ((a[[p]] - ls * b[[p]]) %*% xs))
    }
  }
  return(dl)
}



gradGeigenEvec <- function(theta, a, b, ind) {
  nthe <- length(theta)
  nind <- length(ind)
  nobj <- nrow(a[[1]])
  aa <- makeMatrix(theta, a)
  bb <- makeMatrix(theta, b)
  eig <- myGeigen(aa, bb)
  xmat <- eig$vec
  labd <- eig$val
  dx <- array(0, c(nobj, nind, nthe))
  for (s in 1:nind) {
    si <- ind[s]
    xs <- xmat[, si]
    ls <- labd[si]
    for (p in 1:nthe) {
      aa <- crossprod(xmat, (a[[p]] - ls * b[[p]]) %*% xs)
      bb <- 1 / ((labd - ls) + ei(si, nobj)) - ei(si, nobj)
      cc <- aa * bb
      dd <- drop(crossprod(xs, b[[p]] %*% xs))
      dx[, s, p] <- -xmat %*% cc - dd * xs / 2.0
    }
  }
  return(dx)
}


hessianGeigenEval <- function(theta, a, b, ind) {
  nthe <- length(theta)
  nobj <- nrow(a[[1]])
  nind <- length(ind)
  aa <- makeMatrix(theta, a)
  bb <- makeMatrix(theta, b)
  eig <- myGeigen(aa, bb)
  xmat <- eig$vec
  labd <- eig$val
  dl <- gradGeigenEval(theta, a, b, ind)
  dx <- gradGeigenEvec(theta, a, b, ind)
  ddl <- array(0, c(nthe, nthe, nind))
  for (s in 1:nind) {
    si <- ind[s]
    xs <- xmat[, si]
    ls <- labd[si]
    for (p in 1:nthe) {
      del <- (a[[p]] - ls * b[[p]]) %*% xs
      dal <- sum(xs * (b[[p]] %*% xs))
      for (r in 1:nthe) {
        ddl[p, r, s] <- 2 * sum(dx[, s, r] * del) - dl[s, r] * dal
      }
    }
  }
  return(ddl)
}


hessianGeigenEvec <- function(theta, a, b, ind) {
  nthe <- length(theta)
  nobj <- nrow(a[[1]])
  nind <- length(ind)
  aa <- makeMatrix(theta, a)
  bb <- makeMatrix(theta, b)
  eig <- myGeigen(aa, bb)
  xmat <- eig$vec
  labd <- eig$val
  dl <- gradGeigenEval(theta, a, b, 1:nobj)
  dx <- gradGeigenEvec(theta, a, b, 1:nobj)
  ddx <- array(0, c(nthe, nthe, nobj, nind))
  for (k in 1:nobj) {
    for (s in 1:nind) {
      si <- ind[s]
      xs <- xmat[, si]
      for (r in 1:nthe) {
        carr <- a[[r]] - labd[si] * b[[r]]
        for (u in 1:nthe) {
          dxs <- dx[, si, u]
          hess <- 0.0
          for (t in 1:nobj) {
            if (t == si) {
              next
            }
            xt <- xmat[, t]
            dxt <- dx[, t, u]
            dil <- labd[t] - labd[si]
            hess <- hess + dxt[k] * sum(xt * (carr %*% xs)) / dil
            fac <- sum(dxt * (carr %*% xs)) + sum(dxs * (carr %*% xt))
            fac <- fac - dl[si, u] * sum(xs * (b[[r]] %*% xt))
            hess <- hess + (fac / dil) * xt[k]
            fac <- xt[k] * sum(xs * (carr %*% xt))
            fac <- fac * (dl[t, u] - dl[si, u])
            hess <- hess - fac / (dil^2)
          }
          hess <- hess + sum(dxs * (b[[r]] %*% xs)) * xs[k]
          hess <- hess + sum(xs * (b[[r]] %*% xs)) * dxs[k] / 2
          ddx[u, r, k, s] <- -hess
        }
      }
    }
  }
  return(ddx)
}
