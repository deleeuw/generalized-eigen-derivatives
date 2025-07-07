
perturbGeigen <- function(a, b, da, db, p) {
  n <- nrow(a)
  m <- length(p)
  dl <- rep(0, m)
  dx <- matrix(0, n, m)
  eig <- myGeigen(a, b)
  xmat <- eig$vec
  labd <- eig$val
  for (s in p) {
    xs <- xmat[, s]
    ls <- labd[s]
    dl[s] <- sum(xs * ((da - ls * db) %*% xs))
    aa <- crossprod(xmat, (da - ls * db) %*% xs)
    bb <- 1 / ((labd - ls) + ei(s, n)) - ei(s, n)
    cc <- aa * bb
    dd <- drop(crossprod(xs, db %*% xs))
    dx[, s] <- -xmat %*% cc - dd * xs / 2.0
  }
  return(list(dl = dl, dx = dx))
}

perturbCheck <- function(a, b, da, db, p, eps) {
  eab <- myGeigen(a, b)
  edd <- myGeigen(a + eps * da, b + eps * db)
  ept <- perturbGeigen(a, b, da, db, p)
  ecp <- edd$val[p]
  eck <- eab$val[p] + eps * ept$dl
  xcp <- edd$vec[,p]
  xck <- eab$vec[,p] + eps * ept$dx
  return(list(ecp = ecp, eck = eck, xcp = xcp, xck = xck))
}
