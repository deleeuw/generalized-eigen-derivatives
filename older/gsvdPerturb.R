perturbGsvdComp <- function(f, g, h, df, dg, dh, eps) {
  n <- nrow(f)
  m <- ncol(f)
  a <- matrix(0, n + m, n + m)
  a[1:n, n + 1:m] <- f
  a <- a + t(a)
  b <- matrix(0, n + m, n + m)
  b[1:n, 1:n] <- g
  b[n + 1:m, n + 1:m] <- h
  eab <- myGeigen(a, b)
  dl <- rep(0, m)
  for (s in 1:m) {
    xs <- eab$vec[1:n, s]
    ys <- eab$vec[n + 1:m, s]
    ls <- eab$val[s]
    dl[s] <- 2 * sum(xs * (df %*% ys)) - ls * (sum(xs * (dg %*% xs)) + sum(ys * (dh %*% ys)))
    for (t in 1:m) {
      if (t == s) {
        next
      }
    }
  }
  pab <- perturbGsvdCheck(f, g, h, df, dg, dh, eps)
  return(list(eab$val[1:m] + eps * dl, pab$val[1:m]))
}

perturbGsvdCheck <- function(f, g, h, df, dg, dh, eps) {
  n <- nrow(f)
  m <- ncol(f)
  fp <- f + eps * df
  gp <- g + eps * dg
  hp <- h + eps * dh
  a <- matrix(0, n + m, n + m)
  a[1:n, n + 1:m] <- fp
  a <- a + t(a)
  b <- matrix(0, n + m, n + m)
  b[1:n, 1:n] <- gp
  b[n + 1:m, n + 1:m] <- hp
  eab <- myGeigen(a, b)
  return(eab)
}