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
