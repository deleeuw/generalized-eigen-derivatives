hessian <- function(f, x, eps = 1e-5) {
  n <- length(x)
  fx <- f(x)
  dh <- matrix(0, n, n)
  for (i in 1:n) {
    dif <- f(x + eps * ei(i, n)) - fx
    for (j in 1:n) {
      dijf <- f(x + eps * (ei(i, n) + ei(j, n))) - fx
      djf <- f(x + eps * ei(j, n)) - fx
      dh[i, j] <- (dijf - (dif + djf)) / eps^2
    }
  }
  return(dh)
}

f <- function(x) {
  return((x[1] * x[2]^2 + x[2] * x[1]^2)/2)
}

# hessian(f, c(10, 15))


ei <- function(i, n) {
  return(ifelse(i == 1:n, 1, 0))
}
 