set.seed(12345)
data(cattell, package = "psych")
cmat <- cattell
a <- .5 * cbind(matrix(rnorm(24), 12, 2), diag(12))
pattern <- cbind(matrix(1, 12, 2), diag(12))
library(numDeriv)


mdfaDerivatives <- function(a, cmat, pattern) {
  n <- nrow(cmat)
  aca <- crossprod(a, cmat %*% a)
  hca <- eigen(aca)
  hvec <- hca$vectors[, 1:n]
  hval <- diag(1 / sqrt(hca$values[1:n]))
  ac2 <- hvec %*% hval %*% t(hvec)
  g <- pattern * (cmat %*% a %*% ac2)
  return(g)
}

mdfaNumDerivatives <- function(a, cmat) {
  p <- ncol(a) - nrow(a)
  aa <- c(as.vector(a[, 1:p]), diag(a[, -(1:p)]))
  n <- nrow(cmat)
  p <- (length(aa) / n) - 1
  theFunc <- function(aa) {
    aa <- cbind(matrix(aa[1:(n * p)], n , p), diag(aa[-(1:(n * p))]))
    ad <- crossprod(aa, cmat %*% aa)
    lb <- abs(eigen(ad)$values)
    return(sum(sqrt(lb)))
  }
  g <- grad(theFunc, aa)
  h <- hessian(theFunc, aa)
  g <- cbind(matrix(g[1:(n * p)], n , p), diag(g[-(1:(n * p))]))
  return(list(g = g, h = h))
}


