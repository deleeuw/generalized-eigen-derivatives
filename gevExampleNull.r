k1 <- c(1, 1) / sqrt(2)
k2 <- c(1, -1) / sqrt(2)

e <- diag(4:1)
aa <- k %*% d %*% t(k)
bb <- k %*% e %*% t(k)
aintercept <- diag(2)
bintercept <- diag(2)

theA <- function(theta) {
  return(aintercept + theta[1] * outer(k1, k1) + theta[2] * outer(k2, k2))
}

theB <- function(theta) {
  return(bintercept + theta[1] * outer(k2, k2) + theta[2] * outer(k1, k1))
}

dsA <- function(theta, s) {
  if (s == 1) {
    return(outer(k1, k1))
  } else {
    return(outer(k2, k2))
  }
}

dsB <- function(theta, s) {
  if (s == 1) {
    return(outer(k2, k2))
  } else {
    return(outer(k1, k1))
  }}

dstA <- function(theta, s, t) {
  return(0)
}

dstB <- function(theta, s, t) {
  return(0)
}