

set.seed(12345)
aarray <- array(0, c(4, 4, 3))
barray <- array(0, c(4, 4, 3))
aarray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aarray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 1] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 2] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
barray[, , 3] <- crossprod(matrix(rnorm(400), 100, 4)) / 100
aintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
bintercept <- crossprod(matrix(rnorm(400), 100, 4)) / 100
thetaa <- c(1, 1, 1)
thetab <- c(1, 1, 1)

myGeigen <- function(a, b) {
  h <- eigen(solve(b, a))
  lbd <- h$values
  x <- h$vectors
  xbx <- apply(x, 2, function(z)
    sum(z * (b %*% z)))
  x <- x %*% diag(1 / sqrt(xbx))
  return(list(values = lbd, vectors = x))
}

gevLinGrad <- function(thetaa,
                       thetab,
                       aarray,
                       barray,
                       aintercept = 0,
                       bintercept = 0) {
  pa <- length(thetaa) # number of parameters for a, s = 1,...,pa
  pb <- length(thetab) # number of parameters for b, s = 1,...,pb
  n <- nrow(aarray[, , 1]) # n is the order of the matrix, i=1,,,n
  a <- matrix(aintercept, n, n)
  for (s in 1:pa) {
    a <- a + thetaa[s] * aarray[, , s]
  }
  b <- matrix(bintercept, n, n)
  if (pb > 0) {
    for (s in 1:pb) {
      b <- b + thetab[s] * barray[, , s]
    }
  }
  e <- myGeigen(a, b)
  x <- e$vectors
  l <- e$values
  dl <- matrix(0, n, pa + pb)
  dx <- array(0, c(n, pa + pb, n))
  for (i in 1:n) {
    # eigen loop, for each i dl_i is a p-vector and dx_i is an n * p matrix
    xi <- x[, i]
    li <- l[i]
    mpl <- l - li
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    w <- x %*% mpl %*% t(x)
    for (s in 1:pa) {
      dfar <- aarray[, , s]
      dl[i, s] <- sum(xi * (dfar %*% xi))
      dx[, s, i] <- -w %*% dfar %*% xi
    }
    for (s in (pa + 1):(pa + pb)) {
      sp <- s - pa
      bsp <- barray[, , sp]
      dfar <- -li * bsp
      dl[i, s] <- sum(xi * (dfar %*% xi))
      dx[, s, i] <- -w %*% dfar %*% xi - 0.5 * sum(xi * (bsp %*% xi)) * xi
    }
  }
  return(list(
    values = l,
    vectors = x,
    dvalues = dl,
    dvectors = dx
  ))
}


