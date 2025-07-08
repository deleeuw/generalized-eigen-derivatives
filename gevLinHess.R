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
theta <- c(thetaa, thetab)

myGeigen <- function(a, b) {
  h <- eigen(solve(b, a))
  lbd <- h$values
  x <- h$vectors
  xbx <- apply(x, 2, function(z)
    sum(z * (b %*% z)))
  x <- x %*% diag(1 / sqrt(xbx))
  return(list(values = lbd, vectors = x))
}

gevLinHess <- function(thetaa,
                       thetab,
                       aarray,
                       barray,
                       aintercept = 0,
                       bintercept = 0) {
  pa <- length(thetaa) # number of parameters for a, s = 1,...,pa
  pb <- length(thetab) # number of parameters for b, s = 1,...,pb
  n <- nrow(aarray[, , 1]) # n is the order of the matrix, i=1,,,n
  p <- pa + pb
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
  h <- myGeigen(a, b)
  x <- h$vectors
  l <- h$values
  ddl <- array(0, c(p, p, n))
  for (i in 1:n) {
    xi <- x[, i]
    li <- l[i]
    mpl <- l - li
    mpl <- diag(ifelse(mpl == 0, 0, 1 / mpl))
    wi <- x %*% mpl %*% t(x)
    for (s in 1:pa) {
      aas <- aarray[, , s]
      xas <- drop(aas %*% xi)
      for (t in 1:pa) {
        aat <- aarray[, , t]
        xat <- drop(aat %*% xi)
        ddl[s, t, i] = -2 * sum(xas * (wi %*% xat))
      }
    }
    for (s in 1:pa) {
      aas <- aarray[, , s]
      xas <- drop(aas %*% xi)
      xxs <- sum(xi * xas)
      for (t in (pa + 1):(pa + pb)) {
        bat <- barray[, , t - pa]
        xbt <- drop(bat %*% xi)
        ddl[s, t, i] = 2 * li * sum(xas * (wi %*% xbt)) - sum(xi * xbt) * xxs
        ddl[t, s, i] = ddl[s, t, i]
      }
    }
    for (s in (pa + 1):(pa + pb)) {
      bas <- barray[, , s - pa]
      xbs <- drop(bas %*% xi)
      xxs <- sum(xi * xbs)
      for (t in (pa + 1):(pa + pb)) {
        bat <- barray[, , t - pa]
        xbt <- drop(bat %*% xi)
        ddl[s, t, i] = -2 * (li ^ 2) * sum(xbs * (wi %*% xbt)) + 2 * li * sum(xi * xbt) * xxs
      }
    }
  }
  return(ddl)
}

checkFun <- function(theta,
                     aarrayl,
                     barrayl,
                     ainterceptl,
                     binterceptl,
                     indil) {
  n <- nrow(aarrayl[, , 1])
  pal <- dim(aarrayl)[3]
  pbl <- dim(barrayl)[3]
  a <- matrix(ainterceptl, n, n)
  for (s in 1:pal) {
    a <- a + theta[s] * aarrayl[, , s]
  }
  b <- matrix(binterceptl, n, n)
  if (pbl > 0) {
    for (s in (pal + 1):(pal + pbl)) {
      b <- b + theta[s] * barrayl[, , s - pal]
    }
  }
  e <- myGeigen(a, b)
  l <- e$values
  return(l[indil])
}

for (indi in 1:4) {
  h <- hessian(
    func = checkFun,
    x = theta,
    aarrayl = aarray,
    barrayl = barray,
    ainterceptl = 0,
    binterceptl = 0,
    indil = indi
  )
  print(h)
}