
set.seed(12345)
theta <- c(.1, .2, .3, .4, .5, .5)
nr <- 4
nc <- 3
n <- 7
p <- 6
hessianl <- FALSE
hessianx <- FALSE


theF <- function(theta) {
  return(outer(theta[1:4], 1:3, "^"))
}

theG <- function(theta) {
  g <- diag(nr) + theta[5] * (1- diag(nr))
  return(g)
}

theH <- function(theta) {
  h <- diag(nc) + theta[6] * (1 -diag(nc))
  return(h)
}

dF <- function(theta, s) {
  df <- matrix(0, nr, nc)
  if (s > 4) {
    return(df)
  } else {
    df[s, ] <- c(1, 2 * theta[s], 3 * theta[s]^2)
    return(df)
  }
}

dG <- function(theta, s) {
  dg <- matrix(0, nr, nr)
  if (s != 5) {
    return(dg)
  } else {
    return(1- diag(nr))
  }
}

dH <- function(theta, s) {
  dh <- matrix(0, nc, nc)
  if (s != 6) {
    return(dh)
  } else {
    return(1 - diag(nc))
  }
}

ddF <- function(theta, s, t) {
  ddf <- matrix(0, nr, nc)
  if ((s != t) || (s > 4) || (t > 4)) {
    return(ddf)
  }
  else {
    ddf[s, ] <- c(0, 2 , 6 * theta[s])
    return(ddf)
  }
}

ddG <- function(theta, s, t) {
  return(matrix(0, nr, nr))
}

ddH <- function(theta, s, t) {
  return(matrix(0, nc, nc))
}

theA <- function(theta) {
  a <- matrix(0, n, n)
  a[1:nr, (nr + 1):n] <- theF(theta)
  return(a + t(a))
}

theB <- function(theta) {
  b <- matrix(0, n, n)
  b[1:nr, 1:nr] <- theG(theta)
  b[(nr + 1):n, (nr + 1):n] <- theH(theta)
  return(b)
}

dA <- function(theta, s) {
  a <- matrix(0, n, n)
  a[1:nr, (nr + 1):n] <- dF(theta, s)
  return(a + t(a))
}

dB <- function(theta, s) {
  b <- matrix(0, n, n)
  b[1:nr, 1:nr] <- dG(theta, s)
  b[(nr + 1):n, (nr + 1):n] <- dH(theta, s)
  return(b)
}

ddA <- function(theta, s, t) {
  a <- matrix(0, n, n)
  a[1:nr, (nr + 1):n] <- ddF(theta, s, t)
  return(a + t(a))
}

ddB <- function(theta, s, t) {
  b <- matrix(0, n, n)
  b[1:nr, 1:nr] <- ddG(theta, s, t)
  b[(nr + 1):n, (nr + 1):n] <- ddH(theta, s, t)
  return(b)
}