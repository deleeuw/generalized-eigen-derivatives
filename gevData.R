da <- outer(ei(1, 4), ei(1, 4))
db <- matrix(0, 4, 4)

a <- list(
  a1 = matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3),
  a2 = matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), 3, 3),
  a3 = matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3),
  a4 = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), 3, 3),
  a5 = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3),
  a6 = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1), 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3)
)

b <- list(
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  matrix(0, 3, 3),
  b1 = matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3),
  b2 = matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), 3, 3),
  b3 = matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3),
  b4 = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), 3, 3),
  b5 = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3),
  b6 = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1), 3, 3)
)

theta <- c(1, 2, 3, 4, 5, 6, -1, -1, -1, 3, 3, 3)
aa <- makeMatrix(theta, a)
bb <- makeMatrix(theta, b)
eab <- myGeigen(aa, bb)

