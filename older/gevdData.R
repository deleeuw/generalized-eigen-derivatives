source("gevdUtils.R")

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

ar <- list(
  ar1 = matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3),
  ar2 = matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), 3, 3),
  ar3 = matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3),
  ar4 = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), 3, 3),
  ar5 = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3),
  ar6 = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1), 3, 3)) 

br <- list(
  br1 = matrix(0, 3, 3),
  br2 = matrix(0, 3, 3),
  br3 = matrix(0, 3, 3),
  br4 = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), 3, 3),
  br5 = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3),
  br6 = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1), 3, 3))

tr <- c(1, 2, 3, 4, 5, 6)

as <- list(
  as1 = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3),
  as2 = matrix(c(0, -1, -1, -1, 0, -1, -1, -1, 0), 3, 3),
  as3 = matrix(0, 3, 3),
  as4 = matrix(0, 3, 3))

bs <- list(
  bs1 = matrix(0, 3, 3),
  bs2 = matrix(0, 3, 3),
  bs3 = matrix(c(2, -1, 0, -1, 2, 0, 0, 0, 0), 3, 3),
  bs4 = matrix(c(0, 0, 0, 0, 2, -1, 0, -1 ,2), 3, 3))

ts <- c(3, 1, 1, 1) 