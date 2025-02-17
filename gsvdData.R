f <- 1 / outer(1:4, 1:3, "+")
g <- diag(4)
h <- diag(3)
sfgh <- svd(f)
a <- matrix(0, 7, 7)
a[1:4, 5:7] <- f
a <- a + t(a)
b <- diag(7)
eab <- eigen(a)
df <- matrix(1, 4, 3)
dg <- matrix(1, 4, 4)
dh <- matrix(1, 3, 3)

