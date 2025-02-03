

myGeigen <- function(a, b) {
  e <- eigen(solve(b, a))
  val <- e$values
  vec <- e$vectors
  d <- diag(crossprod(vec, b %*% vec))
  vec <- vec %*% diag(1 / sqrt(d))
  return(list(val = val, vec = vec))
}

ei <- function(i, n)
  ifelse(i == 1:n, 1, 0)

makeMatrix <- function(x, a) {
  m <- length(x)
  n <- nrow(a[[1]])
  aa <- matrix(0, n , n)
  for (k in 1:m) {
    aa <- aa + x[k] * a[[k]]
  }
  return(aa)
}

matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}
