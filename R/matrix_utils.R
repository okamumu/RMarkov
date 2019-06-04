eyeM <- function(n) Matrix::sparseMatrix(i=1:n, j=1:n, x=rep.int(x=1, times=n))

zeroM <- function(m, n = m, dims = c(m, n)) Matrix(0, dims[1], dims[2])

fillM <- function(v, m, n = m, dims = c(m, n)) Matrix(v, dims[1], dims[2])

vone <- function(n, byrow = FALSE, elem=1:n) {
  if(byrow) {
    v <- Matrix(0, 1, n)
    v[elem] <- 1
    v
  } else {
    v <- Matrix(0, n, 1)
    v[elem] <- 1
    v
  }
}
