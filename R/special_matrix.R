#' Special sparce matrices
#'
#' Identity matrix, Zero matrix, A matrix whose elements are filled
#'
#' @name special.matrix
#' @param m An integer. The number of rows.
#' @param n An integer. The number of rows.
#' @param dims A vector for the dimension of matrix.
#' @param v A numeric value. A value to fill the matrix/vector.
#' @return A special matrix/vector.
NULL

#' @describeIn special.matrix An n-by-n identity matrix.
#' @export

eyeM <- function(n) as(sparseMatrix(i=1:n, j=1:n, x=rep.int(x=1, times=n)), "dgCMatrix")

#' @describeIn special.matrix An m-by-n matrix whose elements are filled by v
#' @export

fillM <- function(v, m, n = m, dims = c(m, n)) as(Matrix(v, dims[1], dims[2]), "dgCMatrix")

#' @describeIn special.matrix An m-by-n zero matrix.
#' @export

zeroM <- function(m, n = m, dims = c(m, n)) fillM(0, m=m, n=n)

#' @describeIn special.matrix A vector whose elements are 1.
#' @param byrow A logical, indicating whether the vector is a row vector.
#' @param elem An integer vector, indicating the elements to be fiiled by 1.
#' @export

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
