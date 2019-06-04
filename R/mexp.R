#' Matrix exponential with uniformization
#'
#' This is the function to compute the following matrix exponential with uniformization:
#' \deqn{\exp(At) x}
#'
#'
#' @param A A square matrix. This is one of matrix, dgeMatrix, dgCMatrix, dgRMatrix and dgTMatrix
#' @param x A vector or matrix.
#' @param t A value.
#' @param transpose A logical. If TRUE, the matrix A is transposed. The default is FALSE.
#' @param ufact A numeric value for the uniformization factor.
#' @param eps A numeric value for tolerance error.
#' @param rmax An integer indicating the maximum value of right bound of Poisson distribution. If \code{t} or the absolute
#' diagonal of \code{A} is too large, the right bound may exceed \code{rmax}. This can be relaxed as large as memory allows.
#' @param matrix.class A string indicating the matrix object which is used in the computation. If this is NULL,
#' the matrix class of \code{A} is directly used.
#' @export

mexpAx <- function(A, x, t = 1.0, transpose = FALSE,
                   eps = sqrt(.Machine$double.eps), ufact = 1.01,
                   rmax = 1000, matrix.class = NULL) {
  if (!is.null(matrix.class)) {
    A <- as(A, matrix.class)
  } else if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }

  if (any(Matrix::diag(A) == 0))
    A <- diag.padding.zero(A)

  if (missing(x))
    x <- diag(1, dim(A)[1L])

  if (is.matrix(x)) {
    Cmexp_unif_mat(transpose, A, x, t, ufact, eps, rmax)
  } else if (is.vector(x)) {
    Cmexp_unif_vec(transpose, A, x, t, ufact, eps, rmax)
  } else {
    stop("v should be vector or matrix")
  }
}

#' Integral of matrix exponential with uniformization
#'
#' This is the function to compute the following matrix exponential with uniformization:
#' \deqn{y = \exp(Au) x} and \deqn{cy = \int_0^t \exp(Au) du x}
#'
#'
#' @param A A square matrix. This is one of matrix, dgeMatrix, dgCMatrix, dgRMatrix and dgTMatrix
#' @param x A vector or matrix.
#' @param t A value.
#' @param transpose A logical. If TRUE, the matrix A is transposed. The default is FALSE.
#' @param ufact A numeric value for the uniformization factor.
#' @param eps A numeric value for tolerance error.
#' @param rmax An integer indicating the maximum value of right bound of Poisson distribution. If \code{t} or the absolute
#' diagonal of \code{A} is too large, the right bound may exceed \code{rmax}. This can be relaxed as large as memory allows.
#' @param matrix.class A string indicating the matrix object which is used in the computation. If this is NULL,
#' the matrix class of \code{A} is directly used.
#' @export

cmexpAx <- function(A, x, t = 1.0, transpose = FALSE,
                    eps = sqrt(.Machine$double.eps), ufact = 1.01,
                    rmax = 1000, matrix.class = NULL) {
  if (!is.null(matrix.class)) {
    A <- as(A, matrix.class)
  } else if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }

  if (any(Matrix::diag(A) == 0))
    A <- diag.padding.zero(A)

  if (missing(x))
    x <- diag(1, dim(A)[1L])

  if (is.matrix(x)) {
    Cmexpint_unif_mat(transpose, A, x, t, ufact, eps, rmax)
  } else if (is.vector(x)) {
    Cmexpint_unif_vec(transpose, A, x, t, ufact, eps, rmax)
  } else {
    stop("v should be vector or matrix")
  }
}
