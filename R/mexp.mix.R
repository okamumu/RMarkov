#' Integral of a matrix exponential with a density function
#'
#' Compute a mixture of matrix exponential with a density function.
#'
#' \code{mexpAx.mix} computes the mixture of matrix exponential of \code{A} with a density function \code{f},
#' multiplied by a vector or matrix \code{x}, i.e., \eqn{\int \exp(A u) f(u) du x}.
#' \code{mexpAx.mix} uses the uniformization for the computation of the matrix exponential.
#' Also, the default numerical integration is the doubly-exponential formulas (DE formulas),
#' which is provided by the package deformula.
#'
#' The argument \code{inte.func} provides the points and weights for numerical integration. This function
#' determines the domain of integration as well as the numerical quadrature. The default is \code{deformula::deformula.zeroinf}
#' that is the deformula on the domain [0, infinity). The function of \code{inte.func} should returns a list with the elements
#' \code{x}, \code{w} and \code{h}. The element \code{x} indicates a vector of points for integration which should be a non-decreasing sequence.
#' \code{w} is a vector of weights corresponding the points. \code{h} is an additional weight. That is, the exact integration should be computed
#' \eqn{h \sum_i w_i f(x_i)}.
#'
#' @param A A square matrix.
#' @param x A numeric vector or matrix.
#' @param f A density function. The parameters of \code{f} are given in the argument \code{...}.
#' @param ... A list for parameters of \code{f}.
#' @param transpose A logical; if TRUE, the matrix A is transposed.
#' @param ufact A numeric value for the uniformization.
#' @param rmax An integer, indicating the maximum number of rightbound of poisson distribution in uniformization.
#' @param eps A numeric value of tolarence error in uniformization.
#' @param matrix.class A string character, indicating the class of Matrix for the computation.
#' @param inte.func A function that provides the points and weights for numerical integration.
#' @return A vector or matrix of integral of matrix exponential.
#' @examples
#' ### integration with Gamma distribution with shape=2 and scale=1 on domain [0, infinity)
#' A <- rbind(
#'  c(-2,2,0),
#'  c(3,-5,2),
#'  c(0,1,-1))
#' mexpAx.mix(A=A, f=dgamma, shape=2, scale=1)
#' @export

mexpAx.mix <- function(A, x, f, ...,
                       transpose = FALSE,
                       eps = sqrt(.Machine$double.eps),
                       ufact = 1.01,
                       rmax = 1000,
                       matrix.class = NULL,
                       inte.func = deformula::deformula.zeroinf) {
  ires <- inte.func(f=f, ...)

  if (!is.null(matrix.class)) {
    A <- as(A, matrix.class)
  } else if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }

  if (any(Matrix::diag(A) == 0))
    A <- diag.padding.zero(A)

  if (missing(x))
    x <- diag(1, dim(A)[1L])

  dt <- diff(c(0,ires$x))

  if (is.matrix(x)) {
    Cmexp_mix_unif_mat(transpose, A, x, ires$w, dt, ufact, eps, rmax) * ires$h
  } else if (is.vector(x)) {
    Cmexp_mix_unif_vec(transpose, A, x, ires$w, dt, ufact, eps, rmax) * ires$h
  } else {
    stop("v should be vector or matrix")
  }
}

#' Integral of a matrix exponential with a density function
#'
#' Compute a mixture of matrix exponential with a density function
#'
#' \code{cmexpAx.mix} computes a cumulative value \eqn{\int \int_0^z exp(transpose(A) u) du f(z) dz}.
#'
#' \code{cmexpAx.mix} uses the uniformization for the computation of the matrix exponential.
#' Also, the default numerical integration is the doubly-exponential formulas (DE formulas),
#' which is provided by the package deformula.
#'
#' The argument \code{inte.func} provides the points and weights for numerical integration. This function
#' determines the domain of integration as well as the numerical quadrature. The default is \code{deformula::deformula.zeroinf}
#' that is the deformula on the domain [0, infinity). The function of \code{inte.func} should returns a list with the elements
#' \code{x}, \code{w} and \code{h}. The element \code{x} indicates a vector of points for integration which should be a non-decreasing sequence.
#' \code{w} is a vector of weights corresponding the points. \code{h} is an additional weight. That is, the exact integration should be computed
#' \eqn{h \sum_i w_i f(x_i)}.
#'
#' @param A A square matrix.
#' @param x A numeric vector or matrix.
#' @param f A density function. The parameters of \code{f} are given in the argument \code{...}.
#' @param ... A list for parameters of \code{f}.
#' @param transpose A logical; if TRUE, the matrix A is transposed.
#' @param ufact A numeric value for the uniformization.
#' @param rmax An integer, indicating the maximum number of rightbound of poisson distribution in uniformization.
#' @param eps A numeric value of tolarence error in uniformization.
#' @param matrix.class A string character, indicating the class of Matrix for the computation.
#' @param inte.func A function that provides the points and weights for numerical integration.
#' @return A list with components
#' \item{y}{\eqn{\int \exp(A u) f(u) du x}.}
#' \item{cy}{\eqn{\int \int_0^z \exp(A u) du f(z) dz x}}
#'
#' @examples
#' ### integration with Gamma distribution with shape=2 and scale=1 on domain [0, infinity)
#' A <- rbind(
#'  c(-2,2,0),
#'  c(3,-5,2),
#'  c(0,1,-1))
#' cmexpAv.mix(A=A, f=dgamma, shape=2, scale=1)
#' @export

cmexpAx.mix <- function(A, x, f, ...,
                        transpose = FALSE,
                        eps = sqrt(.Machine$double.eps),
                        ufact = 1.01,
                        rmax = 1000,
                        matrix.class = NULL,
                        inte.func = deformula::deformula.zeroinf) {
  ires <- inte.func(f=f, ...)

  if (!is.null(matrix.class)) {
    A <- as(A, matrix.class)
  } else if (is.matrix(A)) {
    A <- as(A, "dgeMatrix")
  }

  if (any(Matrix::diag(A) == 0))
    A <- diag.padding.zero(A)

  if (missing(x))
    x <- diag(1, dim(A)[1L])

  dt <- diff(c(0,ires$x))

  if (is.matrix(x)) {
    res <- Cmexpint_mix_unif_mat(transpose, A, x, ires$w, dt, ufact, eps, rmax)
  } else if (is.vector(x)) {
    res <- Cmexpint_mix_unif_vec(transpose, A, x, ires$w, dt, ufact, eps, rmax)
  } else {
    stop("v should be vector or matrix")
  }
  lapply(res, function(elem) elem * ires$h)
}
