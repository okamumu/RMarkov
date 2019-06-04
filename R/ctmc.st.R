#' GTH algorithm for the stationary solution of continuous-time Markov chain (CTMC)
#'
#' Compute the stationary probability vector of CTMC with GTH
#' \deqn{x Q = 0, \quad x1 = 1}
#'
#' GTH algorithm is a variation of Guss elimination, and thus this requires n^2 memory.
#'
#' @param Q A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix, dgRMatrix, dgCMatrix and dgTMatrix are allowed.
#' @param nmax The maximum size of the matrix Q.
#' @return Stationary vector
#' @examples
#' Q <- rbind(
#'   c(-10, 0, 3, 7),
#'   c(4, -5, 1, 0),
#'   c(0, 2, -3, 1),
#'   c(0, 0, 1, -1)
#' )
#' ctmc.st.gth(Q)
#' @export

ctmc.st.gth <- function(Q, nmax = 500) {
  if (dim(Q)[1L] > nmax)
    stop("The dimension of Q should be up to nmax.")
  if (any(Matrix::diag(Q) == 0))
    stop("The diagonal elements of Q should not be zero.")

  if (is.matrix(Q)) {
    Q <- as(Q, "dgeMatrix")
  }
  Cmarkovst_gth(Q)
}

#' Power method for the stationary solution of continuous-time Markov chain (CTMC)
#'
#' Compute the stationary probability vector of CTMC with power method
#' \deqn{x Q = 0, \quad x1 = 1}
#'
#' The power method can be applied to obtaining the quasi stationary of CTMC.
#'
#' @param Q A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix, dgRMatrix, dgCMatrix and dgTMatrix are allowed.
#' @param x0 A numeric vector. The inital vector for poiwer method.
#' @param ufact A numeric value for the uniformization. The default is 1.01.
#' @param steps An integer. The algoritm checks the convergence every 'steps'.
#' @param rtol A numeric value. The algorithm stops when the relative error is less than 'rtol'.
#' @param maxiter An integer. The maximum number of iterations.
#' @param matrix.class A string to indicate the matrix type which is used in the computation.
#' @examples
#' Q <- rbind(
#'   c(-10, 0, 3, 7),
#'   c(4, -5, 1, 0),
#'   c(0, 2, -3, 1),
#'   c(0, 0, 1, -1)
#' )
#' ctmc.st.power(Q)
#' @export

ctmc.st.power <- function(Q, x0, ufact = 1.01, maxiter = 5000L, steps = 50L,
                       rtol = sqrt(.Machine$double.eps),
                       matrix.class = NULL) {

  if (missing(x0) || is.null(x0)) {
    n <- dim(Q)[1L]
    x0 <- rep(1/n, n)
  }
  if (!is.null(matrix.class)) {
    Q <- as(Q, matrix.class)
  } else if (is.matrix(Q)) {
    Q <- as(Q, "dgeMatrix")
  }
  Cmarkovst_power(Q, x0=x0, ufact=ufact, steps=steps, rtol=rtol, maxiter=maxiter)
}

#' Gauss-Seidal (GS) algorithm for the stationary solution of continuous-time Markov chain (CTMC)
#'
#' Compute the stationary probability vector of CTMC with GS
#' \deqn{x Q = 0, \quad x1 = 1}
#'
#' In general, GS is an iterative method to solve the linear equation.
#' Since the GS maintains the sparse storage of matrix, it is suitable for solving the large-sized CTMC.
#'
#' TODO: SOR (successive over relaxation) should be implemented.
#'
#' @param Q A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix and dgCMatrix are allowed.
#' @param x0 A numeric vector. The inital vector for GS method.
#' @param steps An integer. The algoritm checks the convergence every 'steps'.
#' @param rtol A numeric value. The algorithm stops when the relative error is less than 'rtol'.
#' @param maxiter An integer. The maximum number of iterations.
#' @param matrix.class A string to indicate the matrix type which is used in the computation.
#' @examples
#' Q <- rbind(
#'   c(-10, 0, 3, 7),
#'   c(4, -5, 1, 0),
#'   c(0, 2, -3, 1),
#'   c(0, 0, 1, -1)
#' )
#' ctmc.st.gs(Q, matrix.class="CsparseMatrix")
#' @export

ctmc.st.gs <- function(Q, x0, maxiter = 5000L, steps = 50L,
                       rtol = sqrt(.Machine$double.eps),
                       matrix.class = NULL) {

  if (missing(x0) || is.null(x0)) {
    n <- dim(Q)[1L]
    x0 <- rep(1/n, n)
  }
  if (!is.null(matrix.class)) {
    Q <- as(Q, matrix.class)
  } else if (is.matrix(Q)) {
    Q <- as(Q, "dgeMatrix")
  }
  Cmarkovst_gs(Q, x0=x0, steps=steps, rtol=rtol, maxiter=maxiter)
}

#
# ctmc.st.sor.cbfunc <- function(iter) {
#   cat('#')
#   if (iter < 100) {
#     steps <- 5
#   } else if (iter < 500) {
#     steps <- 10
#   } else if (iter < 1000) {
#     steps <- 20
#   } else {
#     steps <- 50
#   }
#   steps
# }
#
# ctmc.st.sor.updateomega <- function(iter, dx1, dx2, dx3, x) {
#   # d1 <- max(abs(x1))
#   # d2 <- max(abs(x2))
#   # if (iter >= 2000 && d2 != 0.0) {
#   #   lambda2 <- (d1/d2 + omega - 1.0) / (omega * sqrt(d1/d2))
#   #   if (abs(lambda2) < 1.0) {
#   #     omega <- 2.0 / (1.0 + sqrt(1.0 - lambda2^2))
#   #   } else {
#   #     omega <- 1.0
#   #   }
#   #   cat(gettextf("omega update %f\n", omega))
#   # } else if (omega != 1.0) {
#   #   omega <- 1.0
#   # }
# }
#
# ctmc.st.sor <- function(Q, x0, maxiter = 5000L, steps = 50L, rtol = sqrt(.Machine$double.eps),
#     omega = 1.0, matrix.class = NULL, cbfunc = ctmc.st.sor.cbfunc, update.omega = ctmc.st.sor.updateomega, ...) {
#
#     if (missing(x0) || is.null(x0)) {
#         n <- dim(Q)[1L]
#         x0 <- rep(1/n, n)
#     }
#
#     if (!is.null(matrix.class)) {
#         Q <- as(Q, matrix.class)
#     } else if (is.matrix(Q)) {
#         Q <- as(Q, "dgeMatrix")
#     }
#     stopifnot(all(Matrix::diag(Q) != 0))
# #      Q <- diag.padding.zero(Q)
#
#     res <- switch(class(Q),
#         "dgeMatrix" = .Call(markov_st_sor_dense, Q, x0, maxiter, rtol, omega, steps,
#           body(cbfunc), body(update.omega), environment()),
#         "dgRMatrix" = .Call(markov_st_sor_csr, Q, x0, maxiter, rtol, omega, steps,
#           body(cbfunc), body(update.omega), environment()),
#         "dgCMatrix" = .Call(markov_st_sor_csc, Q, x0, maxiter, rtol, omega, steps,
#           body(cbfunc), body(update.omega), environment()),
#         "dgTMatrix" = {
#             Q <- as(Q, "dgCMatrix")
#             .Call(markov_st_sor_csc, Q, x0, maxiter, rtol, omega, steps,
#               body(cbfunc), body(update.omega), environment())
#         },
#         stop("invalid Matrix class:", class(Q)))
#
#     list(x=res[[1]], convergence=(res[[2]]==0), iter=res[[3]], rerror=res[[4]], omega=res[[5]])
# }
