#' Gauss-Seidal (GS) algorithm for the quasi-stationary solution of continuous-time Markov chain (CTMC)
#'
#' Compute the quasi-stationary probability vector of CTMC with GS
#' \deqn{x Q = -gam x, \quad x1 = 1}
#'
#' In general, GS is an iterative method to solve the linear equation.
#' Since the GS maintains the sparse storage of matrix, it is suitable for solving the large-sized CTMC.
#'
#' @param Q A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix and dgCMatrix are allowed.
#' @param x0 A numeric vector. The inital vector for GS method.
#' @param steps An integer. The algoritm checks the convergence every 'steps'.
#' @param rtol A numeric value. The algorithm stops when the relative error is less than 'rtol'.
#' @param maxiter An integer. The maximum number of iterations.
#' @param matrix.class A string to indicate the matrix type which is used in the computation.
#' @examples
#' Q <- rbind(
#'   c(-15, 0, 3, 7),
#'   c(4, -5, 1, 0),
#'   c(0, 2, -3, 1),
#'   c(0, 0, 1, -1)
#' )
#' ctmc.qst.gs(Q, matrix.class="CsparseMatrix")
#' @export

ctmc.qst.gs <- function(Q, x0, maxiter = 5000L, steps = 50L,
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
  xi <- -as.vector(Q %*% vone(dim(Q)[1L]))
  Cmarkovqst_gs(Q, xi=xi, x0=x0, steps=steps, rtol=rtol, maxiter=maxiter)
}

