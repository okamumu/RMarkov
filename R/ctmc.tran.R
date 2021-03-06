#' Uniformization for transient analysis of countinuous-time Markov chains (CTMCs)
#'
#' Compute the state probability vectors and the cumulative values of them at arbitrary time point.
#'
#' The state probability vector is given by \deqn{x(t) = \exp(Q' t) x} and the cumulative value becomes
#' \deqn{cx(t) = \int_0^t \exp(Q' u) du x}, where the prime indicates the transpose operator.
#'
#' \code{ctmc.tran.unif} uses the uniformization to compute the transient solution.
#'
#' \code{rmax} is the maximum number of rightbound of Poisson distribution. The default is 1000.
#' The required value depends on the length of time interval and the magnitude of diagonal element of the
#' infinitesimal generator Q. If the program stops when the required rightbound exceeds \code{rmax},
#' it can be changed so that \code{rmax} exceeds the requirement value of rightbound.
#'
#' \code{ctmc.tran.unif} gives the state probability vectors and the cumulative values of them at every time points.
#'
#' @param Q A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix, dgRMatrix, dgCMatrix and dgTMatrix are allowed.
#' @param x An initial probability vector or a matrix for initial probability vectors.
#' @param t A nemeric vector, indicating the time sequence, which should be a non-increasing.
#' @param cx An initial probability vector or a matrix for initial cumulative vectors.
#' @param transpose A logical, indicating the matrix Q is transposed or not. The default is TRUE.
#' @param ufact A numeric value for the uniformization. The default is 1.01.
#' @param rmax An integer, indicating the maximum number of rightbound of Poisson distribution in uniformization (see details).
#' @param eps A numeric value of tolarence error in both uniformization, which determines the rightboud of Poisson distribution.
#' @param matrix.class A string character, indicating the class of Matrix for the computation.
#'
#' @return The function returs the following list.
#' \item{t}{A numeric vector of time points.}
#' \item{x}{A numeric vector of state probability vector at every time point.}
#' \item{cx}{A numeric vector of cumulative state probability vector at every time point.}
#'
#' @examples
#' A <- rbind(
#'  c(-2,2,0),
#'  c(3,-5,2),
#'  c(0,1,-1))
#' ctmc.tran.unif(Q=A, x=c(1,0,0), t=seq(0,1,length.out=10))
#'
#' @export

ctmc.tran.unif <- function(Q, x, t, cx, transpose = TRUE,
                           eps = sqrt(.Machine$double.eps), ufact = 1.01,
                           rmax = 1000, matrix.class = NULL) {
  if (!is.null(matrix.class)) {
    Q <- as(Q, matrix.class)
  } else if (is.matrix(Q)) {
    Q <- as(Q, "dgeMatrix")
  }

  if (any(Matrix::diag(Q) == 0))
    Q <- diag.padding.zero(Q)

  dt <- diff(c(0,t))
  if (any(dt < 0))
    stop("Time (argument t) should be a non-decreasing sequence.")

  if (is.matrix(x)) {
    if (missing(cx)) {
      cx <- matrix(0.0, nrow(x), ncol(x))
    }
    res <- Ctran_unif_mat(transpose, Q, x, cx, dt, ufact, eps, rmax)
    list(t=t, x=array(res$x, dim=c(nrow(x), ncol(x), length(t))), cx=array(res$cx, dim=c(nrow(x), ncol(x), length(t))))
  } else {
    if (missing(cx)) {
      cx <- numeric(length(x))
    }
    res <- Ctran_unif_vec(transpose, Q, x, cx, dt, ufact, eps, rmax)
    list(t=t, x=array(res$x, dim=c(length(x), length(t))), cx=array(res$cx, dim=c(length(x), length(t))))
  }
}

#' Uniformization for transient analysis of countinuous-time Markov chains (CTMCs) with rewards
#'
#' Compute the expected instantaneous and cumulative rewards at arbitrary time point.
#'
#' The instantaneous reward is given by \deqn{irwd(t) = r \exp(Q' t) x} and the cumulative reward becomes
#' \deqn{crwd(t) = r \int_0^t \exp(Q' u) du x}, where the prime indicates the transpose operator, x and r are
#' initial probabilities and reward vector or matrix.
#'
#' \code{ctmc.tran.rwd.unif} uses the uniformization to compute the transient solution.
#'
#' \code{rmax} is the maximum number of rightbound of Poisson distribution. The default is 1000.
#' The required value depends on the length of time interval and the magnitude of diagonal element of the
#' infinitesimal generator Q. If the program stops when the required rightbound exceeds \code{rmax}, it can be
#' changed so that \code{rmax} exceeds the requirement value of rightbound.
#'
#' \code{transient_unif} and \code{transient_rwd_unif} are written by Rcpp, which are called from \code{ctmc.tran.unif}.
#'
#' If the argument \code{r} is missing, \code{ctmc.tran.unif} gives the state probability vectors and the cumulative values of them
#' at every time points. If the argument \code{r} is given, the return value provides the state probability vector and the cumulative one
#' at the last time points. In addition, the instantaneous and cumulative rewards (\code{irwd} and \code{crwd}) at every time points are provided.
#'
#' The argument \code{direction} corresponds to the transpose operator of Q. When \code{x} is a matrix and \code{direction} is \code{forward},
#' the state probability vector is computed as \deqn{x(t) = \exp(t(Q) t) x.} Then the instantaneous reward becomes \deqn{irwd(t) = t(r) \exp(t(Q) t) x.}
#' When \code{x} is a matrix and \code{direction} is \code{backward},
#' the state probability vector is computed as \deqn{x(t) = \exp(Q t) x.} Then the instantaneous reward becomes \deqn{irwd(t) = t(r) \exp(Q t) x.}
#' Therefore, it should be noted that the roles of \code{x} and \code{r} are opposite when \code{direction} is \code{backward}.
#'
#' @param Q A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix, dgRMatrix, dgCMatrix and dgTMatrix are allowed.
#' @param x An initial probability vector or a matrix for initial probability vectors.
#' @param t A nemeric vector, indicating the time sequence, which should be a non-increasing.
#' @param r A numeric vector or matrix indicating the reward.
#' @param cx An initial cumulative vector or a matrix.
#' @param transpose A logical, indicating the matrix Q is transposed or not. The default is TRUE.
#' @param ufact A numeric value for the uniformization. The default is 1.01.
#' @param rmax An integer, indicating the maximum number of rightbound of Poisson distribution in uniformization (see details).
#' @param eps A numeric value of tolarence error in both uniformization, which determines the rightboud of Poisson distribution.
#' @param matrix.class A string character, indicating the class of Matrix for the computation.
#'
#' @return The function returs the following list.
#' \item{t}{A numeric vector of time points.}
#' \item{x}{A numeric vector of state probability vector at the last time point.}
#' \item{cx}{A numeric vector of cumulative state probability vector at the last time point.}
#' \item{irwd}{A numeric vector of instantaneous rewards at every time point.}
#' \item{crwd}{A numeric vector of cumulative rewards at every time point.}
#'
#' @examples
#' A <- rbind(
#'  c(-2,2,0),
#'  c(3,-5,2),
#'  c(0,1,-1))
#' ctmc.tran.rwd.unif(Q=A, x=c(1,0,0), r=c(1,1,0), t=seq(0,1,length.out=10))
#'
#' @export

ctmc.tran.rwd.unif <- function(Q, x, r, t, cx, transpose = TRUE,
                               eps = sqrt(.Machine$double.eps), ufact = 1.01,
                               rmax = 1000, matrix.class = NULL) {
  if (!is.null(matrix.class)) {
    Q <- as(Q, matrix.class)
  } else if (is.matrix(Q)) {
    Q <- as(Q, "dgeMatrix")
  }

  if (any(Matrix::diag(Q) == 0))
    Q <- diag.padding.zero(Q)

  dt <- diff(c(0,t))
  if (any(dt < 0))
    stop("Time (argument t) should be a non-decreasing sequence.")

  if (is.matrix(x) && is.matrix(r)) {
    if (missing(cx)) {
      cx <- matrix(0, nrow(x), ncol(x))
    }
    res <- Ctran_unif_rwd_mat_mat(transpose, Q, x, cx, r, dt, ufact, eps, rmax)
    list(t=t, x=res$x, cx=res$cx,
         irwd=array(res$irwd, dim=c(ncol(r), ncol(x), length(t))),
         crwd=array(res$crwd, dim=c(ncol(r), ncol(x), length(t))))
  } else if (is.vector(x) && is.matrix(r)) {
    if (missing(cx)) {
      cx <- numeric(length(x))
    }
    res <- Ctran_unif_rwd_vec_mat(transpose, Q, x, cx, r, dt, ufact, eps, rmax)
    list(t=t, x=res$x, cx=res$cx,
         irwd=array(res$irwd, dim=c(ncol(r), length(t))),
         crwd=array(res$crwd, dim=c(ncol(r), length(t))))
  } else if (is.vector(x) && is.vector(r)) {
    if (missing(cx)) {
      cx <- numeric(length(x))
    }
    res <- Ctran_unif_rwd_vec_vec(transpose, Q, x, cx, r, dt, ufact, eps, rmax)
    list(t=t, x=res$x, cx=res$cx,
         irwd=res$irwd,
         crwd=res$crwd)
  } else if (is.matrix(x) && is.vector(r)) {
    if (missing(cx)) {
      cx <- matrix(0, nrow(x), ncol(x))
    }
    r <- matrix(r, ncol=1)
    res <- Ctran_unif_rwd_mat_mat(transpose, Q, x, cx, r, dt, ufact, eps, rmax)
    list(t=t, x=res$x, cx=res$cx,
         irwd=array(res$irwd, dim=c(ncol(r), ncol(x), length(t))),
         crwd=array(res$crwd, dim=c(ncol(r), ncol(x), length(t))))
  }
}
