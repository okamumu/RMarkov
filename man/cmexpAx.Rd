% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mexp.R
\name{cmexpAx}
\alias{cmexpAx}
\title{Integral of matrix exponential with uniformization}
\usage{
cmexpAx(A, x, t = 1, cx, transpose = FALSE,
  eps = sqrt(.Machine$double.eps), ufact = 1.01, rmax = 1000,
  matrix.class = NULL)
}
\arguments{
\item{A}{A square matrix. This is one of matrix, dgeMatrix, dgCMatrix, dgRMatrix and dgTMatrix}

\item{x}{A vector or matrix.}

\item{t}{A value.}

\item{cx}{A vector or matrix.}

\item{transpose}{A logical. If TRUE, the matrix A is transposed. The default is FALSE.}

\item{eps}{A numeric value for tolerance error.}

\item{ufact}{A numeric value for the uniformization factor.}

\item{rmax}{An integer indicating the maximum value of right bound of Poisson distribution. If \code{t} or the absolute
diagonal of \code{A} is too large, the right bound may exceed \code{rmax}. This can be relaxed as large as memory allows.}

\item{matrix.class}{A string indicating the matrix object which is used in the computation. If this is NULL,
the matrix class of \code{A} is directly used.}
}
\description{
This is the function to compute the following matrix exponential with uniformization:
\deqn{y = \exp(Au) x} and \deqn{cy = cx + \int_0^t \exp(Au) du x}
}
