% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmc.st.R
\name{ctmc.st.power}
\alias{ctmc.st.power}
\title{Power method for the stationary solution of continuous-time Markov chain (CTMC)}
\usage{
ctmc.st.power(Q, x0, ufact = 1.01, maxiter = 5000L, steps = 50L,
  rtol = sqrt(.Machine$double.eps), matrix.class = NULL)
}
\arguments{
\item{Q}{A n-by-n matrix. A CTMC kernel. matrix, dgeMatrix, dgRMatrix, dgCMatrix and dgTMatrix are allowed.}

\item{x0}{A numeric vector. The inital vector for poiwer method.}

\item{ufact}{A numeric value for the uniformization. The default is 1.01.}

\item{maxiter}{An integer. The maximum number of iterations.}

\item{steps}{An integer. The algoritm checks the convergence every 'steps'.}

\item{rtol}{A numeric value. The algorithm stops when the relative error is less than 'rtol'.}

\item{matrix.class}{A string to indicate the matrix type which is used in the computation.}
}
\description{
Compute the stationary probability vector of CTMC with power method
\deqn{x Q = 0, \quad x1 = 1}
}
\details{
The power method can be applied to obtaining the quasi stationary of CTMC.
}
\examples{
Q <- rbind(
  c(-10, 0, 3, 7),
  c(4, -5, 1, 0),
  c(0, 2, -3, 1),
  c(0, 0, 1, -1)
)
ctmc.st.power(Q)
}
