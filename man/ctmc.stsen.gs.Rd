% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmc.stsen.R
\name{ctmc.stsen.gs}
\alias{ctmc.stsen.gs}
\title{Sensitivity function of continuous-time Markov chain}
\usage{
ctmc.stsen.gs(Q, pis, b, sinit, maxiter = 5000L, steps = 10L,
  rtol = sqrt(.Machine$double.eps), matrix.class = NULL)
}
\arguments{
\item{Q}{An infinitesimal generator of CTMC. A square matrix.}

\item{pis}{A vector for the stationary distribution.}

\item{b}{A vector which is obtained from the derivatives of Q.
For example, in the case of the first derivative, b = pis * dQ where dQ is the first derivateive of Q
with respect to a parameter.}

\item{sinit}{A numeric vector which is used as an initial vector for the iterative algorithm.}

\item{maxiter}{An integer, indicating the maximum number of itrations.}

\item{steps}{An integer, indicating the number of itrations for checking the convergence.}

\item{rtol}{A numeric value, indicating the tolerance relative error to stop the iterative algorithm.}

\item{matrix.class}{A string character, indicating the class of Matrix for the computation.}
}
\value{
A list with components
\item{x}{The sensitivity function.}
\item{convergence}{A logical, indicating the algorithm converges or not.}
\item{iter}{The number of iterations.}
\item{rerror}{The relative error when the algorithm is stopped.}
}
\description{
Compute the sensitivity function, the first or higher derivative of stationary distribution in terms of
a parameter \eqn{\theta}. For example, the first derivative of stationary distribution \eqn{s} of
continuous-time Markov chain with the infinitesimal generator \eqn{Q} such that
\deqn{s Q + \pi_s \frac{d}{d\theta} Q = 0, \, s 1 = 0.}
}
\details{
\code{ctmc.stsen.gs} uses Gauss-Seidel algorithm to solve it.
This is appropriate in the case of a large sparse Markov kernel.
}
\examples{
A <- rbind(
 c(-2,2,0),
 c(3,-5,2),
 c(0,1,-1))
Adash <- rbind(
 c(-1,1,0),
 c(0,0,0),
 c(0,0,0))
pis <- ctmc.st.gth(A)
ctmc.stsen.gs(A, pis, pis \%*\% Adash)
}
