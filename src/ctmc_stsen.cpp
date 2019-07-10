#include "marlib_Rcpp.h"
using namespace Rcpp;

/**
 * @param y A vector. Ex. in the case of first derivative, b = pi dQ/dtheta
 */

template <typename T1, typename MatT>
List Cmarkovstsen_gs(T1 Q, NumericVector x0, NumericVector b, NumericVector pis,
                     int steps, double rtol, int maxiter, MatT) {
  const int m = marlib::nrow(Q, MatT());
  const int n = marlib::ncol(Q, MatT());
  if (m != n) {
    stop("Matrix Q should be a square matrix.");
  }
  if (n != x0.size()) {
    stop("Vector x0 should be the same dimension of Q.");
  }
  NumericVector x = clone(x0);
  marlib::marlib_params params;
  params.rtol = rtol;
  params.steps = steps;
  params.maxiter = maxiter;
  marlib::ctmc_stsen_gs(Q, x, b, pis, params, [](marlib::marlib_params){R_CheckUserInterrupt();}, MatT());
  return List::create(
    Named("x")=x,
    Named("convergence")=(params.info==0),
    Named("iter")=params.iter,
    Named("rerror")=params.rerror
  );
}

// [[Rcpp::export]]
List Cmarkovstsen_gs(S4 Q, NumericVector x0, NumericVector y,
                     NumericVector pis, int steps, double rtol, int maxiter) {
  std::string classname = as<std::string>(Q.attr("class"));
  if (classname == "dgeMatrix") {
    return Cmarkovstsen_gs(Q, x0, y, pis, steps, rtol, maxiter, marlib::DenseMatrixT());
  } else if (classname == "dgCMatrix") {
    return Cmarkovstsen_gs(Q, x0, y, pis, steps, rtol, maxiter, marlib::CSCMatrixT());
  } else {
    stop("markovstsen_gs requires either dgeMatrix or dgCMatrix class: %s", classname);
  }
}
