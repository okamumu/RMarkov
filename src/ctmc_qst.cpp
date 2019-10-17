#include "marlib_Rcpp.h"
using namespace Rcpp;

template <typename T1, typename MatT>
List Cmarkovqst_gs(const T1& Q, NumericVector xi, NumericVector x0, int steps, double rtol, int maxiter, MatT) {
  const int m = marlib::nrow(Q, MatT());
  const int n = marlib::ncol(Q, MatT());
  if (m != n) {
    stop("Matrix Q should be a square matrix.");
  }
  if (n != x0.size() && n != xi.size()) {
    stop("Vectors xi and x0 should be the same dimension of Q.");
  }
  NumericVector x = clone(x0);
  marlib::marlib_params params;
  params.rtol = rtol;
  params.steps = steps;
  params.maxiter = maxiter;
  double gam = marlib::ctmc_qst_gs(Q, xi, x, params, [](marlib::marlib_params){R_CheckUserInterrupt();}, MatT());
  return List::create(
    Named("x")=x,
    Named("gam")=gam,
    Named("convergence")=(params.info==0),
    Named("iter")=params.iter,
    Named("rerror")=params.rerror
  );
}

// [[Rcpp::export]]
List Cmarkovqst_gs(S4 Q, NumericVector xi, NumericVector x0, int steps, double rtol, int maxiter) {
  std::string classname = as<std::string>(Q.attr("class"));
  if (classname == "dgeMatrix") {
    return Cmarkovqst_gs(Q, xi, x0, steps, rtol, maxiter, marlib::DenseMatrixT());
  } else if (classname == "dgCMatrix") {
    return Cmarkovqst_gs(Q, xi, x0, steps, rtol, maxiter, marlib::CSCMatrixT());
  } else {
    stop("markovst_gs requires either dgeMatrix or dgCMatrix class: %s", classname);
  }
}
