#include <Rcpp.h>
using namespace Rcpp;

#include <marlib.h>

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
  NumericVector prevx(n);
  int iter = 0;
  int info = 1;
  double rerror;

  while(1) {
    marlib::dcopy(x, prevx);
    for (int i=0; i<steps; i++) {
      marlib::gsstep(marlib::TRANS(), -1.0, Q, 0.0, 1.0, b, x, MatT(), marlib::ArrayT());
      double tmp = marlib::dsum(x);
      marlib::daxpy(-tmp, pis, x);
    }
    marlib::daxpy(-1.0, x, prevx);
    iter += steps;
    rerror = Rcpp::max(Rcpp::abs(prevx/x));
    if (rerror < rtol) {
      info = 0;
      break;
    }
    if (iter >= maxiter) {
      info = -1;
      break;
    }
    R_CheckUserInterrupt();
  }
  return List::create(
    Named("x")=x,
    Named("convergence")=(info==0),
    Named("iter")=iter,
    Named("rerror")=rerror
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
