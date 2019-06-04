#include <Rcpp.h>
using namespace Rcpp;

#include <marlib.h>

using namespace Rcpp;

template <typename T1, typename MatT>
NumericVector Cmarkovst_gth(T1 Q, MatT) {
  const int m = marlib::nrow(Q, MatT());
  const int n = marlib::ncol(Q, MatT());
  if (m != n) {
    stop("Matrix Q should be a square matrix.");
  }
  NumericMatrix A(n,n);
  NumericVector x(n);

  dcopy(Q, A, MatT(), marlib::DenseMatrixT());
  marlib::gth_impl(A, x);
  return x;
}

template <typename T1, typename MatT>
List Cmarkovst_power(const T1& Q, NumericVector x0, double ufact, int steps, double rtol, int maxiter, MatT) {
  const int m = marlib::nrow(Q, MatT());
  const int n = marlib::ncol(Q, MatT());
  if (m != n) {
    stop("Matrix Q should be a square matrix.");
  }
  if (n != x0.size()) {
    stop("Vector x0 should be the same dimension of Q.");
  }
  T1 P = clone(Q);
  NumericVector x = clone(x0);
  NumericVector tmpv(n);
  NumericVector prevx(n);
  marlib::unif(P, ufact, tmpv, MatT());
  int iter = 0;
  int info = 1;
  double rerror;

  while(1) {
    marlib::dcopy(x, prevx);
    for (int i=0; i<steps; i++) {
      marlib::dcopy(x, tmpv);
      marlib::dgemv(marlib::TRANS(), 1.0, P, tmpv, 0.0, x, MatT());
      double tmp = marlib::dasum(x);
      marlib::dscal(1.0/tmp, x);
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

template <typename T1, typename MatT>
List Cmarkovst_gs(const T1& Q, NumericVector x0, int steps, double rtol, int maxiter, MatT) {
  const int m = marlib::nrow(Q, MatT());
  const int n = marlib::ncol(Q, MatT());
  if (m != n) {
    stop("Matrix Q should be a square matrix.");
  }
  if (n != x0.size()) {
    stop("Vector x0 should be the same dimension of Q.");
  }
  NumericVector b(n);
  NumericVector x = clone(x0);
  NumericVector prevx(n);
  int iter = 0;
  int info = 1;
  double rerror;

  while(1) {
    marlib::dcopy(x, prevx);
    for (int i=0; i<steps; i++) {
      marlib::gsstep(marlib::TRANS(), 1.0, Q, 0.0, 1.0, b, x, MatT(), marlib::ArrayT());
      double tmp = marlib::dasum(x);
      marlib::dscal(1.0/tmp, x);
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
NumericVector Cmarkovst_gth(S4 Q) {
  std::string classname = as<std::string>(Q.attr("class"));
  if (classname == "dgeMatrix") {
    return Cmarkovst_gth(Q, marlib::DenseMatrixT());
  } else if (classname == "dgRMatrix") {
    return Cmarkovst_gth(Q, marlib::CSRMatrixT());
  } else if (classname == "dgCMatrix") {
    return Cmarkovst_gth(Q, marlib::CSCMatrixT());
  } else if (classname == "dgTMatrix") {
    return Cmarkovst_gth(Q, marlib::COOMatrixT());
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Cmarkovst_power(S4 Q, NumericVector x0, double ufact, int steps, double rtol, int maxiter) {
  std::string classname = as<std::string>(Q.attr("class"));
  if (classname == "dgeMatrix") {
    return Cmarkovst_power(Q, x0, ufact, steps, rtol, maxiter, marlib::DenseMatrixT());
  } else if (classname == "dgRMatrix") {
    return Cmarkovst_power(Q, x0, ufact, steps, rtol, maxiter, marlib::CSRMatrixT());
  } else if (classname == "dgCMatrix") {
    return Cmarkovst_power(Q, x0, ufact, steps, rtol, maxiter, marlib::CSCMatrixT());
  } else if (classname == "dgTMatrix") {
    return Cmarkovst_power(Q, x0, ufact, steps, rtol, maxiter, marlib::COOMatrixT());
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Cmarkovst_gs(S4 Q, NumericVector x0, int steps, double rtol, int maxiter) {
  std::string classname = as<std::string>(Q.attr("class"));
  if (classname == "dgeMatrix") {
    return Cmarkovst_gs(Q, x0, steps, rtol, maxiter, marlib::DenseMatrixT());
  } else if (classname == "dgCMatrix") {
    return Cmarkovst_gs(Q, x0, steps, rtol, maxiter, marlib::CSCMatrixT());
  } else {
    stop("markovst_gs requires either dgeMatrix or dgCMatrix class: %s", classname);
  }
}
