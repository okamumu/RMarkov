#include "marlib_Rcpp.h"
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

  marlib::ctmc_st_gth(Q, A, x, MatT());
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
  unif(P, ufact, tmpv, MatT());
  marlib::ctmc_st_params params(rtol, steps, maxiter);
  marlib::ctmc_st_power(P, x, params, [](marlib::ctmc_st_params){R_CheckUserInterrupt();}, MatT());
  return List::create(
    Named("x")=x,
    Named("convergence")=(params.info==0),
    Named("iter")=params.iter,
    Named("rerror")=params.rerror
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
  NumericVector x = clone(x0);
  marlib::ctmc_st_params params(rtol, steps, maxiter);
  marlib::ctmc_st_gs(Q, x, params, [](marlib::ctmc_st_params){R_CheckUserInterrupt();}, MatT());
  return List::create(
    Named("x")=x,
    Named("convergence")=(params.info==0),
    Named("iter")=params.iter,
    Named("rerror")=params.rerror
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