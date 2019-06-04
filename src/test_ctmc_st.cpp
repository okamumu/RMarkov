#include <marlib.h>

using namespace Rcpp;

// general

template <typename T1, typename MatT>
NumericVector C_gth(T1 Q) {
  const int n = marlib::nrow(Q, MatT());
  NumericMatrix A(n,n);
  marlib::dcopy(Q, A, MatT(), marlib::DenseMatrixT());
  NumericVector x(n);
  marlib::gth_impl(A, x);
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_gth_dense(NumericMatrix Q) {
  return C_gth<NumericMatrix, marlib::DenseMatrixT>(Q);
}

//' @export
// [[Rcpp::export]]
NumericVector C_gth_csr(S4 Q) {
  return C_gth<S4, marlib::CSRMatrixT>(Q);
}

//' @export
// [[Rcpp::export]]
NumericVector C_gth_csc(S4 Q) {
  return C_gth<S4, marlib::CSCMatrixT>(Q);
}

//' @export
// [[Rcpp::export]]
NumericVector C_gth_coo(S4 Q) {
  return C_gth<S4, marlib::COOMatrixT>(Q);
}