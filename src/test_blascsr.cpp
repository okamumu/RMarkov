#include <Rcpp.h>
using namespace Rcpp;

#include <marlib.h>

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvN_csr(double alpha, S4 A, NumericVector x,
                      double beta, NumericVector y) {
  marlib::dgemv(marlib::NOTRANS(), alpha, A, x, beta, y, marlib::CSRMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvT_csr(double alpha, S4 A, NumericVector x,
                           double beta, NumericVector y) {
  marlib::dgemv(marlib::TRANS(), alpha, A, x, beta, y, marlib::CSRMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
S4 C_dgerN_csr(double alpha, NumericVector x, NumericVector y, S4 A) {
  marlib::dger(marlib::NOTRANS(), alpha, x, y, A, marlib::CSRMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
S4 C_dgerT_csr(double alpha, NumericVector x, NumericVector y, S4 A) {
  marlib::dger(marlib::TRANS(), alpha, x, y, A, marlib::CSRMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN_csr(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::CSRMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN_csr(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::CSRMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT_csr(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::CSRMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT_csr(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::CSRMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN2_csr(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSRMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN2_csr(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSRMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT2_csr(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSRMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT2_csr(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSRMatrixT());
  return C;
}
