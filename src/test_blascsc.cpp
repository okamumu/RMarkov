#include "marlib_Rcpp.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvN_csc(double alpha, S4 A, NumericVector x,
                      double beta, NumericVector y) {
  marlib::dgemv(marlib::NOTRANS(), alpha, A, x, beta, y, marlib::CSCMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvT_csc(double alpha, S4 A, NumericVector x,
                           double beta, NumericVector y) {
  marlib::dgemv(marlib::TRANS(), alpha, A, x, beta, y, marlib::CSCMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
S4 C_dgerN_csc(double alpha, NumericVector x, NumericVector y, S4 A) {
  marlib::dger(marlib::NOTRANS(), alpha, x, y, A, marlib::CSCMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
S4 C_dgerT_csc(double alpha, NumericVector x, NumericVector y, S4 A) {
  marlib::dger(marlib::TRANS(), alpha, x, y, A, marlib::CSCMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN_csc(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::CSCMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN_csc(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::CSCMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT_csc(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::CSCMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT_csc(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::CSCMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN2_csc(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSCMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN2_csc(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSCMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT2_csc(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSCMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT2_csc(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::CSCMatrixT());
  return C;
}
