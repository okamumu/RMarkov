#include <Rcpp.h>
using namespace Rcpp;

#include <marlib.h>

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvT(double alpha, NumericMatrix A, NumericVector x,
                      double beta, NumericVector y) {
  marlib::dgemv(marlib::TRANS(), alpha, A, x, beta, y, marlib::DenseMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvN(double alpha, NumericMatrix A, NumericVector x,
                       double beta, NumericVector y) {
  marlib::dgemv(marlib::NOTRANS(), alpha, A, x, beta, y, marlib::DenseMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgerT(double alpha, NumericVector x, NumericVector y, NumericMatrix A) {
  marlib::dger(marlib::TRANS(), alpha, x, y, A, marlib::DenseMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgerN(double alpha, NumericVector x, NumericVector y, NumericMatrix A) {
  marlib::dger(marlib::NOTRANS(), alpha, x, y, A, marlib::DenseMatrixT());
  return A;
}
