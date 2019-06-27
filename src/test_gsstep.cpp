#include "marlib_Rcpp.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int C_sparse_nnz(NumericMatrix A) {
  return marlib::sparse_nnz(A);
}

//' @export
// [[Rcpp::export]]
NumericVector C_nrow_dense(NumericMatrix A) {
  int m = marlib::nrow(A, marlib::DenseMatrixT());
  int n = marlib::ncol(A, marlib::DenseMatrixT());
  return NumericVector::create(m, n);
}

//' @export
// [[Rcpp::export]]
NumericVector C_nrow_csr(S4 A) {
  int m = marlib::nrow(A, marlib::CSRMatrixT());
  int n = marlib::ncol(A, marlib::CSRMatrixT());
  return NumericVector::create(m, n);
}

//' @export
// [[Rcpp::export]]
NumericVector C_nrow_csc(S4 A) {
  int m = marlib::nrow(A, marlib::CSCMatrixT());
  int n = marlib::ncol(A, marlib::CSCMatrixT());
  return NumericVector::create(m, n);
}

//' @export
// [[Rcpp::export]]
NumericVector C_nrow_coo(S4 A) {
  int m = marlib::nrow(A, marlib::COOMatrixT());
  int n = marlib::ncol(A, marlib::COOMatrixT());
  return NumericVector::create(m, n);
}

//// CSR

//' @export
// [[Rcpp::export]]
void C_create_csr(NumericMatrix A, NumericVector value, IntegerVector ptr, IntegerVector ind) {
  marlib::dense_to_sparse(A, value, ptr, ind, 0, marlib::CSRMatrixT());
}

//' @export
// [[Rcpp::export]]
void C_sparse_to_dense_csr(S4 A, NumericMatrix B) {
  marlib::sparse_to_dense(A, B, marlib::CSRMatrixT());
}

//// CSC

//' @export
// [[Rcpp::export]]
void C_create_csc(NumericMatrix A, NumericVector value, IntegerVector ptr, IntegerVector ind) {
  marlib::dense_to_sparse(A, value, ptr, ind, 0, marlib::CSCMatrixT());
}

//' @export
// [[Rcpp::export]]
void C_sparse_to_dense_csc(S4 A, NumericMatrix B) {
  marlib::sparse_to_dense(A, B, marlib::CSCMatrixT());
}

//// COO

//' @export
// [[Rcpp::export]]
void C_create_coo(NumericMatrix A, NumericVector value, IntegerVector row, IntegerVector col) {
  marlib::dense_to_sparse(A, value, row, col, 0, marlib::COOMatrixT());
}

//' @export
// [[Rcpp::export]]
void C_sparse_to_dense_coo(S4 A, NumericMatrix B) {
  marlib::sparse_to_dense(A, B, marlib::COOMatrixT());
}

///// diag

//' @export
// [[Rcpp::export]]
NumericVector C_diag_dense(NumericMatrix A) {
  NumericVector x(std::min(A.nrow(), A.ncol()));
  marlib::diag_get(A, x, marlib::DenseMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_diag_csr(S4 A) {
  int m = Rcpp::as<Rcpp::IntegerVector>(A.slot("Dim"))[0];
  int n = Rcpp::as<Rcpp::IntegerVector>(A.slot("Dim"))[1];
  NumericVector x(std::min(m, n));
  marlib::diag_get(A, x, marlib::CSRMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_diag_csc(S4 A) {
  int m = Rcpp::as<Rcpp::IntegerVector>(A.slot("Dim"))[0];
  int n = Rcpp::as<Rcpp::IntegerVector>(A.slot("Dim"))[1];
  NumericVector x(std::min(m, n));
  marlib::diag_get(A, x, marlib::CSCMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_diag_coo(S4 A) {
  int m = Rcpp::as<Rcpp::IntegerVector>(A.slot("Dim"))[0];
  int n = Rcpp::as<Rcpp::IntegerVector>(A.slot("Dim"))[1];
  NumericVector x(std::min(m, n));
  marlib::diag_get(A, x, marlib::COOMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
void C_diag_set_dense(NumericMatrix A, NumericVector x) {
  marlib::diag_set(A, x, marlib::DenseMatrixT());
}

//' @export
// [[Rcpp::export]]
void C_diag_set_csr(S4 A, NumericVector x) {
  marlib::diag_set(A, x, marlib::CSRMatrixT());
}

//' @export
// [[Rcpp::export]]
void C_diag_set_csc(S4 A, NumericVector x) {
  marlib::diag_set(A, x, marlib::CSCMatrixT());
}

//' @export
// [[Rcpp::export]]
void C_diag_set_coo(S4 A, NumericVector x) {
  marlib::diag_set(A, x, marlib::COOMatrixT());
}
