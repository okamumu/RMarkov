#include <Rcpp.h>
using namespace Rcpp;

#include "blas/blas1.h"
#include "blas/blas2.h"
#include "blas/blas3.h"

#include "blas/sparse.h"
#include "blas/blas_csr.h"

#include "type/sparse_matrix.h"

// [[Rcpp::export]]
double C_ddot(NumericVector x, NumericVector y) {
  return marlib::ddot(x,y);
}

// [[Rcpp::export]]
double C_ddot2(S4 x, S4 y) {
  return marlib::ddot(x,y);
}

// [[Rcpp::export]]
double C_dasum(NumericVector x) {
  return marlib::dasum(x);
}

// [[Rcpp::export]]
double C_dasum2(S4 x) {
  return marlib::dasum(x);
}

// [[Rcpp::export]]
double C_idamax(NumericVector x) {
  return marlib::idamax(x);
}

// [[Rcpp::export]]
double C_idamax2(S4 x) {
  return marlib::idamax(x);
}

// [[Rcpp::export]]
NumericVector C_dcopy(NumericVector x, NumericVector y) {
  marlib::dcopy(x, y);
  return y;
}

// [[Rcpp::export]]
NumericVector C_dcopy2(S4 x, NumericVector y) {
  marlib::dcopy(x, y);
  return y;
}

// [[Rcpp::export]]
NumericVector C_dscal(double alpha, NumericVector x) {
  marlib::dscal(alpha, x);
  return x;
}

// [[Rcpp::export]]
NumericVector C_daxpy(double alpha, NumericVector x, NumericVector y) {
  marlib::daxpy(alpha, x, y);
  return y;
}

// [[Rcpp::export]]
NumericVector C_dfill(NumericVector x, double v) {
  marlib::dfill(x, v, marlib::ArrayT());
  return x;
}

// [[Rcpp::export]]
NumericVector C_dgemv(bool trans, double alpha, NumericMatrix A, NumericVector x,
                      double beta, NumericVector y) {
  if (trans) {
    marlib::dgemv(marlib::TRANS(), alpha, A, x, beta, y, marlib::DenseMatrixT());
  } else {
    marlib::dgemv(marlib::NOTRANS(), alpha, A, x, beta, y, marlib::DenseMatrixT());
  }
  return y;
}

// [[Rcpp::export]]
NumericMatrix C_dger(bool trans, double alpha, NumericVector x, NumericVector y,
                     NumericMatrix A) {
  if (trans) {
    marlib::dger(marlib::TRANS(), alpha, x, y, A, marlib::DenseMatrixT());
  } else {
    marlib::dger(marlib::NOTRANS(), alpha, x, y, A, marlib::DenseMatrixT());
  }
  return A;
}

// [[Rcpp::export]]
NumericMatrix C_dgemm(bool transA, bool trnasB, double alpha, NumericMatrix A, NumericMatrix B,
                      double beta, NumericMatrix C) {
  if (transA == FALSE) {
    marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT());
  } else {
    marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT());
  }
  return C;
}

// [[Rcpp::export]]
List C_create_csr(NumericMatrix A) {
  int nnz = marlib::sparse_nnz(A);
  NumericVector value(nnz);
  IntegerVector rowptr(A.nrow()+1);
  IntegerVector colind(nnz);
  marlib::dense_to_sparse(A, value, rowptr, colind, 1, marlib::CSRMatrixT());
  return List::create(Named("value")=value, Named("rowptr")=rowptr, Named("colind")=colind);
}

NumericMatrix as_matrix(S4 A) {
  IntegerVector dim = A.slot("Dim");
  //  std::string classname = as<std::string>(A.attr("class"));
  NumericVector val = A.slot("x");
  val.attr("dim") = Dimension(dim[0], dim[1]);
  return as<NumericMatrix>(val);
}

marlib::csr_matrix as_csr(S4 A) {
  IntegerVector dim = A.slot("Dim");
//  std::string classname = as<std::string>(A.attr("class"));
  NumericVector val = A.slot("x");
  IntegerVector ptr = A.slot("p");
  IntegerVector ind = A.slot("j");
  return marlib::csr_matrix(dim[0], dim[1], val.length(), 0, &val[0], &ptr[0], &ind[0]);
}

// [[Rcpp::export]]
NumericMatrix C_create_dense(S4 A) {
  IntegerVector dim = A.slot("Dim");
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    NumericVector val = A.slot("x");
    val.attr("dim") = Dimension(dim[0], dim[1]);
    NumericMatrix m = as<NumericMatrix>(val);
    return m;
  } else if (classname == "dgRMatrix") {
    // NumericVector val = A.slot("x");
    // IntegerVector ptr = A.slot("p");
    // IntegerVector ind = A.slot("j");
    // marlib::csr_matrix spA(dim[0], dim[1], val.length(), 0, &val[0], &ptr[0], &ind[0]);
    NumericMatrix m(dim[0], dim[1]);
    sparse_to_dense(A, m, marlib::CSRMatrixT());
    return m;
  } else {
    return NumericMatrix(1,1);
  }
}

// [[Rcpp::export]]
NumericVector C_dgemv2(bool trans, double alpha, S4 A, NumericVector x,
                      double beta, NumericVector y) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    return NumericVector(1);
  } else if (classname == "dgRMatrix") {
    // marlib::csr_matrix m = as_csr(A);
    if (trans) {
      marlib::dgemv(marlib::TRANS(), alpha, A, x, beta, y, marlib::CSRMatrixT());
    } else {
      marlib::dgemv(marlib::NOTRANS(), alpha, A, x, beta, y, marlib::CSRMatrixT());
    }
    return y;
  } else {
    return NumericVector(1);
  }
}
