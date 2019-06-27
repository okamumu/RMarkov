#ifndef MARLIB_RCPP_H
#define MARLIB_RCPP_H

#include <Rcpp.h>
#include <marlib.h>

namespace marlib {

template <>
struct vector_traits<Rcpp::NumericVector> {
  static int size(const Rcpp::NumericVector& v) { return v.size(); }
  static const double* value(const Rcpp::NumericVector& v) { return &v[0]; }
  static double* value(Rcpp::NumericVector& v) { return &v[0]; }
  static int inc(const Rcpp::NumericVector& v) { return 1; }
};

template <>
struct vector_traits<Rcpp::NumericMatrix> {
  static int size(const Rcpp::NumericMatrix& v) { return v.size(); }
  static const double* value(const Rcpp::NumericMatrix& v) { return &v[0]; }
  static double* value(Rcpp::NumericMatrix& v) { return &v[0]; }
  static int inc(const Rcpp::NumericMatrix& v) { return 1; }
};

template <>
struct dense_matrix_traits<Rcpp::NumericMatrix> {
  static int nrow(const Rcpp::NumericMatrix& m) { return m.nrow(); }
  static int ncol(const Rcpp::NumericMatrix& m) { return m.ncol(); }
  static const double* value(const Rcpp::NumericMatrix& m) { return &m[0]; }
  static double* value(Rcpp::NumericMatrix& m) { return &m[0]; }
  static int ld(const Rcpp::NumericMatrix& m) { return nrow(m); }
};

template <>
struct vector_traits<Rcpp::S4> {
  static int size(const Rcpp::S4& v) { return Rcpp::as<Rcpp::NumericVector>(v.slot("x")).length(); }
  static const double* value(const Rcpp::S4& v) { return &Rcpp::as<Rcpp::NumericVector>(v.slot("x"))[0]; }
  static double* value(Rcpp::S4& v) { return &Rcpp::as<Rcpp::NumericVector>(v.slot("x"))[0]; }
  static int inc(const Rcpp::S4& v) { return 1; }
};

template <>
struct dense_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static int ld(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
};

template <>
struct csr_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static int nnz(const Rcpp::S4& m) { return Rcpp::as<Rcpp::NumericVector>(m.slot("x")).length(); }
  static int base(const Rcpp::S4& m) { return 0; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const int* rowptr(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("p"))[0]; }
  static const int* colind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("j"))[0]; }
};

template <>
struct csc_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static int nnz(const Rcpp::S4& m) { return Rcpp::as<Rcpp::NumericVector>(m.slot("x")).length(); }
  static int base(const Rcpp::S4& m) { return 0; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const int* colptr(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("p"))[0]; }
  static const int* rowind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("i"))[0]; }
};

template <>
struct coo_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static int nnz(const Rcpp::S4& m) { return Rcpp::as<Rcpp::NumericVector>(m.slot("x")).length(); }
  static int base(const Rcpp::S4& m) { return 0; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const int* rowind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("i"))[0]; }
  static const int* colind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("j"))[0]; }
};

}

#endif
