#pragma once

namespace marlib {

  // tags
  struct TRANS{};
  struct NOTRANS{};

  struct ArrayT{};
  struct DenseMatrixT{};
  struct CSRMatrixT{};
  struct CSCMatrixT{};
  struct COOMatrixT{};

  struct KRONMatrixT{};

  //   template <class T>
  //   struct get_category;
  //
  //   template <>
  //   struct get_category<> {
  //   using type = constant_value_tag;
  // };
  //
  // /**
  // * @brief A meta function to get a type as the tag.
  // * The return value is given by 'type' type.
  // * This is an implementation for double.
  // */
  // template <>
  // struct get_category<double> {
  //   using type = constant_value_tag;
  // };
  //

  // traits: std::vector & NumericVector in Rcpp
  template <class T>
  struct vector_traits {
    static int size(const T& v) { return v.size(); }
    static const double* value(const T& v) { return &v[0]; }
    static double* value(T& v) { return &v[0]; }
    static int inc(const T& v) { return 1; }
  };

  // traits: NumericMatrix in Rcpp
  template <class T>
  struct dense_matrix_traits {
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static const double* value(const T& m) { return &m[0]; }
    static double* value(T& m) { return &m[0]; }
    static int ld(const T& m) { return m.nrow(); }
  };

  template <class T>
  struct csr_matrix_traits {
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static int nnz(const T& m) { return m.nnz(); }
    static int base(const T& m) { return m.base(); }
    static double* value(T& m) { return &m[0]; }
    static const double* value(const T& m) { return &m[0]; }
    static const int* rowptr(const T& m) { return &m.rowptr(0); }
    static const int* colind(const T& m) { return &m.colind(0); }
  };

  // traits

  template <>
  struct vector_traits<array> {
    static int size(const array& v) { return v.size(); }
    static const double* value(const array& v) { return &v[0]; }
    static double* value(array& v) { return &v[0]; }
    static int inc(const array& v) { return v.inc(); }
  };

  // specialize for double*
  template <>
  struct vector_traits<double*> {
    static const double* value(const double* v) { return &v[0]; }
    static double* value(double* v) { return &v[0]; }
    static int inc(const double* v) { return 1; }
  };

  template <>
  struct vector_traits<const double*> {
    static const double* value(const double* v) { return &v[0]; }
    static int inc(const double* v) { return 1; }
  };

  template <class T>
  struct csc_matrix_traits {
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static int nnz(const T& m) { return m.nnz(); }
    static int base(const T& m) { return m.base(); }
    static double* value(T& m) { return &m[0]; }
    static const double* value(const T& m) { return &m[0]; }
    static const int* colptr(const T& m) { return &m.colptr(0); }
    static const int* rowind(const T& m) { return &m.rowind(0); }
  };

  template <class T>
  struct coo_matrix_traits {
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static int nnz(const T& m) { return m.nnz(); }
    static int base(const T& m) { return m.base(); }
    static double* value(T& m) { return &m[0]; }
    static const double* value(const T& m) { return &m[0]; }
    static const int* rowind(const T& m) { return &m.rowind(0); }
    static const int* colind(const T& m) { return &m.colind(0); }
  };

  // specialize for S4 class in Matrix package
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

  // kron

  // template <class T>
  // struct kron_matrix_traits {
  //   using MatrixT = DenseMatrixT;
  //   static int nrow(const T& m) { return m.nrow(); }
  //   static int ncol(const T& m) { return m.ncol(); }
  //   static int left(const T& m) { return m.left(); }
  //   static int right(const T& m) { return m.right(); }
  //   static double* value(T& m) { return &m[0]; }
  //   static const double* value(const T& m) { return &m[0]; }
  //   static const const_dense_matrix& gen(const T& m) {
  //
  //   }
  // };

    // /**
  //  * @brief A meta function of if-then-else. If Cond is true, the return
  //  * value (type) becomes a type given by Then.
  //  */
  // template<bool Cond, typename Then, typename Else>
  // struct if_ {
  //   using type = Then;
  // };
  //
  // /**
  // * @brief A meta function of if-then-else. If Cond is false, the return
  // * value (type) becomes a type given by Else.
  // */
  // template<typename Then, typename Else>
  // struct if_<false, Then, Else> {
  //   using type = Else;
  // };
  //
}
