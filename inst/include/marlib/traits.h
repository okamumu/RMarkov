#ifndef MARLIB_TRAITS_H
#define MARLIB_TRAITS_H

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

  template <class T>
  struct vector_traits {
    static int size(const T& v);
    static const double* value(const T& v);
    static double* value(T& v);
    static int inc(const T& v);
  };

  template <class T>
  struct dense_matrix_traits {
    static int nrow(const T& m);
    static int ncol(const T& m);
    static const double* value(const T& m);
    static double* value(T& m);
    static int ld(const T& m);
  };

  template <class T>
  struct csr_matrix_traits {
    static int nrow(const T& m);
    static int ncol(const T& m);
    static int nnz(const T& m);
    static int base(const T& m);
    static double* value(T& m);
    static const double* value(const T& m);
    static const int* rowptr(const T& m);
    static const int* colind(const T& m);
  };

  template <class T>
  struct csc_matrix_traits {
    static int nrow(const T& m);
    static int ncol(const T& m);
    static int nnz(const T& m);
    static int base(const T& m);
    static double* value(T& m);
    static const double* value(const T& m);
    static const int* colptr(const T& m);
    static const int* rowind(const T& m);
  };

  template <class T>
  struct coo_matrix_traits {
    static int nrow(const T& m);
    static int ncol(const T& m);
    static int nnz(const T& m);
    static int base(const T& m);
    static double* value(T& m);
    static const double* value(const T& m);
    static const int* rowind(const T& m);
    static const int* colind(const T& m);
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

#endif
