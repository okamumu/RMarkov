#ifndef MARLIB_BLAS2_H
#define MARLIB_BLAS2_H

namespace marlib {

template<typename T1, typename T2, typename T3>
void dgemv(TRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const int ld = traits1::ld(A);
  const int inc1 = traits2::inc(x);
  const int inc2 = traits3::inc(y);
  f77dgemv('T', m, n, alpha, traits1::value(A), ld,
              traits2::value(x), inc1, beta, traits3::value(y), inc2);
}

template<typename T1, typename T2, typename T3>
void dgemv(NOTRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const int ld = traits1::ld(A);
  const int inc1 = traits2::inc(x);
  const int inc2 = traits3::inc(y);
  f77dgemv('N', m, n, alpha, traits1::value(A), ld,
              traits2::value(x), inc1, beta, traits3::value(y), inc2);
}

template<typename T1, typename T2, typename T3>
void dger(NOTRANS, double alpha, const T1& x, const T2& y, T3& A, DenseMatrixT) {
  using traits1 = vector_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int inc1 = traits1::inc(x);
  const int inc2 = traits2::inc(y);
  const int m = traits3::nrow(A);
  const int n = traits3::ncol(A);
  const int ld = traits3::ld(A);
  f77dger(m, n, alpha, traits1::value(x), inc1, traits2::value(y), inc2,
             traits3::value(A), ld);
}

template<typename T1, typename T2, typename T3>
void dger(TRANS, double alpha, const T1& x, const T2& y, T3& A, DenseMatrixT) {
  using traits1 = vector_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int inc1 = traits1::inc(x);
  const int inc2 = traits2::inc(y);
  const int m = traits3::nrow(A);
  const int n = traits3::ncol(A);
  const int ld = traits3::ld(A);
  f77dger(m, n, alpha, traits2::value(y), inc2, traits1::value(x), inc1,
             traits3::value(A), ld);
}

}

#endif
