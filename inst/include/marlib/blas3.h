#ifndef BLAS3_H
#define BLAS3_H

namespace marlib {

template<typename T1, typename T2, typename T3>
void dgemm(NOTRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int k = traits1::ncol(A);
  const int ld1 = traits1::ld(A);
  const int ld2 = traits2::ld(B);
  const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  const int ld3 = traits3::ld(C);
  f77dgemm('N', 'N', m, n, k, alpha, traits1::value(A), ld1,
              traits2::value(B), ld2, beta, traits3::value(C), ld3);
}

template<typename T1, typename T2, typename T3>
void dgemm(TRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int k = traits1::nrow(A);
  const int ld1 = traits1::ld(A);
  const int ld2 = traits2::ld(B);
  const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  const int ld3 = traits3::ld(C);
  f77dgemm('T', 'N', m, n, k, alpha, traits1::value(A), ld1,
              traits2::value(B), ld2, beta, traits3::value(C), ld3);
}

template<typename T1, typename T2, typename T3>
void dgemm(NOTRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int k = traits1::ncol(A);
  const int ld1 = traits1::ld(A);
  const int ld2 = traits2::ld(B);
  const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  const int ld3 = traits3::ld(C);
  f77dgemm('N', 'T', m, n, k, alpha, traits1::value(A), ld1,
              traits2::value(B), ld2, beta, traits3::value(C), ld3);
}

template<typename T1, typename T2, typename T3>
void dgemm(TRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int k = traits1::nrow(A);
  const int ld1 = traits1::ld(A);
  const int ld2 = traits2::ld(B);
  const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  const int ld3 = traits3::ld(C);
  f77dgemm('T', 'T', m, n, k, alpha, traits1::value(A), ld1,
              traits2::value(B), ld2, beta, traits3::value(C), ld3);
}

}

#endif

