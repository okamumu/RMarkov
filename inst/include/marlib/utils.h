#pragma once

#include "traits.h"

namespace marlib {

template <typename T1>
int sparse_nnz(const T1& A) {
  using traits1 = dense_matrix_traits<T1>;
  const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const double* value = traits1::value(A);
  const int lda = traits1::ld(A);
  int nnz = 0;
  for (int j=0, cn=0; j<n; j++, cn+=lda) {
    for (int i=0; i<m; i++) {
      if (std::fpclassify(value[i+cn]) != FP_ZERO) {
        nnz++;
      }
    }
  }
  return nnz;
}

/// matrix copy

template<typename T1, typename T2, typename MatT>
void dcopy(const T1& x, T2& y, MatT, DenseMatrixT) {
  sparse_to_dense(x, y, MatT());
}

template<typename T1, typename T2>
void dcopy(const T1& x, T2& y, DenseMatrixT, DenseMatrixT) {
  dcopy(x,y);
}

/// matrix dim

template<typename T1>
int nrow(const T1& X, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  return traits1::nrow(X);
}

template<typename T1>
int nrow(const T1& X, CSRMatrixT) {
  using traits1 = csr_matrix_traits<T1>;
  return traits1::nrow(X);
}

template<typename T1>
int nrow(const T1& X, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  return traits1::nrow(X);
}

template<typename T1>
int nrow(const T1& X, COOMatrixT) {
  using traits1 = coo_matrix_traits<T1>;
  return traits1::nrow(X);
}

template<typename T1>
int ncol(const T1& X, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  return traits1::ncol(X);
}

template<typename T1>
int ncol(const T1& X, CSRMatrixT) {
  using traits1 = csr_matrix_traits<T1>;
  return traits1::ncol(X);
}

template<typename T1>
int ncol(const T1& X, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  return traits1::ncol(X);
}

template<typename T1>
int ncol(const T1& X, COOMatrixT) {
  using traits1 = coo_matrix_traits<T1>;
  return traits1::ncol(X);
}

/// diag

template <typename T1, typename T2>
void diag_get(const T1& A, T2& x, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  const int n = std::min(traits1::nrow(A), traits1::ncol(A));
  const double* value = traits1::value(A);
  const int lda = traits1::ld(A);
  using traits2 = vector_traits<T2>;
  double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int i=0; i<n; i++, value+=lda+1, valueX+=incx) {
    *valueX = *value;
  }
}

template <typename T1, typename T2>
void diag_get(const T1& A, T2& x, CSRMatrixT) {
  using traits1 = csr_matrix_traits<T1>;
  const int m = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const double* value = traits1::value(A);
  const int* rowptr = traits1::rowptr(A);
  const int* colind = traits1::colind(A);
  const int base = traits1::base(A);
  using traits2 = vector_traits<T2>;
  double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int i=0; i<m; i++, valueX+=incx) {
    *valueX = 0.0;
    for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
      int j = colind[z] - base;
      if (i > j) {
        continue;
      } else if (i == j) {
        *valueX = value[z];
        break;
      } else if (i < j) {
        break;
      }
    }
  }
}

template <typename T1, typename T2>
void diag_get(const T1& A, T2& x, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  // const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const double* value = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const int base = traits1::base(A);
  using traits2 = vector_traits<T2>;
  double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int j=0; j<n; j++, valueX+=incx) {
    *valueX = 0.0;
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      if (i < j) {
        continue;
      } else if (i == j) {
        *valueX = value[z];
        break;
      } else if (i > j) {
        break;
      }
    }
  }
}

template <typename T1, typename T2>
void diag_get(const T1& A, T2& x, COOMatrixT) {
  using traits1 = coo_matrix_traits<T1>;
  // const int m = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const int nnz = traits1::nnz(A);
  const double* value = traits1::value(A);
  const int* rowind = traits1::rowind(A);
  const int* colind = traits1::colind(A);
  const int base = traits1::base(A);
  using traits2 = vector_traits<T2>;
  double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  dfill(x, 0.0);
  for (int z=0; z<nnz; z++) {
    if (rowind[z] == colind[z]) {
      valueX[(rowind[z]-base)*incx] = value[z];
    }
  }
}

/// diag2

template <typename T1, typename T2>
void diag_set(T1& A, const T2& x, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  const int n = std::min(traits1::nrow(A), traits1::ncol(A));
  double* value = traits1::value(A);
  const int lda = traits1::ld(A);
  using traits2 = vector_traits<T2>;
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int i=0; i<n; i++, value+=lda+1, valueX+=incx) {
    *value = *valueX;
  }
}

template <typename T1, typename T2>
void diag_set(T1& A, const T2& x, CSRMatrixT) {
  using traits1 = csr_matrix_traits<T1>;
  const int m = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  double* value = traits1::value(A);
  const int* rowptr = traits1::rowptr(A);
  const int* colind = traits1::colind(A);
  const int base = traits1::base(A);
  using traits2 = vector_traits<T2>;
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int i=0; i<m; i++, valueX+=incx) {
    for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
      int j = colind[z] - base;
      if (i > j) {
        continue;
      } else if (i == j) {
        value[z] = *valueX;
        break;
      } else if (i < j) {
        break;
      }
    }
  }
}

template <typename T1, typename T2>
void diag_set(T1& A, const T2& x, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  // const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  double* value = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const int base = traits1::base(A);
  using traits2 = vector_traits<T2>;
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int j=0; j<n; j++, valueX+=incx) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      if (i < j) {
        continue;
      } else if (i == j) {
        value[z] = *valueX;
        break;
      } else if (i > j) {
        break;
      }
    }
  }
}

template <typename T1, typename T2>
void diag_set(T1& A, const T2& x, COOMatrixT) {
  using traits1 = coo_matrix_traits<T1>;
  // const int m = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const int nnz = traits1::nnz(A);
  double* value = traits1::value(A);
  const int* rowind = traits1::rowind(A);
  const int* colind = traits1::colind(A);
  const int base = traits1::base(A);
  using traits2 = vector_traits<T2>;
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);

  for (int z=0; z<nnz; z++) {
    if (rowind[z] == colind[z]) {
      value[z] = valueX[(rowind[z]-base)*incx];
    }
  }
}

}


