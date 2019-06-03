#pragma once

#include "traits.h"

namespace marlib {

  template <typename T1, typename T2, typename T3, typename T4>
  int dense_to_sparse(const T1& A, T2& spA, T3& rowptr, T4& colind, int base, CSRMatrixT) {
    using traits1 = dense_matrix_traits<T1>;
    const int m = traits1::nrow(A);
    const int n = traits1::ncol(A);
    const double* valueA = traits1::value(A);
    const int lda = traits1::ld(A);
    int z = 0, zz = base;
    rowptr[0] = zz;
    for (int i=0; i<m; i++) {
      for (int j=0, idx=base, p=0; j<n; j++, idx++, p+=lda) {
        if (std::fpclassify(valueA[i+p]) != FP_ZERO) {
          spA[z] = valueA[i+p];
          colind[z] = idx;
          z++;
          zz++;
        }
      }
      rowptr[i+1] = zz;
    }
    return z;
  }

  template <typename T1, typename T2>
  void sparse_to_dense(const T1& X, T2& Y, CSRMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = dense_matrix_traits<T2>;
    const int base = traits1::base(X);
    const int m = traits1::nrow(X);
    const int n = traits1::ncol(X);
    const double* valueX = traits1::value(X);
    const int* rowptr = traits1::rowptr(X);
    const int* colind = traits1::colind(X);
    double* valueY = traits2::value(Y);
    const int ld = traits2::ld(Y);
    int z = 0, zz = base;
    for (int i=0, ii=base; i<m; i++, ii++) {
      for (int j=0, p=0, jj=base; j<n; j++, jj++, p+=ld) {
        if (rowptr[i] <= zz && zz < rowptr[i+1] && jj == colind[z]) {
          valueY[i+p] = valueX[z];
          z++;
          zz++;
        } else {
          valueY[i+p] = 0.0;
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemv(NOTRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, CSRMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = vector_traits<T2>;
    using traits3 = vector_traits<T3>;
    const int base = traits1::base(A);
    const int m = traits1::nrow(A);
    // const int n = traits1::ncol(A);
    const double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const double* valueX = traits2::value(x);
    const int incx = traits2::inc(x);
    double* valueY = traits3::value(y);
    const int incy = traits3::inc(y);

    marlib::dscal(beta, y);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        valueY[i*incy] += alpha * valueA[z] * valueX[j*incx];
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemv(TRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, CSRMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = vector_traits<T2>;
    using traits3 = vector_traits<T3>;
    const int base = traits1::base(A);
    const int m = traits1::nrow(A);
    // const int n = traits1::ncol(A);
    const double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const double* valueX = traits2::value(x);
    const int incx = traits2::inc(x);
    double* valueY = traits3::value(y);
    const int incy = traits3::inc(y);

    marlib::dscal(beta, y);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        valueY[j*incy] += alpha * valueA[z] * valueX[i*incx];
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dger(NOTRANS, double alpha, const T1& x, const T2& y, T3& A, CSRMatrixT) {
    using traits1 = csr_matrix_traits<T3>;
    using traits2 = vector_traits<T1>;
    using traits3 = vector_traits<T2>;
    const int base = traits1::base(A);
    const int m = traits1::nrow(A);
    // const int n = traits1::ncol(A);
    double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const double* valueX = traits2::value(x);
    const int incx = traits2::inc(x);
    const double* valueY = traits3::value(y);
    const int incy = traits3::inc(y);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        valueA[z] += alpha * valueX[i*incx] * valueY[j*incy];
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dger(TRANS, double alpha, const T1& x, const T2& y, T3& A, CSRMatrixT) {
    using traits1 = csr_matrix_traits<T3>;
    using traits2 = vector_traits<T1>;
    using traits3 = vector_traits<T2>;
    const int base = traits1::base(A);
    const int m = traits1::nrow(A);
    // const int n = traits1::ncol(A);
    double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const double* valueX = traits2::value(x);
    const int incx = traits2::inc(x);
    const double* valueY = traits3::value(y);
    const int incy = traits3::inc(y);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        valueA[z] += alpha * valueY[i*incx] * valueX[j*incy];
      }
    }
  }
/*
  inline void dcsrmmNN(int m, int n, int, double alpha,
    const double *A, const int *rowptr, const int *colind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j];
        }
      }
    }
  }

  inline void dcsrmmTN(int m, int n, int k, double alpha,
    const double *A, const int *rowptr, const int *colind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i];
        }
      }
    }
  }

  inline void dcsrmmNT(int m, int n, int, double alpha,
    const double *A, const int *rowptr, const int *colind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j*ldb];
        }
      }
    }
  }

  inline void dcsrmmTT(int m, int n, int k, double alpha,
    const double *A, const int *rowptr, const int *colind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i*ldb];
        }
      }
    }
  }
*/

}
