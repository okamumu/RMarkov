#ifndef BLAS_CSC_H
#define BLAS_CSC_H

namespace marlib {

template <typename T1, typename T2, typename T3, typename T4>
int dense_to_sparse(const T1& A, T2& spA, T3& colptr, T4& rowind, int base, CSCMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const int lda = traits1::ld(A);

  int z = 0, zz = base;
  colptr[0] = zz;
  for (int j=0, p=0; j<n; j++, p+=lda) {
    for (int i=0, idx=base; i<m; i++, idx++) {
      if (std::fpclassify(A[i+p]) != FP_ZERO) {
        spA[z] = A[i+p];
        rowind[z] = idx;
        z++;
        zz++;
      }
    }
    colptr[j+1] = zz;
  }
  return z;
}

template <typename T1, typename T2>
void sparse_to_dense(const T1& X, T2& Y, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  const int base = traits1::base(X);
  const int m = traits1::nrow(X);
  const int n = traits1::ncol(X);
  const double* valueX = traits1::value(X);
  const int* colptr = traits1::colptr(X);
  const int* rowind = traits1::rowind(X);
  double* valueY = traits2::value(Y);
  const int ldy = traits2::ld(Y);

  int z = 0, zz = base;
  for (int j=0, p=0, jj=base; j<n; j++, jj++, p+=ldy) {
    for (int i=0, ii=base; i<m; i++, ii++) {
      if (colptr[j] <= zz && zz < colptr[j+1] && ii == rowind[z]) {
        valueY[i+p] = valueX[z];
        z++;
        zz++;
      } else {
        valueY[i+p] = 0.0;
      }
    }
  }
}

// level 2

template <typename T1, typename T2, typename T3>
void dgemv(NOTRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int base = traits1::base(A);
  // const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);
  double* valueY = traits3::value(y);
  const int incy = traits3::inc(y);

  dscal(beta, y);
  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      valueY[i*incy] += alpha * valueA[z] * valueX[j*incx];
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemv(TRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int base = traits1::base(A);
  // const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);
  double* valueY = traits3::value(y);
  const int incy = traits3::inc(y);

  dscal(beta, y);
  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      valueY[j*incy] += alpha * valueA[z] * valueX[i*incx];
    }
  }
}

template <typename T1, typename T2, typename T3>
void dger(NOTRANS, double alpha, const T1& x, const T2& y, T3& A, CSCMatrixT) {
  using traits1 = csc_matrix_traits<T3>;
  using traits2 = vector_traits<T1>;
  using traits3 = vector_traits<T2>;
  const int base = traits1::base(A);
  // const int m = traits1::nrow(A);
  const int n = traits1::ncol(A);
  double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueX = traits2::value(x);
  const int incx = traits2::inc(x);
  const double* valueY = traits3::value(y);
  const int incy = traits3::inc(y);

  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      valueA[z] += alpha * valueX[i*incx] * valueY[j*incy];
    }
  }
}

template <typename T1, typename T2, typename T3>
void dger(TRANS, double alpha, const T1& x, const T2& y, T3& A, CSCMatrixT) {
  dger(NOTRANS(), alpha, y, x, A, CSCMatrixT());
}

// level 3

template <typename T1, typename T2, typename T3>
void dgemm(NOTRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSCMatrixT, DenseMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int base = traits1::base(A);
  const int k = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueB = traits2::value(B);
  const int ldb = traits2::ld(B);
  // const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<k; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Bptr = valueB;
      double* Cptr = valueC;
      for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
        Cptr[i] += alpha * valueA[z] * Bptr[j];
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemm(TRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSCMatrixT, DenseMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int base = traits1::base(A);
  const int m = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueB = traits2::value(B);
  const int ldb = traits2::ld(B);
  // const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<m; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Bptr = valueB;
      double* Cptr = valueC;
      for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
        Cptr[j] += alpha * valueA[z] * Bptr[i];
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemm(NOTRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSCMatrixT, DenseMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int base = traits1::base(A);
  const int k = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueB = traits2::value(B);
  const int ldb = traits2::ld(B);
  // const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<k; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Bptr = valueB;
      double* Cptr = valueC;
      for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
        Cptr[i] += alpha * valueA[z] * Bptr[j*ldb];
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemm(TRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSCMatrixT, DenseMatrixT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int base = traits1::base(A);
  const int m = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const double* valueB = traits2::value(B);
  const int ldb = traits2::ld(B);
  // const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<m; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Bptr = valueB;
      double* Cptr = valueC;
      for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
        Cptr[j] += alpha * valueA[z] * Bptr[i*ldb];
      }
    }
  }
}

// dcsrmm2

template <typename T1, typename T2, typename T3>
void dgemm(NOTRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSCMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = csc_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  // const int m = traits1::nrow(A);
  // const int k = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int lda = traits1::ld(A);
  // const int k = traits2::nrow(B);
  // const int n = traits2::ncol(B);
  const double* valueB = traits2::value(B);
  const int* colptr = traits2::colptr(B);
  const int* rowind = traits2::rowind(B);
  const int base = traits2::base(B);
  const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Aptr = valueA;
      double* Cptr = valueC;
      for (int v=0; v<m; v++, Aptr+=1, Cptr+=1) {
        Cptr[j*ldc] += alpha * Aptr[i*lda] * valueB[z];
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemm(TRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSCMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = csc_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  // const int k = traits1::nrow(A);
  // const int m = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int lda = traits1::ld(A);
  // const int k = traits2::nrow(B);
  // const int n = traits2::ncol(B);
  const double* valueB = traits2::value(B);
  const int* colptr = traits2::colptr(B);
  const int* rowind = traits2::rowind(B);
  const int base = traits2::base(B);
  const int m = traits3::nrow(C);
  const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Aptr = valueA;
      double* Cptr = valueC;
      for (int v=0; v<m; v++, Aptr+=lda, Cptr+=1) {
        Cptr[j*ldc] += alpha * Aptr[i] * valueB[z];
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemm(NOTRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSCMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = csc_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  // const int m = traits1::nrow(A);
  const int k = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int lda = traits1::ld(A);
  // const int n = traits2::nrow(B);
  // const int k = traits2::ncol(B);
  const double* valueB = traits2::value(B);
  const int* colptr = traits2::colptr(B);
  const int* rowind = traits2::rowind(B);
  const int base = traits2::base(B);
  const int m = traits3::nrow(C);
  // const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<k; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Aptr = valueA;
      double* Cptr = valueC;
      for (int v=0; v<m; v++, Aptr+=1, Cptr+=1) {
        Cptr[i*ldc] += alpha * Aptr[j*lda] * valueB[z];
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void dgemm(TRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSCMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = csc_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  const int k = traits1::nrow(A);
  // const int m = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int lda = traits1::ld(A);
  // const int n = traits2::nrow(B);
  // const int k = traits2::ncol(B);
  const double* valueB = traits2::value(B);
  const int* colptr = traits2::colptr(B);
  const int* rowind = traits2::rowind(B);
  const int base = traits2::base(B);
  const int m = traits3::nrow(C);
  // const int n = traits3::ncol(C);
  double* valueC = traits3::value(C);
  const int ldc = traits3::ld(C);

  dscal(beta, C);
  for (int j=0; j<k; j++) {
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      const double* Aptr = valueA;
      double* Cptr = valueC;
      for (int v=0; v<m; v++, Aptr+=lda, Cptr+=1) {
        Cptr[i*ldc] += alpha * Aptr[j] * valueB[z];
      }
    }
  }
}

}

#endif

