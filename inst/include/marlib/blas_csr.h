#ifndef MARLIB_BLAS_CSR_H
#define MARLIB_BLAS_CSR_H

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

  // level 2

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

    dscal(beta, y);
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

    dscal(beta, y);
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
    dger(NOTRANS(), alpha, y, x, A, CSRMatrixT());
  }

  // level 3

  template <typename T1, typename T2, typename T3>
  void dgemm(NOTRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSRMatrixT, DenseMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = dense_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const int base = traits1::base(A);
    const double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const double* valueB = traits2::value(B);
    const int ldb = traits2::ld(B);
    const int m = traits3::nrow(C);
    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        const double* Bptr = valueB;
        double* Cptr = valueC;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[i] += alpha * valueA[z] * Bptr[j];
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemm(TRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSRMatrixT, DenseMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = dense_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const int base = traits1::base(A);
    const double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const int k = traits2::nrow(B);
    const double* valueB = traits2::value(B);
    const int ldb = traits2::ld(B);
    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        const double* Bptr = valueB;
        double* Cptr = valueC;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[j] += alpha * valueA[z] * Bptr[i];
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemm(NOTRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSRMatrixT, DenseMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = dense_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const int base = traits1::base(A);
    const double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const double* valueB = traits2::value(B);
    const int ldb = traits2::ld(B);
    const int m = traits3::nrow(C);
    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        const double* Bptr = valueB;
        double* Cptr = valueC;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[i] += alpha * valueA[z] * Bptr[j*ldb];
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemm(TRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, CSRMatrixT, DenseMatrixT) {
    using traits1 = csr_matrix_traits<T1>;
    using traits2 = dense_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const int base = traits1::base(A);
    const double* valueA = traits1::value(A);
    const int* rowptr = traits1::rowptr(A);
    const int* colind = traits1::colind(A);
    const int k = traits2::ncol(B);
    const double* valueB = traits2::value(B);
    const int ldb = traits2::ld(B);
//    const int m = traits3::nrow(C);
    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
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
  void dgemm(NOTRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSRMatrixT) {
    using traits1 = dense_matrix_traits<T1>;
    using traits2 = csr_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const double* valueA = traits1::value(A);
    const int lda = traits1::ld(A);
    const int k = traits2::nrow(B);
    const int base = traits2::base(B);
    const double* valueB = traits2::value(B);
    const int* rowptr = traits2::rowptr(B);
    const int* colind = traits2::colind(B);
    const int m = traits3::nrow(C);
//    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        const double* Aptr = valueA;
        double* Cptr = valueC;
        for (int v=0; v<m; v++, Aptr+=1, Cptr+=1) {
          Cptr[j*ldc] += alpha * Aptr[i*lda] * valueB[z];
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemm(TRANS, NOTRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSRMatrixT) {
    using traits1 = dense_matrix_traits<T1>;
    using traits2 = csr_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const double* valueA = traits1::value(A);
    const int lda = traits1::ld(A);
    const int k = traits2::nrow(B);
    const int base = traits2::base(B);
    const double* valueB = traits2::value(B);
    const int* rowptr = traits2::rowptr(B);
    const int* colind = traits2::colind(B);
    const int m = traits3::nrow(C);
    //    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        const double* Aptr = valueA;
        double* Cptr = valueC;
        for (int v=0; v<m; v++, Aptr+=lda, Cptr+=1) {
          Cptr[j*ldc] += alpha * Aptr[i] * valueB[z];
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemm(NOTRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSRMatrixT) {
    using traits1 = dense_matrix_traits<T1>;
    using traits2 = csr_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const double* valueA = traits1::value(A);
    const int lda = traits1::ld(A);
    const int n = traits2::nrow(B);
    const int base = traits2::base(B);
    const double* valueB = traits2::value(B);
    const int* rowptr = traits2::rowptr(B);
    const int* colind = traits2::colind(B);
    const int m = traits3::nrow(C);
    //    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<n; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
        const double* Aptr = valueA;
        double* Cptr = valueC;
        for (int v=0; v<m; v++, Aptr+=1, Cptr+=1) {
          Cptr[i*ldc] += alpha * Aptr[j*lda] * valueB[z];
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void dgemm(TRANS, TRANS, double alpha, const T1& A, const T2& B, double beta, T3& C, DenseMatrixT, CSRMatrixT) {
    using traits1 = dense_matrix_traits<T1>;
    using traits2 = csr_matrix_traits<T2>;
    using traits3 = dense_matrix_traits<T3>;
    const double* valueA = traits1::value(A);
    const int lda = traits1::ld(A);
    const int n = traits2::nrow(B);
    const int base = traits2::base(B);
    const double* valueB = traits2::value(B);
    const int* rowptr = traits2::rowptr(B);
    const int* colind = traits2::colind(B);
    const int m = traits3::nrow(C);
    //    const int n = traits3::ncol(C);
    double* valueC = traits3::value(C);
    const int ldc = traits3::ld(C);

    dscal(beta, C);
    for (int i=0; i<n; i++) {
      for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
        int j = colind[z] - base;
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
