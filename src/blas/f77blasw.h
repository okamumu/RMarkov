/**
  * @file f77blasw.h
  * @brief A wrapper for f77blas
  */

#pragma once

namespace blas {

  // level 1
  void dcopy(int n, const double *x, int incx, double *y, int incy);
  void dscal(int n, double alpha, double *x, int incx);
  void daxpy(int n, double alpha, const double *x, int incx, double *y, int incy);

  double ddot(int n, const double *x, int incx, const double *y, int incy);
  double dasum(int n, const double *x, int incx);
  double dnrm2(int n, const double *x, int incx);
  int idamax(int n, const double *x, int incx);

  // level 2
  void dgemv(char trans, int m, int n, double alpha,
    const double *A, int lda, const double *x, int incx,
    double beta, double *y, int incy);

  void dger(int m, int n, double alpha,
    const double *x, int incx, const double *y, int incy,
    double *A, int lda);

  // level 3
  void dgemm(char transA, char transB,
    int m, int n, int k, double alpha,
    const double *A, int lda, const double *B, int ldb,
    double beta, double *C, int ldc);

  // use for f77blas
  #ifdef F77BLAS

  extern "C" {
    void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
    void dscal_(const int *n, const double *alpha, double *x, const int *incx);
    void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
    double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
    double dasum_(const int *n, const double *x, const int *incx);
    double dnrm2_(const int *n, const double *x, const int *incx);
    int idamax_(const int *n, const double *x, const int *incx);
    void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *lda,
      const double *x, const int *incx, const double *beta, double *y, const int *incy);
    void dger_(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
      const double *y, const int *incy, double *A, const int *lda);
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha,
      const double *A, const int *lda, const double *B, const int *ldb, const double *beta, double *C, const int *ldc);
  }

  #define __DCOPY__ dcopy_
  #define __DSCAL__ dscal_
  #define __DAXPY__ daxpy_
  #define __DDOT__ ddot_
  #define __DASUM__ dasum_
  #define __DNRM2__ dnrm2_
  #define __IDAMAX__ idamax_
  #define __DGEMV__ dgemv_
  #define __DGER__ dger_
  #define __DGEMM__ dgemm_
  #endif

  // level 1

/**
  * @bried copy from x to y
  * @param n The length of vector
  * @param x A vector
  * @param incx An increment of index of vector x
  * @param y A vector
  * @param incy An increment of index of vector y
  */

  inline
  void dcopy(int n, const double *x, int incx, double *y, int incy) {
    __DCOPY__(&n, x, &incx, y, &incy);
  }

/**
  * @brief x = alpha * x
  * @param n The length of vector x
  * @param alpha A value
  * @param x A vector
  * @param incx An increment of index of vector x
  */

  inline
  void dscal(int n, double alpha, double *x, int incx) {
    __DSCAL__(&n, &alpha, x, &incx);
  }

/**
  * @brief y = a * x + y
  * @param n The length of vector
  * @param x A vector
  * @param incx An increment of index of vector x
  * @param y A vector
  * @param incy An increment of index of vector y
  */

  inline
  void daxpy(int n, double alpha, const double *x, int incx, double *y, int incy) {
    __DAXPY__(&n, &alpha, x, &incx, y, &incy);
  }

/**
  * @brief dot product
  * @param n The length of vector
  * @param x A vector
  * @param incx An increment of index of vector x
  * @param y A vector
  * @param incy An increment of index of vector y
  * @return A value of the dot product of x and y
  */

  inline
  double ddot(int n, const double *x, int incx, const double *y, int incy) {
    return __DDOT__(&n, x, &incx, y, &incy);
  }

/**
  * @brief A value of absolute sum of vector as l1norm
  * @param n The length of vector
  * @param x A vector
  * @param incx An increment of index of vector x
  * @return A value of absolute sum
  */

  inline
  double dasum(int n, const double *x, int incx) {
    return __DASUM__(&n, x, &incx);
  }

/**
  * @brief A value of l2norm of x
  * @param n The length of vector
  * @param x A vector
  * @param incx An increment of index of vector x
  * @return A value of l2norm
  */

  inline
  double dnrm2(int n, const double *x, int incx) {
    return __DNRM2__(&n, x, &incx);
  }

/**
  * @brief Get a pointer of the absolute maximum of x as L-inf-norm
  * @param n The length of vector
  * @param x A vector
  * @param incx An increment of index of venctor x
  * @return The pointer of the absolute maximum
  */

  inline
  int idamax(int n, const double *x, int incx) {
    return __IDAMAX__(&n, x, &incx) - 1;
  }

  // level 2

/**
  * @brief y = alpha * trans(A) * x + beta * y
  * @param trans a charactor to indicate the transpose of matrix A. If trans is 'T', the matrix A is transposed. Otherwise, if trans is 'N', the matrix A is not indicated.
  * @param m the number of rows of matrix A
  * @param n the number of columns of matrix A
  * @param alpha a value
  * @param A a matrix
  * @param lda length of column data of matrix A
  * @param x a vector
  * @param incx the increment of index of vector x
  * @param beta a value
  * @param y a vector
  * @param incy increment of index of vector y
  */

  inline
  void dgemv(char trans, int m, int n, double alpha,
    const double *A, int lda, const double *x, int incx,
    double beta, double *y, int incy) {
    __DGEMV__(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

/**
  * @brief A = alpha * x * y + A
  * @param m the number of rows of matrix A
  * @param n the number of columns of matrix A
  * @param alpha a value
  * @param x a vector
  * @param incx the increment of index of vector x
  * @param y a vector
  * @param incy increment of index of vector y
  * @param A a matrix
  * @param lda length of column data of matrix A
  */

  inline
  void dger(int m, int n, double alpha,
    const double *x, int incx, const double *y, int incy,
    double *A, int lda) {
    __DGER__(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  // level 3

/**
  * @brief C = alpha * transa(A) * transb(B) + beta * C
  * @param transa a charactor to indicate the transpose of matrix A. If trans is 'T', the matrix A is transposed. Otherwise, if trans is 'N', the matrix A is not indicated.
  * @param transb a charactor to indicate the transpose of matrix B. If trans is 'T', the matrix B is transposed. Otherwise, if trans is 'N', the matrix B is not indicated.
  * @param m the number of rows of matrix C
  * @param n the number of columns of matrix C
  * @param alpha a value
  * @param A a matrix
  * @param lda length of column data of matrix A
  * @param B a matrix
  * @param ldb length of column data of matrix B
  * @param beta a value
  * @param C a matrix
  * @param ldc length of column data of matrix C
  */

  inline
  void dgemm(char transA, char transB,
    int m, int n, int k, double alpha,
    const double *A, int lda, const double *B, int ldb,
    double beta, double *C, int ldc) {
    __DGEMM__(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  }

}
