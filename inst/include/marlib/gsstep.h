#ifndef GSSTEP_H
#define GSSTEP_H

namespace marlib {

/*
 !  SOR step for solving the following linear equation
 !
 !         alpha * trans(A - sigma I) * x = b
 !
 !         The step computes
 !
 !         notrans
 !         x := (D/omega + L)^(-1) (b/alpha - (U - D (1-omega)/omega - sigma I) * x)
 !
 !         trans
 !         x := (D/omega + tr(U))^(-1) (b/alpha - (tr(L) - D (1-omega)/omega - sigma I) * x)
 !
 !         A: square matrix
 !         x: vector (in; initial vector for the step, out; updated vector)
 !         b: constant vector
 */

template<typename T1, typename T2, typename T3>
void gsstep(NOTRANS, double alpha, const T1& A, double sigma,
            double omega, const T2& b, T3& x, DenseMatrixT, ArrayT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int n = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int lda = traits1::ld(A);
  const double* valueB = traits2::value(b);
  const int incb = traits2::inc(b);
  double* valueX = traits3::value(x);
  const int incx = traits3::inc(x);

  double* xptr0 = valueX;
  for (int i=1; i<=n; i++, valueA+=1, valueB+=incb, xptr0+=incx) {
    double tmpd = 0.0;
    double tmpx = *valueB / alpha;
    const double* Aix = valueA;
    double* xptr1 = valueX;
    for (int j=1; j<=n; j++, Aix+=lda, xptr1+=incx) {
      if (i == j) {
        tmpd = *Aix;
        tmpx += sigma * *xptr1;
      } else {
        tmpx -= *Aix * *xptr1;
      }
    }
    *xptr0 = (omega/tmpd) * tmpx + (1.0-omega) * *xptr0;
  }

}

template<typename T1, typename T2, typename T3>
void gsstep(TRANS, double alpha, const T1& A, double sigma,
            double omega, const T2& b, T3& x, DenseMatrixT, ArrayT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int n = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int lda = traits1::ld(A);
  const double* valueB = traits2::value(b);
  const int incb = traits2::inc(b);
  double* valueX = traits3::value(x);
  const int incx = traits3::inc(x);

  double* xptr0 = valueX;
  for (int j=1; j<=n; j++, valueA+=lda, valueB+=incb, xptr0+=incx) {
    double tmpd = 0.0;
    double tmpx = *valueB / alpha;
    const double* Ajx = valueA;
    double* xptr1 = valueX;
    for (int i=1; i<=n; i++, Ajx+=1, xptr1+=incx) {
      if (i == j) {
        tmpd = *Ajx;
        tmpx += sigma * *xptr1;
      } else {
        tmpx -= *Ajx * *xptr1;
      }
    }
    *xptr0 = (omega/tmpd) * tmpx + (1.0-omega) * *xptr0;
  }

}

template<typename T1, typename T2, typename T3>
void gsstep(NOTRANS, double alpha, const T1& A, double sigma,
            double omega, const T2& b, T3& x, CSRMatrixT, ArrayT) {
  using traits1 = csr_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int n = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* rowptr = traits1::rowptr(A);
  const int* colind = traits1::colind(A);
  const int base = traits1::base(A);
  const double* valueB = traits2::value(b);
  const int incb = traits2::inc(b);
  double* valueX = traits3::value(x);
  const int incx = traits3::inc(x);

  for (int i=0; i<n; i++, valueB+=incb) {
    double tmpd = 0.0;
    double tmpx = *valueB / alpha;
    for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
      int j = colind[z] - base;
      if (i == j) {
        tmpd = valueA[z];
        tmpx += sigma * valueX[j*incx];
      } else {
        tmpx -= valueA[z] * valueX[j*incx];
      }
    }
    valueX[i*incx] = (omega/tmpd) * tmpx + (1.0-omega) * valueX[i*incx];
  }
}

template<typename T1, typename T2, typename T3>
void gsstep(TRANS, double alpha, const T1& A, double sigma,
            double omega, const T2& b, T3& x, CSCMatrixT, ArrayT) {
  using traits1 = csc_matrix_traits<T1>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int n = traits1::nrow(A);
  // const int n = traits1::ncol(A);
  const double* valueA = traits1::value(A);
  const int* colptr = traits1::colptr(A);
  const int* rowind = traits1::rowind(A);
  const int base = traits1::base(A);
  const double* valueB = traits2::value(b);
  const int incb = traits2::inc(b);
  double* valueX = traits3::value(x);
  const int incx = traits3::inc(x);

  for (int j=0; j<n; j++, valueB+=incb) {
    double tmpd = 0.0;
    double tmpx = *valueB / alpha;
    for (int z=colptr[j]-base; z<colptr[j+1]-base; z++) {
      int i = rowind[z] - base;
      if (i == j) {
        tmpd = valueA[z];
        tmpx += sigma * valueX[j*incx];
      } else {
        tmpx -= valueA[z] * valueX[i*incx];
      }
    }
    valueX[j*incx] = (omega/tmpd) * tmpx + (1.0-omega) * valueX[j*incx];
  }
}

template<typename TR, typename T1, typename T2, typename T3, typename T4>
void gsstep(TR, double alpha, const T1& A, double sigma,
            double omega, const T2& b, T3& x, T4, DenseMatrixT) {
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = dense_matrix_traits<T3>;
  // const int m = traits2::nrow(b);
  const int n = traits2::ncol(b);
  const double* valueB = traits2::value(b);
  const int ldb = traits2::ld(b);
  double* valueX = traits3::value(x);
  const int ldx = traits3::ld(x);

  for (int i=0; i<n; i++, valueB+=ldb, valueX+=ldx) {
    gsstep(TR(), alpha, A, sigma, omega, valueB, valueX, T4(), ArrayT());
  }
}

}

#endif
