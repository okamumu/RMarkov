#ifndef MARLIB_CTMC_STSEN_H
#define MARLIB_CTMC_STSEN_H

namespace marlib {

template <typename T1, typename T2, typename T3, typename T4,
          typename Params, typename Func, typename MatT>
void ctmc_stsen_gs(const T1& Q, T2& x, const T3& b, const T4& pis,
                   Params& params, const Func& callback, MatT) {
  using traits1 = vector_traits<T2>;
  int n = traits1::size(x);
  std::vector<double> prevx(n);
  params.iter = 0;
  params.info = 1;
  while(1) {
    dcopy(x, prevx);
    for (int i=0; i<params.steps; i++) {
      gsstep(TRANS(), -1.0, Q, 0.0, 1.0, b, x, MatT(), ArrayT());
      double tmp = dsum(x);
      daxpy(-tmp, pis, x);
    }
    params.rerror = drerr(prevx, x);
    params.iter += params.steps;
    if (params.rerror < params.rtol) {
      params.info = 0;
      break;
    }
    if (params.iter >= params.maxiter) {
      params.info = -1;
      break;
    }
    callback(params);
  }
}

// template <typename T1, typename MatT>
// List Cmarkovst_gs(const T1& Q, NumericVector x0, int steps, double rtol, int maxiter, MatT) {
//   const int m = marlib::nrow(Q, MatT());
//   const int n = marlib::ncol(Q, MatT());
//   if (m != n) {
//     stop("Matrix Q should be a square matrix.");
//   }
//   if (n != x0.size()) {
//     stop("Vector x0 should be the same dimension of Q.");
//   }
//   NumericVector b(n);
//   NumericVector x = clone(x0);
//   NumericVector prevx(n);
//   int iter = 0;
//   int info = 1;
//   double rerror;
//
//   while(1) {
//     marlib::dcopy(x, prevx);
//     for (int i=0; i<steps; i++) {
//       marlib::gsstep(marlib::TRANS(), 1.0, Q, 0.0, 1.0, b, x, MatT(), marlib::ArrayT());
//       double tmp = marlib::dasum(x);
//       marlib::dscal(1.0/tmp, x);
//     }
//     marlib::daxpy(-1.0, x, prevx);
//     iter += steps;
//     rerror = Rcpp::max(Rcpp::abs(prevx/x));
//     if (rerror < rtol) {
//       info = 0;
//       break;
//     }
//     if (iter >= maxiter) {
//       info = -1;
//       break;
//     }
//     R_CheckUserInterrupt();
//   }
//   return List::create(
//     Named("x")=x,
//     Named("convergence")=(info==0),
//     Named("iter")=iter,
//     Named("rerror")=rerror
//   );
// }

}

#endif

