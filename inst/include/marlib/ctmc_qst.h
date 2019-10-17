#ifndef MARLIB_CTMC_QST_H
#define MARLIB_CTMC_QST_H

namespace marlib {

template <typename T1, typename T2, typename Params, typename Func, typename MatT>
double ctmc_qst_gs(const T1& Q, const T2& xi, T2& x, Params& params, const Func& callback, MatT) {
  using traits1 = vector_traits<T2>;
  int n = traits1::size(x);
  std::vector<double> b(n, 0.0);
  std::vector<double> prevx(n);
  double gam;
  params.iter = 0;
  params.info = 1;
  while(1) {
    dcopy(x, prevx);
    for (int i=0; i<params.steps; i++) {
      gam = ddot(x, xi);
      gsstep(TRANS(), 1.0, Q, -gam, 1.0, b, x, MatT(), ArrayT());
      double tmp = dasum(x);
      dscal(1.0/tmp, x);
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
  return gam;
}

}

#endif

