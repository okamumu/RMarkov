#ifndef MRLIB_MEXP_MIX_H
#define MRLIB_MEXP_MIX_H

namespace marlib {

  template <typename TR, typename T1, typename T2, typename T3, typename T4,
            typename Params, typename Func1, typename Func2,
            typename MatT, typename VecT>
  void mexp_mix(TR, T1& A, T2& x, T2& y, const T3& w, const T4& t,
                Params& params, const Func1& callback1, const Func2& callback2,
                MatT, VecT) {
    using traits1 = vector_traits<T3>;
    using traits2 = vector_traits<T4>;
    const int m = traits1::size(w);
    const double* valueW = traits1::value(w);
    const double* valueT = traits2::value(t);
    const int incw = traits1::inc(w);
    const int inct = traits2::inc(t);
    
    marlib::dfill(y, 0.0);
    double qv = unif(A, params.ufact, MatT());
    double maxt = dmax(t);
    params.r = poi::rightbound(qv*maxt, params.eps);
    if (params.r > params.rmax) {
      callback1(params); //stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
    }
    std::vector<double> prob(params.r+1);
    for (int i=0; i<m; i++, valueW+=incw, valueT+=inct) {
      double w = *valueW;
      double t = *valueT;
      int r = poi::rightbound(qv*t, params.eps);
      double weight = poi::pmf(qv*t, 0, r, prob);
      mexp(TR(), A, prob, r, weight, x, x, MatT(), VecT());
      daxpy(w, x, y);
      callback2(params); // R_CheckUserInterrupt();
    }
  }

  template <typename TR, typename T1, typename T2, typename T3, typename T4,
            typename Params, typename Func1, typename Func2,
            typename MatT, typename VecT>
  void mexpint_mix(TR, T1& A, T2& x, T2& cx, T2& y, T2& cy, const T3& w, const T4& t,
                   Params& params, const Func1& callback1, const Func2& callback2,
                   MatT, VecT) {
    using traits1 = vector_traits<T3>;
    using traits2 = vector_traits<T4>;
    const int m = traits1::size(w);
    const double* valueW = traits1::value(w);
    const double* valueT = traits2::value(t);
    const int incw = traits1::inc(w);
    const int inct = traits2::inc(t);
    
    marlib::dfill(cx, 0.0);
    marlib::dfill(y, 0.0);
    marlib::dfill(cy, 0.0);
    double qv = unif(A, params.ufact, MatT());
    double maxt = dmax(t);
    params.r = poi::rightbound(qv*maxt, params.eps) + 1;
    if (params.r > params.rmax) {
      callback1(params); //stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
    }
    std::vector<double> prob(params.r+1);
    std::vector<double> cprob(params.r+1);
    for (int i=0; i<m; i++, valueW+=incw, valueT+=inct) {
      double w = *valueW;
      double t = *valueT;
      int r = poi::rightbound(qv*t, params.eps) + 1;
      double weight = poi::cpmf(qv*t, 0, r, prob, cprob);
      mexpint(TR(), A, prob, cprob, r, weight, qv*weight, x, x, cx, MatT(), VecT());
      daxpy(w, x, y);
      daxpy(w, cx, cy);
      callback2(params); // R_CheckUserInterrupt();
    }
  }
}


#endif

