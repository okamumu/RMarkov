/**
 * @file poisson.cpp
 * @brief Compute p.m.f. of poisson distribution
 * */

#include <cmath>
#include <cassert>

namespace marlib {

  namespace poi {
  	static const double pi = 4.0 * atan(1.0);
  	static const double NormalQ_LowerQ = 3.0;
  	static const double NormalQ_UpperQ = 37.0;
  	// static const double NormalQ_LowerLogP = -689.0;
  	// static const double NormalQ_UpperLogP = -6.6;
  	static const double NormalQ_eps = 1.0e-8;
  	static const double LOG2PIOVER2 = log(2*pi) / 2.0;
  	static const double POISSON_LAMBDA_MIN = 3.0;
  	static const int POISSON_RIGHT_MAX = 23;

  	double normalt(double x) {
  		double x2 = x * x;
  		double tmp = x;
  		double sum = 1.0 / tmp;
  		tmp *= x2;
  		sum -= 1.0 / tmp;
  		tmp *= x2;
  		sum += 3.0 / tmp;
  		tmp *= x2;
  		sum -= 15.0 / tmp;
  		tmp *= x2;
  		sum += 105.0 / tmp;
  		return (log(sum) - x2/2.0 - LOG2PIOVER2);
  	}

  	double normalq(double p) {
  		const double leps = log(p);
  //		assert(leps <= NORMALQ_UPPER_LOGP && leps >= NORMALQ_LOWER_LOGP);
  		double l = NormalQ_LowerQ;
  		double u = NormalQ_UpperQ;
  		double m = (l + u) / 2;
  		double fm = normalt(m) - leps;
  		while (std::abs(fm) > NormalQ_eps) {
  			if (fm > 0) {
  				l = m;
  			} else {
  				u = m;
  			}
  			m = (l + u)/2;
  			fm = normalt(m) - leps;
  		}
  		return m;
  	}

  	int rightbound(double lambda, double eps) {
  	  if (std::fpclassify(lambda) == FP_ZERO) {
  	    return 0;
  	  }
  	  if (lambda < POISSON_LAMBDA_MIN) {
  	    double tmp = std::exp(-lambda);
  	    double total = tmp;
        int right = 0;
        for (int k=1; k<=POISSON_RIGHT_MAX; k++) {
          right += 1;
          tmp *= lambda / right;
          total += tmp;
          if (total + eps >= 1)
            break;
        }
        return right;
  	  } else {
  	    double z = normalq(eps);
  	    double tmp = z + sqrt(4 * lambda - 1);
  			return static_cast<int>(tmp * tmp / 4 + 1);
  	  }
  	}
  }
}
