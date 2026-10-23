/*
  Student t two-sided p-value calculation.

  Historically this used f2c-translated adaptive quadrature code.  All we need
  here is the survival probability for a Student t statistic, which can be
  computed deterministically via the regularized incomplete beta function.
*/

#include "iwpvalue.h"

#include <cmath>
#include <iostream>
#include <limits>

using std::cerr;

namespace {

constexpr int kMaxIterations = 200;
constexpr double kEpsilon = 3.0e-14;
constexpr double kFpMin = std::numeric_limits<double>::min() / kEpsilon;

// Continued fraction for the incomplete beta function. Based on the standard
// Lentz algorithm formulation used in Numerical Recipes.
double
BetaContinuedFraction(double a, double b, double x) {
  const double qab = a + b;
  const double qap = a + 1.0;
  const double qam = a - 1.0;

  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < kFpMin) {
    d = kFpMin;
  }
  d = 1.0 / d;
  double h = d;

  for (int m = 1; m <= kMaxIterations; ++m) {
    const int m2 = 2 * m;

    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < kFpMin) {
      d = kFpMin;
    }
    c = 1.0 + aa / c;
    if (std::fabs(c) < kFpMin) {
      c = kFpMin;
    }
    d = 1.0 / d;
    h *= d * c;

    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < kFpMin) {
      d = kFpMin;
    }
    c = 1.0 + aa / c;
    if (std::fabs(c) < kFpMin) {
      c = kFpMin;
    }
    d = 1.0 / d;
    const double del = d * c;
    h *= del;

    if (std::fabs(del - 1.0) <= kEpsilon) {
      return h;
    }
  }

  cerr << "BetaContinuedFraction:failed to converge for a " << a << " b " << b << " x "
       << x << '\n';
  return h;
}

// Returns I_x(a,b), the regularized incomplete beta function.
double
RegularizedIncompleteBeta(double a, double b, double x) {
  if (x <= 0.0) {
    return 0.0;
  }
  if (x >= 1.0) {
    return 1.0;
  }

  const double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                             a * std::log(x) + b * std::log1p(-x));

  if (x < (a + 1.0) / (a + b + 2.0)) {
    return bt * BetaContinuedFraction(a, b, x) / a;
  }

  return 1.0 - bt * BetaContinuedFraction(b, a, 1.0 - x) / b;
}

}  // namespace

// Return the two-sided p-value for a Student t statistic with `d` degrees of
// freedom. `x` may be signed; only its magnitude matters.
double
iwpvalue(int d, double x) {
  if (d <= 0) {
    cerr << "iwpvalue:invalid degrees of freedom " << d << '\n';
    return 0.0;
  }

  const double t = std::fabs(x);
  if (t == 0.0) {
    return 1.0;
  }

  const double df = static_cast<double>(d);
  const double beta_x = df / (df + t * t);
  return RegularizedIncompleteBeta(0.5 * df, 0.5, beta_x);
}
