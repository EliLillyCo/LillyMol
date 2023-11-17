/*
  Define the function to be integrated for the P distribution
*/

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <assert.h>

#include "gsl/gsl_integration.h"

using namespace std;

#include "iwpvalue.h"

static int dof = 0;
static double double_dof = 0.0;
static double dof_plus_one_divided_by_two = 0.0;
static double prefactor = 0.0;

/*
  In order to avoid dealing with potentially large numbers (but at the expense of numerical
  accuracy, we can compute the ratio of the gamma functions as a series of ratios
*/

static double
gamma_ratio_even (int n)
{
  double rc = 0.50 * sqrt (M_PI);

  for (int i = n - 1; i > 1; i -= 2)
  {
    rc = rc * static_cast<double> (i) / static_cast<double> (i - 1);
  }

  return rc;
}

static double
gamma_ratio_odd (int n)
{
  if (1 == n)
    return 1.0 / sqrt (M_PI);

  double rc = 2.0 / sqrt (M_PI);

  for (int i = n - 1; i > 2; i -= 2)
  {
    rc = rc * static_cast<double> (i) / static_cast<double> (i - 1);
  }

  return rc;
}

static double
gamma_ratio (int n)
{
  if (n == (n / 2) * 2)
    return gamma_ratio_even (n);
  else
    return gamma_ratio_odd (n);
}

/*
  Return the gamma function of n + 0.5
*/

static double
gamma_half (int n)
{
  int n2m1 = 2 * n - 1;

  double numerator = 1.0;
  for (int i = 3; i <= n2m1; i += 2)
  {
    numerator = numerator * static_cast<double> (i);
  }

  double rc = numerator * sqrt (M_PI) / pow (2.0, static_cast<double> (n));

//cerr << "Gamma (" << n << " + 0.5) = " << rc << " numerator " << numerator << endl;

  return rc;
}

void
set_pvalue_degrees_of_freedom (int d)
{
  dof = d;

  double_dof = static_cast<double> (d);

  dof_plus_one_divided_by_two = static_cast<double> (d + 1) / 2.0;

  prefactor = gamma_ratio (d) / sqrt (dof * M_PI);

//cerr << "Dof = " << dof << " numerator " << numerator << " denominator " << denominator << " prefactor " << prefactor << endl;

  return;
}

double
for_p_distribution (double dx, void * v)
{
  assert (dof > 0);

  double d1 = 1.0 + dx * dx / double_dof;

  double rc = prefactor / pow (d1, dof_plus_one_divided_by_two);

//cerr << " x = " << (*x) << " value " << rc << endl;
  return rc;
}

//    SUBROUTINE Q1DA(F,A,B,EPS,R,E,KF,IFLAG)

extern "C" void dq1da_ (double (*) (const double *), 
                  const double *,      // A
                  const double *,      // B
                  const double *,      // EPS
                  double *,            // RESULT
                  double *,            // ABSERR
                  int *,              // NEVAL
                  int *);             // IFLAG

double
iwpvalue (int d, double x)
{
  if (d != dof)
    set_pvalue_degrees_of_freedom (d);

  double a = 0.0;
  double rc;
  double epsabs = 1.0e-01;
  double epsrel = 1.0e-02;
  double abserr;
  size_t neval;

  static gsl_function * f = nullptr;

  if (NULL == f)
  {
    f = new gsl_function;
    f->function = &for_p_distribution;
    f->params = nullptr;
  }

  gsl_integration_qng (f, a, x, epsabs, epsrel, &rc, &abserr, &neval);

  cerr << "Error " << abserr << " neval " << neval <<endl;

// The value from -inf to zero is always 0.5

  rc += 0.50;

  rc = 2.0 * (1.0 - rc);

  return rc;
}
