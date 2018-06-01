#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <iostream>

// Taken from http://formulapages.com/doc/Generating_Normal_Distributed_Random_Numbers

#include "iwrandom.h"

IW_Box_Muller_Normally_Distributed::IW_Box_Muller_Normally_Distributed()
{
  _n2_valid = 0;

  return;
}

IW_Box_Muller_Normally_Distributed::IW_Box_Muller_Normally_Distributed(random_number_seed_t s)
                : MTRand_int32(s)
{
  _n2_valid = 0;

  return;
}

#define REJECTION_SAMPLING_METHOD
//#define TRIG_FUNCTON_FORM
#if defined(REJECTION_SAMPLING_METHOD)

double
IW_Box_Muller_Normally_Distributed::operator() ()
{
  if (_n2_valid)
  {
    _n2_valid = 0;
    return _n2;
  }

// Must generate two numbers such that u1**2 + u2**2 < 1

  double u1, u2, r;

  do
  {
    u1 = 2.0 * MTRand_int32::closed_closed() - 1.0;
    u2 = 2.0 * MTRand_int32::closed_closed() - 1.0;

    r = (u1 * u1 + u2 * u2);
  }
  while (r >= 1.0);

  double t = sqrt(-2.0 * log(r) / r);

  _n2 = u2 * t;
  _n2_valid = 1;

  return u1 * t;
}

#elif defined(TRIG_FUNCTON_FORM)

double
IW_Box_Muller_Normally_Distributed::operator() ()
{
  if (_n2_valid)
  {
    _n2_valid = 0;
    return _n2;
  }

// Must generate two numbers such that u1**2 + u2**2 < 1

  double u1, u2, r;

  do
  {
    u1 = MTRand_int32::open_open();
    u2 = MTRand_int32::open_open();

    r = (u1 * u1 + u2 * u2);
  }
  while (r >= 1.0);

  _n2 = sqrt(-2.0 * log(u1)) * cos(M_PI * 2.0 * u2);
  _n2_valid = 1;

  return sqrt(-2.0 * log(u1)) * sin(M_PI * 2.0 * u2);
}

#else

This actually does not work properly

double
IW_Box_Muller_Normally_Distributed::operator() ()
{
  if (_n2_valid)
  {
    _n2_valid = 0;
    return _n2;
  }

  double u1 = 2.0 * MTRand_int32::open_open() - 1.0;

  double maxu2 = sqrt(1.0 - u1*u1);

  double u2 = 2.0 * maxu2 * MTRand_int32::open_open() - maxu2;
//double u2 = maxu2 * (2.0 * MTRand_int32::open_open() - 1.0);

  double r = u1*u1 + u2 * u2;

  if (r >= 1.0)
  {
    std::cerr << "Out of ragne " << u1 << " and " << u2 << " r = " << r << "\n";
  }

  double t = sqrt(-2.0 * log(r) / r);

  _n2 = u2 * t;
  _n2_valid = 1;

  return u1 * t;
}

#endif


IW_Inverse_Normal_Normally_Distributed::IW_Inverse_Normal_Normally_Distributed()
{
  return;
}

/*
   Hastings approximation

   Feb 2017. Not sure this works...

   Take a look at, which looks better...


   Voutier
   http://arxiv.org/abs/1002.0567v2
*/

double
IW_Inverse_Normal_Normally_Distributed::operator() ()
{
  double x = MTRand_int32::open_open();

  double nx = 1.0 / sqrt(2.0 * M_PI) * exp(-x*x/2.0);

  double t = 1.0 / (1.0 + 0.33267 * x);

  return 1.0 - nx * (t * (0.4361836 + t * (-0.1201676 + t * 0.9372980)));
}

/*
*/

