#ifndef IW_ACCUMULATOR_IMPLEMENTATION
#define IW_ACCUMULATOR_IMPLEMENTATION

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

#include "accumulator.h"

template <typename T>
void
Accumulator<T>::_default_values ()
{
  _n = 0;
  _xsum = 0.0;
  _x2sum = 0.0;
  _minval = 0;
  _maxval = 0;
}

template <typename T>
Accumulator<T>::Accumulator ()
{
  _default_values ();

  return;
}

template <typename T>
Accumulator<T>::~Accumulator ()
{
  if (-87 == _n)
    cerr << "Freeing already deleted Accumulator\n";

  _n = -87;

  return;
}

template <typename T>
int
Accumulator<T>::ok () const
{
  if (_n < 0)
    return 0;

  return 1;
}

template <typename T>
void
Accumulator<T>::reset ()
{
  _default_values ();

  return;
}

template <typename T>
int
Accumulator<T>::extra (T x)
{
  _n++;
  _xsum += x;
  _x2sum += (x * x);

  if (1 == _n)
  {
    _minval = x;
    _maxval = x;
  }
  else if (x < _minval)
  {
    _minval = x;
  }
  else if (x > _maxval)
  {
    _maxval = x;
  }

  return _n;
}

/*
  Add N copies of a value
*/

template <typename T>
int
Accumulator<T>::extra (T x, int n)
{
  _n += n;
  _xsum += n * x;
  _x2sum += n * (x * x);

  if (1 == _n)
  {
    _minval = x;
    _maxval = x;
  }
  else if (x < _minval)
  {
    _minval = x;
  }
  else if (x > _maxval)
  {
    _maxval = x;
  }

  return _n;
}

template <typename T>
int
Accumulator<T>::extra (const T * x, int n)
{
  if (0 == _n)
  {
    _minval = x[0];
    _maxval = x[1];
  }

  for (int i = 0; i < n; i++)
  {
    if (x[i] < _minval)
      _minval = x[i];
    else if (x[i] > _maxval)
      _maxval = x[i];

    _xsum += x[i];
    _x2sum += x[i] * x[i];
  }

  _n += n;

  return _n;
}

template <typename T>
int
Accumulator<T>::extra (const Accumulator<T> & rhs)
{
  if (0 == rhs._n)
    return _n;

  if (0 == _n)
  {
    _n = rhs._n;
    _xsum  = rhs._xsum;
    _x2sum = rhs._x2sum;
    _minval = rhs._minval;
    _maxval = rhs._maxval;

    return _n;
  }

  _n += rhs._n;
  _xsum += rhs._xsum;
  _x2sum += rhs._x2sum;

  if (rhs._minval < _minval)
    _minval = rhs._minval;
  if (rhs._maxval > _maxval)
    _maxval = rhs._maxval;

  return _n;
}

/*
  Strictly speaking, _minval and _maxval are unknown, but we lie...
*/

template <typename T>
int
Accumulator<T>::subtract_data (const Accumulator<T> & rhs,
                               Accumulator<T> & zresult) const
{
  assert (_n >= rhs._n);

  zresult._n     = _n - rhs._n;
  zresult._xsum  = _xsum - rhs._xsum;
  zresult._x2sum = _x2sum - rhs._x2sum;

  if (zresult._x2sum >= 0.0)    // exactly how things should be
    ;
  else if (zresult._x2sum > -1.0e-06)     // probably just roundoff errors
  {
    cerr << "Accumulator::subtract_data: trimming possible roundoff error " << zresult._x2sum << endl;
    zresult._x2sum = 0.0;
  }
  else
  {
    cerr << "Accumulator::subtract_data: invalid x2sum\n";
    abort ();
  }

  if (_minval < rhs._minval)
    zresult._minval = _minval;
  else
    zresult._minval = rhs._minval;

  if (_maxval > rhs._maxval)
    zresult._maxval = _maxval;
  else
    zresult._maxval = rhs._maxval;

  return 1;
}

template <typename T>
double
Accumulator<T>::average () const
{
  assert (_n > 0);

  return static_cast<double> (_xsum) / static_cast<double> (_n);
}

template <typename T>
double
Accumulator<T>::variance ()
{       
  assert (_n > 1);

  double tave = average ();
  double rc = _x2sum - _n * tave * tave;

  if (rc < 0.0)   // presumably some roundoff
  {
//  cerr << "Accumulator::variance: Warning, negative variance intermediate " << tmp << endl;

    return 0.0;
  }

  rc = rc / static_cast<double> (_n - 1);

  return rc;
}       

template <typename T>
double
Accumulator<T>::variance () const
{       
  assert (_n > 1);

  double tave = average ();
  double rc = _x2sum - _n * tave * tave;

  if (rc < 0.0)   // presumably some roundoff
  {
//  cerr << "Accumulator::variance: Warning, negative variance intermediate " << tmp << endl;

    return 0.0;
  }

  rc = rc / static_cast<double> (_n - 1);

  return rc;
}       

template <typename T>
int
Accumulator<T>::variance (double & v)
{
  if (_n < 2)
    return 0;

  v = variance ();

  return 1;
}

template <typename T>
ostream &
operator << (ostream & os, const Accumulator<T> & ac)
{
  assert (ac.n () > 0);

  os << "Accumulator " << ac.n () << " values, average " << setw(8) << ac.average ();
  os << ", variance " << setw (8) << ac.variance () << ",\n";
  return os << "min = " << ac.minval () << " max = " << ac.maxval ();
}

template <typename T>
Accumulator<T> &
Accumulator<T>::operator = (const Accumulator<T> & rhs)
{
  _n = rhs._n;

  _xsum = rhs._xsum;
  _x2sum = rhs._x2sum;
  _minval = rhs._minval;
  _maxval = rhs._maxval;

  return *this;
}

template <typename T>
double
Accumulator<T>::average_if_available_minval_if_not () const
{
  if (_n > 1)
    return average ();

  if (_n > 0)
    return _minval;

  cerr << "Accumulator::average_if_available_minval_if_not: no data!\n";
  return 0.0;
}

#endif
