#ifndef IW_ACCUMULATOR_H
#define IW_ACCUMULATOR_H

#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#include "kahan_sum.h"

template <typename T, typename SUMMER>
class Accumulator_Base
{
  protected:
    unsigned int _n;
    T _minval;
    T _maxval;
    SUMMER _xsum;
    SUMMER _x2sum;

//  protected functions

    void _reset();
    void _default_values ();

  public:
    Accumulator_Base ();

    int ok () const;

    unsigned int n () const { return _n;}

    T minval () const { return _minval;}
    T maxval () const { return _maxval;}
    T range  () const { return _maxval - _minval;}

    void reset ();

    T sum () const { return _xsum;}
    T sum_of_squares () const { return _x2sum;}

    unsigned int extra (T);
    unsigned int extra (T, int);       // add N copies of a value
    Accumulator_Base<T, SUMMER> & operator = (const Accumulator_Base<T, SUMMER> &);
    void operator += (const Accumulator_Base<T, SUMMER> & rhs) { extra (rhs);}

    unsigned int extra (const T *, int);       // add an array of values
    unsigned int extra (const Accumulator_Base<T, SUMMER> &);

    void operator () (T e) { (void) extra(e);}

    int subtract_data (const Accumulator_Base<T, SUMMER> &, Accumulator_Base<T, SUMMER> &) const;

    double average () const;
    int    average (double &);

    double average_if_available_minval_if_not () const;

    double variance () const;
    int    variance (double &);
};

template <typename T>
class Accumulator : public Accumulator_Base<T, KahanSum>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using Accumulator_Base<T, KahanSum>::_n;
    using Accumulator_Base<T, KahanSum>::_minval;
    using Accumulator_Base<T, KahanSum>::_maxval;
    using Accumulator_Base<T, KahanSum>::_xsum;
    using Accumulator_Base<T, KahanSum>::_x2sum;
#endif

  private:
//  double _xsum;
//  double _x2sum;

//  private functions

    void _default_values ();

  public:

};

template <typename T>
class Accumulator_Int : public Accumulator_Base<T, T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using Accumulator_Base<T, T>::_n;
    using Accumulator_Base<T, T>::_minval;
    using Accumulator_Base<T, T>::_maxval;
    using Accumulator_Base<T, T>::_xsum;
    using Accumulator_Base<T, T>::_x2sum;
#endif

  private:

  public:
};

template <typename T>
std::ostream & operator << (std::ostream &, const Accumulator<T> &);
template <typename T>
std::ostream & operator << (std::ostream &, const Accumulator_Int<T> &);

template <typename T>
class Accumulator_with_Missing_Values : public Accumulator<T>
{
  private:
    int _nmissing;

  public:
    Accumulator_with_Missing_Values ();

    void extra_missing_value () { _nmissing++;}

    int number_missing_values () const { return _nmissing;}

    void reset ();
};

#if defined(ACCUMULATOR_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>

template <typename T, typename SUMMER>
void
Accumulator_Base<T, SUMMER>::_default_values ()
{
  _n = 0;
  _minval = static_cast<T>(0);
  _maxval = static_cast<T>(0);

  _xsum = static_cast<T>(0);
  _x2sum = static_cast<T>(0);

  return;
}

template <typename T, typename SUMMER>
Accumulator_Base<T, SUMMER>::Accumulator_Base ()
{
  _default_values();

  return;
}

template <typename T>
void
Accumulator<T>::_default_values ()
{
  Accumulator_Base<T, KahanSum>::_default_values();

  return;
}

template <typename T, typename SUMMER>
int
Accumulator_Base<T, SUMMER>::ok() const
{
  return 1;
}

template <typename T, typename SUMMER>
void
Accumulator_Base<T, SUMMER>::reset ()
{
  _default_values();

  return;
}

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (T x)
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

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (T x, int n)
{
  if (0 == _n)
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

  _n += n;
  _xsum += n * x;
  _x2sum += n * (x * x);

  return _n;
}

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (const T * x, int n)
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

template <typename T, typename SUMMER>
unsigned int
Accumulator_Base<T, SUMMER>::extra (const Accumulator_Base<T, SUMMER> & rhs)
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

template <typename T, typename SUMMER>
int
Accumulator_Base<T, SUMMER>::subtract_data (const Accumulator_Base<T, SUMMER> & rhs,
                                            Accumulator_Base<T, SUMMER> & zresult) const
{
  assert (_n >= rhs._n);

  zresult._n     = _n - rhs._n;
  zresult._xsum  = _xsum - rhs._xsum;
  zresult._x2sum = _x2sum - rhs._x2sum;

  if (zresult._x2sum >= 0.0)    // exactly how things should be
    ;
  else if (zresult._x2sum > -1.0e-06)     // probably just roundoff errors
  {
    std::cerr << "Accumulator::subtract_data: trimming possible roundoff error " << zresult._x2sum << '\n';
    zresult._x2sum = 0.0;
  }
  else
  {
    std::cerr << "Accumulator::subtract_data: invalid x2sum\n";
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

template <typename T, typename SUMMER>
double
Accumulator_Base<T, SUMMER>::average () const
{
  return static_cast<double>(_xsum) / static_cast<double>(_n);
}

/*template <typename T>
double
Accumulator<T>::variance ()
{       
  assert (_n > 1);

  double tave = Accumulator_Base<T, KahanSum>::average ();
  double rc = _x2sum - _n * tave * tave;

  if (rc < 0.0)   // presumably some roundoff
  {
//  std::cerr << "Accumulator::variance: Warning, negative variance intermediate " << tmp << '\n';

    return 0.0;
  }

  rc = rc / static_cast<double>(_n - 1);

  return rc;
}*/

template <typename T, typename SUMMER>
double
Accumulator_Base<T, SUMMER>::variance () const
{       
  assert (_n > 1);

  if (_n < 2) {
    return 0;
  }

  double tave = Accumulator_Base<T, SUMMER>::average();
  double rc = _x2sum - _n * tave * tave;

  if (rc < 0.0)   // presumably some roundoff
  {
//  std::cerr << "Accumulator::variance: Warning, negative variance intermediate " << tmp << '\n';

    return 0.0;
  }

  rc = rc / static_cast<double>(_n - 1);

  return rc;
}       

template <typename T, typename SUMMER>
int
Accumulator_Base<T, SUMMER>::variance (double & v)
{
  if (_n < 2)
    return 0;

  v = variance();

  return 1;
}

template <typename T>
std::ostream &
operator << (std::ostream & os, const Accumulator<T> & ac)
{
  assert (ac.n() > 0);

  os << "Accumulator " << ac.n() << " values, average " << std::setw(8) << ac.average();
  if (ac.n() > 1)
    os << ", variance " << std::setw(8) << ac.variance();
  os << "\n";

  return os << "min = " << ac.minval() << " max = " << ac.maxval();
}

template <typename T>
std::ostream &
operator << (std::ostream & os, const Accumulator_Int<T> & ac)
{
  assert (ac.n() > 0);

  os << "Accumulator_Int " << ac.n() << " values, average " << ac.average();
  if (ac.n() > 1)
    os << ", variance " << ac.variance();
  os << "\n";

  return os << "min = " << ac.minval() << " max = " << ac.maxval();
}

template <typename T, typename SUMMER>
Accumulator_Base<T, SUMMER> &
Accumulator_Base<T, SUMMER>::operator = (const Accumulator_Base<T, SUMMER> & rhs)
{
  _n = rhs._n;

  _xsum = rhs._xsum;
  _x2sum = rhs._x2sum;
  _minval = rhs._minval;
  _maxval = rhs._maxval;

  return *this;
}

template <typename T, typename SUMMER>
double
Accumulator_Base<T, SUMMER>::average_if_available_minval_if_not () const
{
  if (_n > 1)
    return average();

  if (_n > 0)
    return _minval;

  std::cerr << "Accumulator_Base::average_if_available_minval_if_not: no data!\n";
  return 0.0;
}

#endif

#ifdef ACCUMULATOR_W_MISSING_IMPLEMENTATION

template <typename T>
Accumulator_with_Missing_Values<T>::Accumulator_with_Missing_Values ()
{
  _nmissing = 0;

  return;
}

template <typename T>
void
Accumulator_with_Missing_Values<T>::reset ()
{
  Accumulator<T>::reset();

  _nmissing = 0;

  return;
}

#endif

#ifdef ACCUMULATOR_INT_IMPLEMENTATION
#endif

#endif
