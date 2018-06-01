#ifndef IW_RUNNING_AVE_H
#define IW_RUNNING_AVE_H

#include "most_recent.h"

#if (__GNUC__ == 3) && (__GNUC_MINOR__ > 3)
#define IW_TWO_PHASE_TEMPLATES
#endif

template <typename T>
class Running_Average : public IWMost_Recent<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  using IWMost_Recent<T>::_things;
  using IWMost_Recent<T>::_items_added;
  using IWMost_Recent<T>::_elements_allocated;
  using IWMost_Recent<T>::_next_ptr;
#endif

  private:
    T _sum;
    T _sum_of_squares;

// private functions

    void _default_values ();

  public:
    Running_Average ();
    Running_Average (int);
    Running_Average (const Running_Average<T> &);

    int ok () const;
    int debug_print (std::ostream &) const;

    Running_Average & operator = (const Running_Average &);

    void extra (T);

    void clear ();

    int resize (int);

    double average () const;

    double variance () const;
};

#ifdef RUNNING_AVERAGE_IMPLEMENTATION

#include <assert.h>
#include <math.h>

#include "running_average.h"

template <typename T>
void
Running_Average<T>::_default_values ()
{
  _sum = T (0);
  _sum_of_squares = T (0);

  return;
}

template <typename T>
Running_Average<T>::Running_Average ()
{
  _default_values ();

  return;
}

template <typename T>
Running_Average<T>::Running_Average (int initial_size) : IWMost_Recent<T>::IWMost_Recent (initial_size)
{
  _default_values ();

  return;
}

template <typename T>
Running_Average<T>::Running_Average (const Running_Average<T> & rhs)
{
  _default_values ();

  operator = (rhs);

  return;
}

template <typename T>
int
Running_Average<T>::ok () const
{
  if (! IWMost_Recent<T>::ok ())
    return 0;

  if (0 == _items_added)
  {
    return (T (0) == _sum && T (0) == _sum_of_squares);
  }

  return 1;
}

template <typename T>
Running_Average<T> &
Running_Average<T>::operator = (const Running_Average<T> & rhs)
{
  IWMost_Recent<T>::operator = (rhs);

  _sum = rhs._sum;

  _sum_of_squares = rhs._sum_of_squares;

  return *this;
}

template <typename T>
void
Running_Average<T>::clear ()
{
  IWMost_Recent<T>::clear ();

  _sum = T (0);
  _sum_of_squares = T (0);

  return;
}

/*
  The extra function is complicated by the need to keep track of anything
  which is being discarded;
*/

template <typename T>
void
Running_Average<T>::extra (T e)
{
  assert (_elements_allocated > 0);

  if (_items_added >= _elements_allocated)
  {
    T disappear = _things[_next_ptr];

    _sum -= disappear;
    _sum_of_squares -= (disappear * disappear);
  }

  IWMost_Recent<T>::extra (e);

  _sum += e;
  _sum_of_squares += (e * e);

  return;
}

template <typename T>
double
Running_Average<T>::average () const
{
  assert (_items_added > 0);

  double denominator;
  if (_items_added <= _elements_allocated)
    denominator = double (_items_added);
  else
    denominator = double (_elements_allocated);

  return double (_sum) / denominator;
}

template <typename T>
double
Running_Average<T>::variance () const
{
  assert (_items_added > 1);

  double tave = double (_sum) / double (_items_added);     // the average

  double rc = _sum_of_squares - _items_added * tave * tave;

  if (rc < 0.0)   // presumably some roundoff
    return 0.0;

  rc = sqrt (rc / double (_items_added - 1));

  return rc;
}

#endif
#endif
