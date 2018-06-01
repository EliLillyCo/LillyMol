#ifndef iwMINMAX_H
#define iwMINMAX_H

#if (__GNUC__ == 3) && (__GNUC_MINOR__ > 3)
#define IW_TWO_PHASE_TEMPLATES
#endif

#if (__GNUC__ > 3)
#define IW_TWO_PHASE_TEMPLATES
#endif

/*
  Compares a series of values and keeps track of the minimum.
  NOTE, we only return 1 when we encounter a value LESS THAN our existing
  value.
*/

template <typename T>
class iwmin
{
  protected:
    T _value;
  public:
    iwmin (T x) { _value = x;}

    void initialise (T x) { _value = x;}

    void extra (T x) { if (x < _value) _value = x;}
    int  try_this (T x) { if (x < _value) { _value = x; return 1;} else return 0;}

    T minval () const { return _value;}
};

template <typename T>
class iwmax
{
  protected:
    T _value;

  public:
    iwmax (T x) { _value = x;}

    void initialise (T x) {_value = x;}

    void extra (T x) { if (x > _value) _value = x;}
    int  try_this (T x) { if (x > _value) { _value = x; return 1;} else return 0;}

    T maxval () const { return _value;}
};

template <typename T>
class iwminmax
{
  private:
    T _minval;
    T _maxval;
  public:
    iwminmax (T, T);

    void extra (T x) { if (x < _minval) _minval = x;else if (x > _maxval) _maxval = x;}

    int try_this (T x) { if (x < _minval)\
                           {_minval = x; return -1;}\
                         else if (x > _maxval) \
                           {_maxval = x; return 1;}\
                         return 0;}

    T minval () const { return _minval;}
    T maxval () const { return _maxval;}
};

/*
  Often we need to keep track of the minimum or maximum, together with
  some other piece of data - the identity of the item which has the
  min or max value;

  These need a value in the constructor. But I ran into problems
  when I really didn't know what the range of values would be at
  the time of the constructor, so there is a SET () method
*/

template <typename T, typename W>
class iwminid : public iwmin<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  using iwmin<T>::_value;
#endif

  private:
    W _id_of_min;

  public:
    iwminid (T, W);

    void set (T t, W w) { _value = t; _id_of_min = w;}

    int try_this (T x, W w);

    W which_is_min () const { return _id_of_min;}
};

template <typename T, typename W>
class iwmaxid : public iwmax<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  using iwmax<T>::_value;
#endif

  private:
    W _id_of_max;

  public:
    iwmaxid (T, W);

    void set (T t, W w) { _value = t; _id_of_max = w;}

    int try_this (T x, W w);

    W which_is_max () const { return _id_of_max;}
};

template <typename T>
class iwmincount : public iwmin<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  using iwmin<T>::_value;
#endif

  private:
    int _count_of_number_of_occurrences_of_min;

  public:
    iwmincount () {_count_of_number_of_occurrences_of_min = 0;}

    int try_this (T);

    int how_many_times_did_min_occur () const { return _count_of_number_of_occurrences_of_min;}
};

template <typename T>
class iwmaxcount : public iwmax<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  using iwmax<T>::_value;
#endif

  private:
    int _count_of_number_of_occurrences_of_max;

  public:
    iwmaxcount () {_count_of_number_of_occurrences_of_max = 0;}

    int try_this (T);

    int how_many_times_did_max_occur () const { return _count_of_number_of_occurrences_of_max;}
};

#ifdef IWMINMAX_IMPLEMENTATION

#include <assert.h>

template <typename T>
iwminmax<T>::iwminmax (T init_low, T init_high)
{
  _minval = init_low;
  _maxval = init_high;
}


template <typename T, typename W>
iwminid<T, W>::iwminid (T init_low, W id_low) : iwmin<T> (init_low)
{
  _id_of_min = id_low;
}

template <typename T, typename W>
iwmaxid<T, W>::iwmaxid (T init_low, W id_low) : iwmax<T> (init_low)
{
  _id_of_max = id_low;
}

template <typename T, typename W>
int
iwminid<T, W>::try_this (T x, W w)
{
  if (iwmin<T>::try_this (x))
  {
    _id_of_min = w;
    return 1;
  }
  else
    return 0;
}

template <typename T, typename W>
int
iwmaxid<T, W>::try_this (T x, W w)
{
  if (iwmax<T>::try_this (x))
  {
    _id_of_max = w;
    return 1;
  }
  else
    return 0;
}

/*
  This relies on the fact that iwmin::try_this will return 1 only for
  values which are actually lower than its current value. We then intercept
  any which are equal to what we already have.
*/

template <typename T>
int
iwmincount<T>::try_this (T x)
{
  if (iwmin<T>::try_this (x))     // got a new min value
  {
    _count_of_number_of_occurrences_of_min = 1;
    return 0;
  }
  else if (x == _value)
  {
    _count_of_number_of_occurrences_of_min++;
    return 1;                     // note that we return 1 on an equal comparison
  }
  else
    return 0;
}

template <typename T>
int
iwmaxcount<T>::try_this (T x)
{
  if (iwmax<T>::try_this (x))     // got a new max value
  {
    _count_of_number_of_occurrences_of_max = 1;
    return 0;
  }
  else if (x == _value)
  {
    _count_of_number_of_occurrences_of_max++;
    return 1;                     // note that we return 1 on an equal comparison
  }
  else
    return 0;
}

#endif
#endif
