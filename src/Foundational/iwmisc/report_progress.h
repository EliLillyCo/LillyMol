#ifndef REPORT_PROGRESS_H
#define REPORT_PROGRESS_H

#include <iostream>

#include <sys/types.h>

/*
  Many programmes have a need to report their progress.
  This provides a compact means of embedding that capability
  into a programme
*/

class Command_Line;

template <typename T>
class Report_Progress_Template
{
  private:
    T _times_called;
    T _report_every;
    T _report_next;
    time_t _tzero;
    time_t _tprev;

  public:
    Report_Progress_Template ();

    template <typename C> int initialise (C & cl, char, int verbose);

    int operator ()();

    T times_called() const { return _times_called;}

    int active () const { return _report_every > 0;}

    int set_report_every(T);

    int report (const char * leading, const char * trailing, std::ostream &);

    void initialise (const Report_Progress_Template<T> & rhs);
};

typedef Report_Progress_Template<unsigned int> Report_Progress;

#ifdef REPORT_PROGRESS_IMPLEMENTATION

#include <limits>

#include "cmdline.h"

template <typename T>
Report_Progress_Template<T>::Report_Progress_Template()
{
  _times_called = 0;
  _report_every = 0;
  _report_next = std::numeric_limits<T>::max();
  _tzero = static_cast<time_t>(0);
  _tprev = static_cast<time_t>(0);

  return;
}

/*
  Return 1 if it is time to report
*/

template <typename T>
int
Report_Progress_Template<T>::operator ()()
{
  if (0 == _report_every)
    return 0;

  _times_called++;

  if (_times_called <= _report_next)
    return 0;

  _report_next += _report_every;

  return 1;
}

template<typename T>
int
Report_Progress_Template<T>::report (const char * leading,
                                     const char * trailing,
                                     std::ostream & output)
{
  if (! operator()())
    return 0;

  if (NULL != leading)
    output << leading;

  output << _times_called;

  if (0 != _tzero)
  {
    time_t tnow = time(NULL);
    output << " t=" << (tnow - _tzero) << " (" << (tnow - _tprev) << ")";
    _tprev = tnow;
  }

  if (NULL != trailing)
    output << trailing;

  return 1;
}

template <typename T> template<typename C>
int
Report_Progress_Template<T>::initialise (C & cl, char flag, int verbose)
{
  const_IWSubstring s;

  for (int i = 0; cl.value(flag, s, i); ++i)
  {
    if ("time" == s)
    {
      _tzero = time(NULL);
      _tprev = _tzero;
    }
    else if (! cl.value(flag, _report_every) || 0 == _report_every)
    {
      cerr << "Report_Progress::initialise:the report every option (-" << flag << ") must be a whole +ve number\n";
      return 0;
    }
  }

  _report_next = _report_every;

  if (verbose)
    cerr << "Will report progress every " << _report_every << " items\n";

  return 1;
}

template <typename T>
int
Report_Progress_Template<T>::set_report_every (T s)
{
  _report_every = s;

  if (0 == _times_called)
    _report_next = s;
  else
  {
    _report_next = s + s;

    while (_report_next <= _times_called)
    {
      _report_next += s;
    }
  }

  return 1;
}

/*
  We do not transfer all the properties, 
*/

template <typename T>
void
Report_Progress_Template<T>::initialise (const Report_Progress_Template<T> & rhs)
{
  _times_called = 0;
  _report_every = rhs._report_every;
  _report_next = _report_every;

  _tzero = static_cast<time_t>(0);
  _tprev = static_cast<time_t>(0);

  return;
}

#endif

#endif
