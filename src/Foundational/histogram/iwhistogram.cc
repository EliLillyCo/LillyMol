#include <stdlib.h>

#include <math.h>
#include <limits.h>

#include "Foundational/iwstring/iwstring.h"

#include "iwhistogram.h"

using std::cerr;
using std::endl;

IWHistogram::IWHistogram ()
{
  _min = _max = _dx = 0.0;

  _nsamples = 0;

  _nbuckets = 0;

  _count = nullptr;

  _put_out_of_range_values_in_first_or_last_bucket = 0;

  return;
}

IWHistogram::~IWHistogram ()
{
  if (nullptr != _count)
    delete [] _count;

// set an invalid state

  _min = 0.0;
  _max = -1.0;
  _dx = 5.0;

  return;
}

int
IWHistogram::_initialise ()
{
  _nbuckets = int (rint ((_max - _min) / _dx)) + 1;

  _count = new unsigned int[_nbuckets];

  if (nullptr == _count)
  {
    cerr << "IWHistogram::initialise: cannot allocate array of size " << _nbuckets << endl;
    return 0;
  }

  for (int i = 0; i < _nbuckets; i++)
  {
    _count[i] = 0;
  }

//cerr << "initialise: " << _min << " to " << _max << " nb = " << _nbuckets << endl;

  return 1;
}

#define IWDST_SEPARATOR ','

int
IWHistogram::initialise (const const_IWSubstring & s)
{
  assert (0 == _nsamples);

  if (3 != s.nwords (IWDST_SEPARATOR) || s.starts_with (IWDST_SEPARATOR) || s.ends_with (IWDST_SEPARATOR))
  {
    cerr << "IWHistogram::initialise: initialiser string must be of the form 'min,max,dx', '" <<s << "' is invalid\n";
    return 0;
  }

  const_IWSubstring token;
  int j = 0;
  for (int i = 0; i < 3; i++)
  {
    if (! s.nextword (token, j, IWDST_SEPARATOR))
    {
      cerr << "IWHistogram::initialise: cannot extract token " << i << " from '" << s << "'\n";
      return 0;
    }

    double tmp;
    if (! token.numeric_value (tmp))
    {
      cerr << "IWHistogram::initialise: invalid specifier '" << token << "' in '" << s << "'\n";
      return 0;
    }

    if (0 == i)
      _min = tmp;
    else if (1 == i)
      _max = tmp;
    else if (2 == i)
      _dx = tmp;
    else      // how could this happen
      abort ();
  }

  if (_min >= _max)
  {
    cerr << "IWHistogram::initialise: invalid min max dx " << _min << ',' << _max << ',' << _dx << endl;
    return 0;
  }

  assert (0.0 != _dx);

  return _initialise ();
}

int
IWHistogram::initialise (double mn, double mx, double dx)
{
  assert (0 == _nsamples);

  if (mn >= mx)
  {
    cerr << "IWHistogram::initialise: invalid min and max " << mn << ' ' << mx << endl;
    return 0;
  }

  if (mn + dx > mx)
  {
    cerr << "IWHistogram::initialise: invalid combination " << mn << ' ' << mx << " delta " << dx << endl;
    return 0;
  }

  _min = mn;
  _max = mx;
  _dx = dx;

  return _initialise ();
}

int
IWHistogram::ok () const
{
  if (_nsamples > 0 && nullptr == _count)
    return 0;

  if (_min == _max && _max == _dx)     // probably not initialised
    return 1;

  if (_min >= _max)
    return 0;

  if (_min + _dx > _max)
    return 0;

  return 1;
}

int
IWHistogram::debug_print (std::ostream & os) const
{
  os << "Histogram from " << _min << " to " << _max << " by " << _dx << " with " << _nsamples << " samples\n";

  if (0 == _nsamples)
    return os.good ();

  int populated = 0;
  unsigned int min_pop = UINT_MAX;
  unsigned int max_pop = 0;

  int items = 0;
  for (int i = 0; i < _nbuckets; i++)
  {
    if (_count[i])
    {
      populated++;
      items += _count[i];
    }

    if (_count[i] > max_pop)
      max_pop = _count[i];

    if (_count[i] < min_pop)
      min_pop = _count[i];

    double tmp = _min + i * _dx;
    os << tmp << " " << _count[i] << " values\n";
  }

  os << populated << " of " << _nbuckets << " buckets occupied with " << items << " items. Populations between " << min_pop << " and " << max_pop << endl;

  return os.good ();
}

int
IWHistogram::extra (double e)
{
  int b;
  if (! which_bucket (e, b))
    return 0;

  _count[b]++;
  _nsamples++;

  return 1;
}

int
IWHistogram::which_bucket (double e, int & b) const
{
  if (_min == _max || 0.0 == _dx)
  {
    cerr << "IWHistogram::which_bucket: not initialised\n";
    return 0;
  }

  if (e < _min)
  {
    if (! _put_out_of_range_values_in_first_or_last_bucket)
      return 0; 

    b = 0;

    return 1;
  }

  if (e > _max)
  {
    if (! _put_out_of_range_values_in_first_or_last_bucket)
      return 0;

    b = _nbuckets - 1;

    return 1;
  }

  b = int (rint ((e - _min) / _dx));

  return 1;
}


int
IWHistogram::write (std::ostream & os) const
{
  assert (ok ());

  os << "Histogram with " << _nsamples << " samples\n";

  double cumulative = 0.0;

  for (int i = 0; i < _nbuckets; i++)
  {
    double xstart = _min + i * _dx;
    double fraction;
    if (_nsamples)
    {
      fraction = double (_count[i]) / double (_nsamples);
      cumulative += fraction;
    }
    else
      fraction = 0.0;

    os << _count[i] << " samples (" << static_cast<float> (fraction) << " cumulative " << static_cast<float> (cumulative) << ") between " << xstart << " and " << (xstart + _dx) << endl;
  }

  return os.good ();
}

int
IWHistogram::write_terse(std::ostream & os,
                         int write_zero_values,
                         char output_separator) const
{
  assert (ok ());

  if (0 == _nsamples)
    return os.good ();

  for (int i = 0; i < _nbuckets; i++)
  {
    float xstart = _min + i * _dx;

    if (! write_zero_values && 0 == _count[i])
      ;
    else
      os << xstart << output_separator << _count[i] << endl;
  }

  return os.good ();
}

int
IWHistogram::reset ()
{
  assert (ok ());

  if (0 == _nsamples)
    return 1;

  if (nullptr == _count)
    return 1;

  for (int i = 0; i < _nbuckets; i++)
  {
    _count[i] = 0;
  }

  return 1;
}

IWHistogram &
IWHistogram::add (const IWHistogram & rhs)
{
  assert (rhs._min == _min);
  assert (rhs._max == _max);
  assert (rhs._dx == _dx);

  for (int i = 0; i < _nbuckets; i++)
  {
    _count[i] += rhs._count[i];
  }

  _nsamples += rhs._nsamples;

  return *this;
}

int
IWHistogram::number_samples_less_than (double threshold) const
{
  int nstop;
  if (! which_bucket (threshold, nstop))
  {
    cerr << "IWHistogram::number_samples_less_than: " << threshold << " is out of range\n";
    return 0;
  }

  int rc = 0;
  for (int i = 0; i < nstop; i++)
  {
    rc += _count[i];
  }

  return rc;
}

int
IWHistogram::highest_filled_bucket () const
{
  int rc = -1;

  for (int i = 0; i < _nbuckets; i++)
  {
    if (_count[i] > 0)
      rc = i;
  }

  return rc;
}
