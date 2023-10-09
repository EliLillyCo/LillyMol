#include <stdlib.h>

#include "Foundational/iwmisc/misc.h"

#include "distance_scaling.h"

using std::cerr;
using std::endl;

Distance_Scaling::Distance_Scaling()
{
  _conv = nullptr;

  return;
}

Distance_Scaling::~Distance_Scaling()
{
  if (nullptr != _conv)
    delete [] _conv;
}

int
Distance_Scaling::build (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Distance_Scaling::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

int
Distance_Scaling::build (iwstring_data_source & input)
{
  if (nullptr == _conv)
  {
    _conv = new float[1001];
    set_vector(_conv, 1001, static_cast<float>(-1.0));
  }

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;

    if (! _parse_scaling_record (buffer))
    {
      cerr << "Distance_Scaling::build:invalid data '" << buffer << "'\n";
      return 0;
    }
  }
  
  int rc = 1;
  for (int i = 0; i < 1001; i++)
  {
    if (_conv[i] < static_cast<float>(0.0))
    {
      cerr << "Distance_Scaling::build:negative scaled distance, i = " << i << " value " << _conv[i] << endl;
      rc = 0;
    }
  }

  if (0 == rc)
    return 0;

  for (int i = 1; i < 1001; i++)
  {
    if (_conv[i] < _conv[i - 1])
    {
      cerr << "Distance_Scaling::build:distances out of order\n";
      cerr << " i = " << (i - 1) << " dist " << _conv[i - 1] << endl;
      cerr << " i = " << i << " dist " << _conv[i] << endl;
      rc = 0;
    }
  }

  return rc;
}

int
Distance_Scaling::_parse_scaling_record (const const_IWSubstring & buffer)
{
  int i = 0;
  const_IWSubstring token;

  if (! buffer.nextword(token, i))
    return 0;

  int ndx;

  if (! token.numeric_value(ndx) || ndx < 0 || ndx > 1000)
  {
    cerr << "Distance_Scaling::_parse_scaling_record:invalid index '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "Distance_Scaling::_parse_scaling_record:not enough tokens\n";
    return 0;
  }

  if (! token.numeric_value(_conv[ndx]) || _conv[ndx] < static_cast<float>(0.0))
  {
    cerr << "Distance_Scaling::_parse_scaling_record:negative distance, i = " << i << " '" << token << "'\n";
    return 0;
  }

  return 1;
}

float
Distance_Scaling::convert (float v) const
{
  int ndx = static_cast<int>(v * static_cast<float>(1000.0) + 0.49999);

  if (ndx < 0 || ndx > 1000)
    cerr << "Yipes, index for distance " << v << " out of range, " << ndx << endl;
  assert (ndx >= 0 && ndx <= 1000);

  return _conv[ndx];
}
