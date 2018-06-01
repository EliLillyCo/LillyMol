#include <stdlib.h>
#include <limits>

#include "cmdline.h"
#include "report_progress.h"

Report_Progress::Report_Progress()
{
  _times_called = 0;
  _report_every = 0;
  _report_next = std::numeric_limits<unsigned int>::max();

  return;
}

/*
  Return 1 if it is time to report
*/

int
Report_Progress::operator ()()
{
  if (0 == _report_every)
    return 0;

  _times_called++;

  if (_times_called <= _report_next)
    return 0;

  _report_next += _report_every;

  return 1;
}

template <typename C>
int
Report_Progress::initialise (C & cl, char flag, int verbose)
{
  if (! cl.value(flag, _report_every) || 0 == _report_every)
  {
    cerr << "Report_Progress::initialise:the report every option (-" << flag << ") must be a whole +ve number\n";
    return 0;
  }

  _report_next = _report_every;

  if (verbose)
    cerr << "Will report progress every " << _report_every << " items\n";

  return 1;
}

int
Report_Progress::set_report_every (int s)
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
