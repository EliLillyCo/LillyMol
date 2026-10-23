#ifndef REPORT_PROGRESS_H
#define REPORT_PROGRESS_H

#include <iostream>

#include <sys/types.h>

/*
  Many programmes have a need to report their progress.
  This provides a compact means of embedding that capability
  into a programme
*/

namespace report_progress_internal {
  void DisplayHelpMessage(std::ostream& output);
}

class Command_Line;

template <typename T>
class Report_Progress_Template {
 private:
  T _times_called;
  T _report_every;
  T _report_next;
  time_t _tzero;
  time_t _tprev;

  // If set, write counts with commas.
  bool _human;

 public:
  Report_Progress_Template();

  template <typename C>
  int initialise(C& cl, char, int verbose);

  int operator()();

  T
  times_called() const {
    return _times_called;
  }

  int
  active() const {
    return _report_every > 0;
  }

  bool human() const {
    return _human;
  }

  T
  report_every() const {
    return _report_every;
  }

  // In two phases of a programme, the first phase might read some molecules
  // into a hash, and the progress reporter would be used then. In the
  // second phase, the molecules in the hash get processed, and we want to
  // keep track of progress there.
  void
  reset() {
    _times_called = 0;
  }

  int set_report_every(T);

  int report(const char* leading, const char* trailing, std::ostream&);

  void initialise(const Report_Progress_Template<T>& rhs);
};

typedef Report_Progress_Template<uint64_t> Report_Progress;

#ifdef REPORT_PROGRESS_IMPLEMENTATION

#include <limits>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/string_with_commas.h"

template <typename T>
Report_Progress_Template<T>::Report_Progress_Template() {
  _times_called = 0;
  _report_every = 0;
  _report_next = std::numeric_limits<T>::max();
  _tzero = static_cast<time_t>(0);
  _tprev = static_cast<time_t>(0);
  _human = false;

  return;
}

/*
  Return 1 if it is time to report
*/

template <typename T>
int
Report_Progress_Template<T>::operator()() {
  if (0 == _report_every) {
    return 0;
  }

  _times_called++;

  if (_times_called <= _report_next) [[likely]] {
    return 0;
  }

  _report_next += _report_every;

  return 1;
}

template <typename T>
int
Report_Progress_Template<T>::report(const char* leading, const char* trailing,
                                    std::ostream& output) {
  if (!operator()()) {
    return 0;
  }

  if (nullptr != leading) {
    output << leading;
  }

  if (_human) {
    output << iwstring::with_commas(_times_called);
  } else {
    output << _times_called;
  }

  if (0 != _tzero) {
    time_t tnow = time(NULL);
    output << " t=" << (tnow - _tzero) << " (" << (tnow - _tprev) << ")";
    _tprev = tnow;
  }

  if (nullptr != trailing) {
    output << trailing;
  }

  return 1;
}

template <typename T>
template <typename C>
int
Report_Progress_Template<T>::initialise(C& cl, char flag, int verbose) {
  const_IWSubstring s;

  for (int i = 0; cl.value(flag, s, i); ++i) {
    if ("time" == s) {
      _tzero = time(NULL);
      _tprev = _tzero;
    } else if (s == "human") {
      _human = true;
    } else if (s == "help") {
      report_progress_internal::DisplayHelpMessage(std::cerr);
      ::exit(0);
    } else if (!cl.value(flag, _report_every) || 0 == _report_every) {
      std::cerr << "Report_Progress::initialise:the report every option (-" << flag
                << ") must be a whole +ve number\n";
      return 0;
    }
  }

  _report_next = _report_every;

  if (verbose) {
    std::cerr << "Will report progress every " << _report_every << " items\n";
  }

  return 1;
}

template <typename T>
int
Report_Progress_Template<T>::set_report_every(T s) {
  _report_every = s;

  if (0 == _times_called) {
    _report_next = s;
  } else {
    _report_next = s + s;

    while (_report_next <= _times_called) {
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
Report_Progress_Template<T>::initialise(const Report_Progress_Template<T>& rhs) {
  _times_called = 0;
  _report_every = rhs._report_every;
  _report_next = _report_every;

  _tzero = static_cast<time_t>(0);
  _tprev = static_cast<time_t>(0);

  return;
}

#endif

#endif
