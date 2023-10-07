#include <iostream>
#include "assert.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Foundational/iwmisc/misc.h"
#include "demerit.h"

using std::cerr;

static int demerit_reason_contains_individual_demerits = 0;

void
set_demerit_reason_contains_individual_demerits (int s)
{
  demerit_reason_contains_individual_demerits = s;

  return;
}

static int store_demerit_reasons_like_tsubstructure = 0;

void
set_store_demerit_reasons_like_tsubstructure (int s)
{
  store_demerit_reasons_like_tsubstructure = s;

  return;
}

static int _rejection_threshold = DEFAULT_REJECTION_THRESHOLD;

void
set_rejection_threshold (int s)
{
  _rejection_threshold = s;
}

int
rejection_threshold()
{
  return _rejection_threshold;
}

Demerit::Demerit()
{
  _score = 0;

  _append_demerit_value_with_reason = 0;

  _atomic_demerit = nullptr;

  _need_to_delete_atomic_demerit = 0;

  return;
}

Demerit::Demerit(Demerit&& rhs) {
  _demerit = std::move(rhs._demerit);
  _score = rhs._score;
  _append_demerit_value_with_reason = rhs._append_demerit_value_with_reason;
  _atomic_demerit = rhs._atomic_demerit;
  _need_to_delete_atomic_demerit = rhs._need_to_delete_atomic_demerit;
  rhs._need_to_delete_atomic_demerit = 0;
}

Demerit::~Demerit()
{
  if (nullptr == _atomic_demerit)
    ;
  else if (_need_to_delete_atomic_demerit)
    delete [] _atomic_demerit;

  return;
}

int
Demerit::initialise_atomic_demerits(const int matoms)
{
  assert (nullptr == _atomic_demerit);

  _atomic_demerit = new_int(matoms);
  _need_to_delete_atomic_demerit = 1;

  return 1;
}

int
Demerit::set_atomic_demerit_array(int * d)
{
  assert (nullptr == _atomic_demerit);

  _atomic_demerit = d;

  _need_to_delete_atomic_demerit = 0;

  return 1;
}

int
Demerit::debug_print (std::ostream & os) const
{
  os << "Demerit: " << _demerit.number_elements() << " demerits, ";
  if (_score >= _rejection_threshold)
    os << "REJECTED";
  else
    os << "total " << _score;

  for (int i = 0; i < _demerit.number_elements(); ++i)
  {
    os << _demerit[i]->demerit() << ' ' << _demerit[i]->reason() << '\n';
  }

  return 1;
}

int
Demerit::rejected () const
{
  return _score >= _rejection_threshold;
}

int
Demerit::_add_demerit(const int increment, const const_IWSubstring & reason)
{
  if (increment <= 0)
  {
    cerr << "Demerit::_increment:zero increment ignored\n";
    return 0;
  }

  _score += increment;

  Demerit_and_Reason * x = new Demerit_and_Reason(increment, reason);

  if (_demerit.empty() || increment >= _demerit.last_item()->demerit())
    return _demerit.add(x);

  for (int i = 0; i < _demerit.number_elements(); ++i)
  {
    if (increment <= _demerit[i]->demerit())
      return _demerit.insert_before(i, x);
  }

  return 1;    // should not come here
}

int
Demerit::extra(int increment, const const_IWSubstring reason)
{
  return _add_demerit(increment, reason);
}

int 
Demerit::reject(const const_IWSubstring reason)
{
  return _add_demerit(_rejection_threshold, reason);
}

void
Demerit::do_sort()
{
  const int n = _demerit.number_elements();

  if (n <= 1)
    return;

  if (_demerit[0]->demerit() < _demerit[1]->demerit())
    _demerit.swap_elements(0, 1);

  if (2 == n)
    return;

  _demerit.iwqsort_lambda([] (const Demerit_and_Reason * dr1, const Demerit_and_Reason * dr2)
  {
    if (dr1->demerit() < dr2->demerit())
      return 1;
    if (dr1->demerit() > dr2->demerit())
      return -1;
    return 0;
  });

  return;
}

#ifdef NO_LONGER_USED
void
Demerit::_add_hit_type(const int increment,
                       const const_IWSubstring & reason)
{
  if (store_demerit_reasons_like_tsubstructure)
  {
    if (_types.length())
      _types.add(' ');

    _types << "(1 matches to 'D" << increment << ' ' << reason << "')";
//  _types << "(1 matches to '" << reason << "')";
  }
  else if (demerit_reason_contains_individual_demerits)
  {
    if (_types.length())
      _types.add(':');

    _types << 'D' << increment << ' ' << reason;
  }
  else if (_append_demerit_value_with_reason)
  {
    _types.append_with_spacer(reason, ':');
    _types << ' ' << increment;
  }
  else
    _types.append_with_spacer(reason, ':');

  return;
}
#endif

int
Demerit::write_in_tdt_form (std::ostream & os) const
{
  if (_score >= _rejection_threshold)
    os << "REJ<1>\n";

  os << "DMRT<" << _score << ">\n";
  os << "NDMRT<" << _demerit.number_elements() << ">\n";
  os << "DMRTYP<";
  for (int i = 0; i < _demerit.number_elements(); ++i)
  {
    if (i > 0)
      os << ':';
    os << _demerit[i]->reason();
  }
  os << ">\n";

  return 1;
}

/*int
Demerit::reject (const IWString & reason)
{
  _increment(_rejection_threshold);

  _add_hit_type(reason);

  return 1;
}

int 
Demerit::reject (const char * reason)
{
  _increment(_rejection_threshold);

  _add_hit_type(reason);

  return 1;
}
int
Demerit::extra (int increment, const char * reason)
{
  _increment(increment);

  _add_hit_type(reason);

  return 1;
}

int
Demerit::extra (int increment, const IWString & reason)
{
  _increment(increment);

  _add_hit_type(reason);

  return 1;
}

void
Demerit::_add_hit_type (const char * reason)
{
  if (_types.length())
    _types += ": ";

  _types += reason;

  return;
}

void
Demerit::_add_hit_type (const IWString & reason)
{
  if (_types.length())
    _types += ": ";

  _types += reason;

  return;
}*/
