#include "ttitem.h"

#include <stdlib.h>

#include <random>

void
TTItem::_default_values()
{
  _times_in_training_set = 0;

  _in_training_set = 0;

  std::random_device rd;

  _rng.seed(rd());

  _rnum = _u(_rng);

  return;
}

TTItem::TTItem()
{
  _default_values();

  return;
}

TTItem::TTItem(const const_IWSubstring& s) : _id(s)
{
  _default_values();

  return;
}

void
TTItem::set_in_training_set(int s)
{
  _in_training_set = s;

  if (_in_training_set) {
    _times_in_training_set++;
  }

  return;
}

/*
  When choosing whether or not to volunteer, we need to always try
  to move towards the average. The degree of movement will be a
  function of the "urgency" of getting there.
*/

int
TTItem::volunteer(int splits_remaining, int times_should_have_been_chosen)
{
  int rc;

  if (0 == _times_in_training_set || 0 == times_should_have_been_chosen) {
    rc = 1;
  } else if (_times_in_training_set > times_should_have_been_chosen) {
    rc = _volunteer_lower_probability(splits_remaining, times_should_have_been_chosen);
  } else {
    rc =
        _volunteer_increased_probability(splits_remaining, times_should_have_been_chosen);
  }

  if (rc) {
    set_in_training_set(1);
  }

  return rc;
}

int
TTItem::_volunteer_lower_probability(int splits_remaining,
                                     int times_should_have_been_chosen)
{
  assert(_times_in_training_set > times_should_have_been_chosen);

  double p = static_cast<double>(times_should_have_been_chosen) /
             static_cast<double>(_times_in_training_set);

  if (_u01(_rng) < p) {
    return 1;
  } else {
    return 0;
  }
}

int
TTItem::_volunteer_increased_probability(int splits_remaining,
                                         int times_should_have_been_chosen)
{
  assert(_times_in_training_set <= times_should_have_been_chosen);

  double p = static_cast<double>(_times_in_training_set) /
             static_cast<double>(times_should_have_been_chosen);

  if (_u01(_rng) < p) {
    return 1;
  } else {
    return 0;
  }
}

void
TTItem::switch_training_set_membership()
{
  if (_in_training_set) {
    _in_training_set = 0;
    _times_in_training_set--;
  } else {
    _in_training_set = 1;
    _times_in_training_set++;
  }

  return;
}

int
Times_in_Training_Set_Comparitor::operator()(const TTItem* pt1, const TTItem* pt2) const
{
  int t1 = pt1->times_in_training_set();
  int t2 = pt2->times_in_training_set();

  if (t1 < t2) {
    return -1;
  }
  if (t1 > t2) {
    return 1;
  }
  return 0;
}
