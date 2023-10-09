#ifndef UTILITIES_GENERAL_TTITEM_H
#define UTILITIES_GENERAL_TTITEM_H

#include <random>

#include "Foundational/iwstring/iwstring.h"

class TTItem
{
  private:
    IWString _id;

    int _times_in_training_set;

    int _in_training_set;

    uint32_t _rnum;   // for getting random distributions of things

    std::mt19937 _rng;
    std::uniform_int_distribution<uint32_t> _u;
    std::uniform_real_distribution<float> _u01;

//  private functions

    void _default_values();

    int _volunteer_lower_probability (int splits_remaining, int times_should_have_been_chosen);
    int _volunteer_increased_probability (int splits_remaining, int times_should_have_been_chosen);

  public:
    TTItem();
    TTItem(const const_IWSubstring &);

    const IWString & id() const { return _id;}
    void set_id (const const_IWSubstring & s) { _id = s;}

    uint32_t rnum () const { return _rnum;}

    int times_in_training_set() const { return _times_in_training_set;}

    int in_training_set() const { return _in_training_set;}
    void set_in_training_set (int s);

    int volunteer (int, int);

    void switch_training_set_membership();

    void reset_random_number () {
      _rnum = _u(_rng);
    }
};

class Times_in_Training_Set_Comparitor
{
  private:
  public:
    int operator() (const TTItem *, const TTItem *) const;
};

#endif  // UTILITIES_GENERAL_TTITEM_H
