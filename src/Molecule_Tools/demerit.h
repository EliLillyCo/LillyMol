#ifndef IW_DEMERIT_H
#define IW_DEMERIT_H

#include <iostream>
#include "Foundational/iwstring/iwstring.h"

#define DEFAULT_REJECTION_THRESHOLD 100

class Demerit_and_Reason
{
  private:
    const int _demerit;
    const IWString _reason;

  public:
    Demerit_and_Reason(const int d, const const_IWSubstring & r) : _demerit(d), _reason(r) {};

    int demerit() const { return _demerit;}
    const IWString & reason() const { return _reason;}
};

/*
  The variable _rejected keep track of whether or not a molecule
  has been specifically rejected by any one rule
*/

class Demerit
{
  private:
    resizable_array_p<Demerit_and_Reason> _demerit;
    int _score;

    int _append_demerit_value_with_reason;

//  We can optionally hold demerit values associated with each atom.
//  can be something we allocate, or something given to us

    int * _atomic_demerit;
    int _need_to_delete_atomic_demerit;    // will be true if we allocated it

//  private functions

    int _add_demerit(const int, const const_IWSubstring &);
    void _increment(int, const const_IWSubstring &);
//  void _add_hit_type(const char *);
//  void _add_hit_type(const IWString &);
    void _add_hit_type(int, const const_IWSubstring &);

  public:
    Demerit();
    Demerit(Demerit&& rhs);
    ~Demerit();

    int ok() const;
    int debug_print(std::ostream &) const;

    void set_append_demerit_value_with_reason(const int s) { _append_demerit_value_with_reason = s;}

    int initialise_atomic_demerits(const int matoms);
    int set_atomic_demerit_array(int *);

    int * atomic_demerit() const { return _atomic_demerit;}

    int score() const { return _score;};
    int number_different_demerits_applied() const { return _demerit.number_elements();}

//  int max_demerit_encountered() const { return _max_demerit_encountered;}
//  const IWString & max_demerit_reason() const { return _max_demerit_reason;}

//  const IWString & types() const { return _types;}

//  int extra(int, const char *);
//  int extra(int, const IWString &);
    int extra(int, const const_IWSubstring);
    int reject(const const_IWSubstring);
//  int reject(const char *);
//  int reject(const IWString &);
    int rejected() const;
    int rejected_by_single_rule() const { return 1 == _demerit.number_elements();}

    void do_sort();

    int write_in_tdt_form(std::ostream &) const;

    const resizable_array_p<Demerit_and_Reason> & demerits() const { return _demerit;}
};

extern void set_rejection_threshold(int);
extern int  rejection_threshold();

/*
  Sept 2004. Dan Robertson needs to get all reasons separated with
  their individual demerit values. If that's the case, we build the
  _TYPE variable differently
*/

extern void set_demerit_reason_contains_individual_demerits(int);
extern void set_store_demerit_reasons_like_tsubstructure(int);

#endif
