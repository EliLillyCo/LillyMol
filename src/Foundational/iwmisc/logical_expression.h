#ifndef IW_LOGICAL_EXP_H
#define IW_LOGICAL_EXP_H

class const_IWSubstring;

/*
  In substructure searching we need to be able to evaluate something
  like

  a|b&c

  We store the results (a, b, c) in the base resizable_array<int> and
  the operators in the _operator array

  Calls to the evaluate function return 0 until there are enough
  results for it to produce a result. The result (0 or 1) will go
  in its argument.

  For example, consider
  'a|b&c'

  If only result[0] is known (false), the expression cannot be evaluated.
  Only when all values have been specified can a evaluate succeed.
  
  However, something like  'a&b;c'
  can be evaluated once a and b are known if either a or b is false.
*/

// These numbers must be kept consistent with the char_ops array in the source file

#define IW_LOGEXP_UNDEFINED 0
#define IW_LOGEXP_AND 1
#define IW_LOGEXP_OR 2
#define IW_LOGEXP_XOR 3
#define IW_LOGEXP_LOW_PRIORITY_AND 4

#include <iostream>
using std::cerr;
using std::endl;

#include "iwaray.h"
#include "set_or_unset.h"

/*
  Since the High Priority AND is the tighest binding operator, all
  items in a logical expression are first reduced to sets of AND groups.
  These groups may have just one member
*/

class IW_Logexp_AND_Grouping
{
  private:
    int _n;       // the number of results in our grouping

    int _istart;   // index of first result in the array or results

    Set_or_Unset<int> _result;

    IW_Logexp_AND_Grouping * _next;

//  private functions

    void _mark_remaining_results_not_needed (int * zresults, int istart) const;

  public:
    IW_Logexp_AND_Grouping ();
    ~IW_Logexp_AND_Grouping ();

    int debug_print (std::ostream &) const;

    int istart () const { return _istart;}

    IW_Logexp_AND_Grouping * next () const { return _next;};
    void set_next (IW_Logexp_AND_Grouping * n) { _next = n;}

    int number_results () const { return _n;}

    int initialise (const resizable_array<int> &, int &);

    void reset ();

    void set_not_needed (int * zresults) const;

    int evaluate (int *, int &);
    int evaluate (int *, int &, int);
};

/*
  The low priority and is the lowest priority grouping. We break the
  expression up onto a linked list of these
*/

class IW_Logexp_Low_Priority_and_Grouping
{
  private:

    Set_or_Unset<int> _result;

    int _n;      // the number of values

    int _istart;    // index of where we start in the array of results

    IW_Logexp_AND_Grouping * _first_and_chunk;

    IW_Logexp_Low_Priority_and_Grouping * _next;

//  private functions

    int initialise (const int *, const int *, int, int, int);

  public:
    IW_Logexp_Low_Priority_and_Grouping ();
    ~IW_Logexp_Low_Priority_and_Grouping ();

    int debug_print (std::ostream &) const;

    IW_Logexp_Low_Priority_and_Grouping * next () const { return _next;}

    void reset ();

    int initialise (const resizable_array<int> &, int &);

    int evaluate (int * the_results, int & zresult);
};

class IW_Logical_Expression : private resizable_array<int>
{
  private:
    resizable_array<int> _operator;

//  We may want to evaluate 'A & (!B)'

    resizable_array<int> _unary_operator;

    IW_Logexp_Low_Priority_and_Grouping * _first_low_priority_and_grouping;

//  private functions

    int _evaluate_single_operator (int & zresult);

    int _set_all_operators (int);
    int _all_operators_are (int) const;

    int _initialise ();

  public:
    IW_Logical_Expression ();
    ~IW_Logical_Expression ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int number_results () const { return _number_elements;}

    int evaluate (int &);

    void reset ();

    int add_operator (int);
    int add_operator (char);
    int number_operators () const { return _operator.number_elements ();}
    int set_all_operators (char);
    int set_all_operators (int);
    int all_operators_are (char) const;
    int all_operators_are (int) const;

    int set_unary_operator (int, int);
    int unary_operator (int i) const { return _unary_operator[i];}

    int op (int) const;
    int op (int, const_IWSubstring &) const;

    int result_needed (int) const;
    int set_result (int, int);
    int result (int i) const { return _things[i];}
};

#endif
