#include <stdlib.h>
#include <assert.h>

#include <iostream>

#include "Foundational/iwstring/iwstring.h"
#include "logical_expression.h"
#include "set_or_unset.h"

using std::cerr;

/*
  Results will be stored as 0 and 1's. We need a special value which
  means missing or unknown value. The expression cannot be evaluated
  if there are missing values
*/

#define IW_LOGEXP_UNKNOWN_RESULT -1

static constexpr char kAnd = '&';
static constexpr char kOr = ',';
static constexpr char kXor = '^';
static constexpr char kLPAnd = ';';

/*
  If we are evaluating something like
    '*,*&*'
  once we evaluate the first condition, the whole expression is true.
  But something like
    '*,*&*;*'
  needs to be evaluated further.
  Some results are not needed, and we need a flag to indicate that
*/

#define IW_LOGEXP_RESULT_NOT_NEEDED -2

IW_Logical_Expression::IW_Logical_Expression()
{
  resize(4);

  _operator.resize(3);
  _unary_operator.resize(3);

  _default_state();

  return;
}

void
IW_Logical_Expression::_default_state()
{
// by default, we have a single (as yet unknown) result

  add(IW_LOGEXP_UNKNOWN_RESULT);

  _unary_operator.add(1);

  _first_low_priority_and_grouping = nullptr;

  return;
}

IW_Logical_Expression::~IW_Logical_Expression()
{
  if (nullptr != _first_low_priority_and_grouping)
    delete _first_low_priority_and_grouping;

  assert(ok());

  return;
}

int
IW_Logical_Expression::ok() const
{
  if (1 == _number_elements && 0 == _operator.number_elements())    // just one component
    return 1;

  if (_number_elements != _unary_operator.number_elements())
    return 0;

  if (_number_elements == _operator.number_elements() + 1)
    return 1;

  if (0 == _number_elements && 0 == _operator.number_elements())
    return 1;

#ifdef VERY_CAREFUL_OK
  if (! resizable_array<int>::ok())
    return 0;

  if (! _operator.ok())
    return 0;

  if (! _unary_operator.ok())
    return 0;
#endif

  return 0;
}

static const char char_ops[] = { '?', '&', '|', '^', ';'};

static int
integer_equivalent(char op)
{
  for (int i = 1; i < 5; i++)
  {
    if (char_ops[i] == op)
      return i;
  }

  if (kOr == op)    // the Daylight smarts OR operator
    return 2;

  cerr << "What kind of operator is this '" << op << "'\n";
  return -1;
}

int 
IW_Logical_Expression::debug_print(std::ostream & os) const
{
  os << "Expression with " << _operator.number_elements() << " operators ";

  if (! ok())
  {
    os << "Yipes, object is invalid\n";
  }

  if (0 == _number_elements)
  {
    os << "Object has no values\n";
    if (_operator.number_elements())
      os << "Yipes, it has operators\n";

    return os.good();
  }

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == _unary_operator[i])
      os << '!';

    if (IW_LOGEXP_UNKNOWN_RESULT == _things[i])
      os << '.';
    else
      os << _things[i];

    if (i < _number_elements - 1)    // the last component does not have an operator
      os << char_ops[_operator[i]];
  }

  os << '\n';

  if (nullptr != _first_low_priority_and_grouping)
    _first_low_priority_and_grouping->debug_print(os);

  return os.good();
}

int
IW_Logical_Expression::RemoveAllOperators()
{
  resize_keep_storage(0);
  _operator.resize_keep_storage(0);
  _unary_operator.resize_keep_storage(0);

  _default_state();

  return 1;
}  

int
IW_Logical_Expression::set_unary_operator(int i, int n)
{
  assert (i >= 0);

  if (i >= _unary_operator.number_elements())
  {
    _unary_operator.extend(i + 1, 1);
    resizable_array<int>::extend(i + 1, IW_LOGEXP_UNKNOWN_RESULT);
  }

  _unary_operator[i] = (n != 0);

  return 1;
}

void
IW_Logical_Expression::reset()
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = IW_LOGEXP_UNKNOWN_RESULT;
  }

  if (nullptr != _first_low_priority_and_grouping)
    _first_low_priority_and_grouping->reset();

  return;
}

int
IW_Logical_Expression::set_result(int i, int r)
{
  assert (ok());

  assert (ok_index(i));

// Be careful not to propagate values extended by an OR operation

  if (IW_LOGEXP_UNKNOWN_RESULT != _things[i])
    return 1;

  if (r)
    _things[i] = 1;
  else
    _things[i] = 0;

  if (0 == _unary_operator[i])
    _things[i] = ! _things[i];

  assert (ok());

  return 1;
}

/*
  Sometimes as a result of OR operations, we may no longer need a given value
*/

int
IW_Logical_Expression::result_needed(int i) const
{
  if (ok_index(i))
    return IW_LOGEXP_UNKNOWN_RESULT == _things[i];

  if (0 == _number_elements && 0 == _operator.number_elements())
    return 1;

  cerr << "IW_Logical_Expression::result_needed: result " << i << " is bad\n";
  debug_print(cerr);
  abort();
  return 0;
}

int
IW_Logical_Expression::add_operator(int op)
{
  assert (ok());

  if (IW_LOGEXP_AND == op)
    ;
  else if (IW_LOGEXP_OR == op)
    ;
  else if (IW_LOGEXP_XOR == op)
    ;
  else if (IW_LOGEXP_LOW_PRIORITY_AND == op)
    ;
  else
  {
    cerr << "IW_Logical_Expression::add_operator: invalid operator " << op << '\n';
    abort();
    return 0;
  }

  if (IW_LOGEXP_XOR == op)
  {
    if (_operator.number_elements())
    {
      cerr << "IW_Logical_Expression::add_operator: XOR can only be the first operator\n";
      return 0;
    }
  }
  else if (_operator.empty())    // no need to check
    ;
  else if (IW_LOGEXP_XOR == _operator[0])
  {
    cerr << "IW_Logical_Expression::add_operator: XOR can only be in a simple expression\n";
    return 0;
  }

  _operator.add(op);

  resizable_array<int>::add(IW_LOGEXP_UNKNOWN_RESULT);
  _unary_operator.add(1);

  _first_low_priority_and_grouping = nullptr;

  assert (ok());

  return 1;
}

int
IW_Logical_Expression::add_operator(char op)
{
  int i = integer_equivalent(op);
  
  assert(i > 0);

  return add_operator(i);
}

int
IW_Logical_Expression::_set_all_operators(int op)
{
  assert (ok());

  _operator.set_all(op);

  return 1;
}

int
IW_Logical_Expression::set_all_operators(int op)
{
  if (IW_LOGEXP_AND == op)
    return _set_all_operators(op);
  
  if (IW_LOGEXP_OR == op)
    return _set_all_operators(op);
  
  if (IW_LOGEXP_XOR == op)
  {
    if (1 != _operator.number_elements())
    {
      cerr << "IW_Logical_Expression::set_all_operators: XOR can only be in simple expressions\n";
      return 0;
    }
    return _set_all_operators(op);
  }
  
  if (IW_LOGEXP_LOW_PRIORITY_AND == op)
    return _set_all_operators(op);

  cerr << "IW_Logical_Expression::_set_all_operators: what operator is this " << op << '\n';
  return 0;
  
}

int
IW_Logical_Expression::set_all_operators(char op)
{
  int i = integer_equivalent(op);

  assert (i > 0);

  return set_all_operators(i);
}

int
IW_Logical_Expression::_all_operators_are(int op) const
{
  for (int i = 0; i < _operator.number_elements(); i++)
  {
    if (op != _operator[i])
      return 0;
  }

  return 1;
}

int
IW_Logical_Expression::all_operators_are(char op) const
{
  int i = integer_equivalent(op);

  assert (i > 0);

  return _all_operators_are(i);
}

int
IW_Logical_Expression::all_operators_are(int op) const
{
  if (IW_LOGEXP_AND == op)
    return _all_operators_are(op);
  if (IW_LOGEXP_OR == op)
    return _all_operators_are(op);
  if (IW_LOGEXP_XOR == op)
    return _all_operators_are(op);
  if (IW_LOGEXP_LOW_PRIORITY_AND == op)
    return _all_operators_are(op);

  cerr << "IW_Logical_Expression::all_operators_are: what operator is this " << op << '\n';
  return 0;
}

int
IW_Logical_Expression::op(int i, const_IWSubstring & op) const
{
  assert (_operator.ok_index(i));

  op = char_ops[_operator[i]];

  return 1;
}

int
IW_Logical_Expression::op(int i) const
{
  assert (_operator.ok_index(i));

  return _operator[i];
}

int
IW_Logical_Expression::evaluate(int & zresult)
{
  assert (ok());

  if (1 == _number_elements)     // the most common case
  {
    if (IW_LOGEXP_UNKNOWN_RESULT == _things[0])
      return 0;

//  cerr << "What about " << _unary_operator[0] << " and " << _things[0] << '\n';
    zresult = _things[0];
    return 1;
  }

  if (0 == _number_elements)
  {
    cerr << "IW_Logical_Expression::evaluate: cannot evaluate empty expression\n";
    return 0;
  }

  if (IW_LOGEXP_UNKNOWN_RESULT == _things[0])    // no values known
    return 0;

// The case of just one operator is also handled as a special case
 
  if (2 == _number_elements)
    return _evaluate_single_operator(zresult);

  if (nullptr == _first_low_priority_and_grouping)
  {
    if (! _initialise())
    {
      cerr << "IW_Logical_Expression::evaluate: cannot initialise expression structures\n";
      return 0;
    }
  }

  IW_Logexp_Low_Priority_and_Grouping * current = _first_low_priority_and_grouping;

  do 
  {
    if (! current->evaluate(_things, zresult))
      return 0;

//  If any of our low priority and groupings are false, we are done

    if (0 == zresult)
      return 1;

    current = current->next();

  } while (nullptr != current);

// Each low priority and grouping was true

  zresult = 1;

  return 1;
}

int
IW_Logical_Expression::_initialise()
{
  assert (nullptr == _first_low_priority_and_grouping);

// Maybe we should be more generous with XOR - perhaps enforce this in the low priority and grouping...

  for (int i = 0; i < _operator.number_elements(); i++)
  {
    if (IW_LOGEXP_XOR == _operator[i])
    {
      cerr << "IW_Logical_Expression::_initialise: XOR operator not allowed in complex expressions\n";
      return 0;
    }
  }

  _first_low_priority_and_grouping = new IW_Logexp_Low_Priority_and_Grouping;

  int istart = 0;

  if (! _first_low_priority_and_grouping->initialise(_operator, istart))
  {
    cerr << "IW_Logical_Expression::_initialise: cannot initialise AND grouping\n";
    return 0;
  }

//#define DEBUG_LGEXP_INITIALISE
#ifdef DEBUG_LGEXP_INITIALISE
  cerr << "After initialising\n";
  debug_print(cerr);
#endif

  return 1;
}

int
IW_Logical_Expression::_evaluate_single_operator(int & zresult)
{
  if (IW_LOGEXP_OR == _operator[0])
  {
    if (IW_LOGEXP_UNKNOWN_RESULT == _things[1])     // only the LHS is known
    {
      if (0 == _things[0])     // maybe the rhs will turn out to be true
        return 0;

//    LHS is true. We don't need the rhs

      zresult = 1;
      _things[1] = IW_LOGEXP_RESULT_NOT_NEEDED;

      return 1;
    }

//  both lhs and rhs are known

    zresult = (_things[0] || _things[1]);

    return 1;
  }

// All other operators need the rhs

  if (IW_LOGEXP_UNKNOWN_RESULT == _things[1])
    return 0;

  if (IW_LOGEXP_AND == _operator[0] || IW_LOGEXP_LOW_PRIORITY_AND == _operator[0])
  {
    zresult = (_things[0] && _things[1]);
  }
  else if (IW_LOGEXP_XOR == _operator[0])
  {
    if (_things[0] && _things[1])
      zresult = 0;
    else if (0 == _things[0] && 0 == _things[1])
      zresult = 0;
    else
      zresult = 1;
  }
  else
  {
    cerr << "IW_Logical_Expression::_evaluate_single_operator: what kind of operator is this " << _operator[0] << '\n';
    abort();
  }

  return 1;
}

IW_Logexp_AND_Grouping::IW_Logexp_AND_Grouping()
{
  _n = -1;
  _istart = -1;

  _next = nullptr;

  return;
}

IW_Logexp_AND_Grouping::~IW_Logexp_AND_Grouping ()
{
  if (nullptr != _next)
    delete _next;

  return;
}

int
IW_Logexp_AND_Grouping::debug_print(std::ostream & os) const
{
  os << "AND grouping starting at " << _istart << " with " << _n << " items, result ";

  int zresult;
  if (_result.value(zresult))
    os << zresult << '\n';
  else
    os << "unknown\n";

  if (nullptr == _next)
    return os.good();

  return _next->debug_print(os);
}

void
IW_Logexp_AND_Grouping::reset()
{
  _result.unset();

  if (_next)
    _next->reset();

  return;
}

//#define DEBUG_AND_GROUP_INITIALISE

/*
  Group together one or more results that are part of a high priority
  and grouping. Must be very careful because a single result is considered
  part of an and grouping,

*/

int
IW_Logexp_AND_Grouping::initialise(const resizable_array<int> & the_operators,
                                    int & istart)
{
  assert (_n < 0);    // must not have been initialised before

#ifdef DEBUG_AND_GROUP_INITIALISE
  cerr << "Initialising high priority AND on " << the_operators.number_elements() << " operators, istart = " << istart << '\n';
#endif

  _istart = istart;
  _n = 1;

  int nop = the_operators.number_elements();

  while (istart < nop)
  {
    if (IW_LOGEXP_AND != the_operators[istart])
      break;

    istart++;
    _n++;
  }
  
  istart++;

#ifdef DEBUG_AND_GROUP_INITIALISE
  cerr << "consumed " << _n << " items\n";
#endif

  return 1;
}

//#define DEBUG_AND_GROUP_EVALUATE

int
IW_Logexp_AND_Grouping::evaluate(int * zresults,
                                 int & zresult)
{
#ifdef DEBUG_AND_GROUP_EVALUATE
  cerr << "Evaluating and group starting at " << _istart;
  if (_result.value(zresult))
    cerr << ". Result " << zresult << '\n';
  else
    cerr << ". Result unknown\n";
#endif

  if (_result.value(zresult))
    return 1;
    
  assert (_n > 0 && _istart >= 0);

  for (int i = _istart; i < _istart + _n; i++)
  {
#ifdef DEBUG_AND_GROUP_EVALUATE
    cerr << "Result " << i << " is " << zresults[i] << '\n';
#endif

    if (0 == zresults[i])
    {
      _mark_remaining_results_not_needed(zresults, i + 1);
      _result.set(0);
      zresult = 0;
      return 1;
    }

    if (IW_LOGEXP_UNKNOWN_RESULT == zresults[i])
      return 0;
  }

// All results known and true

  zresult = 1;
  _result.set(1);

  return 1;
}

/*
  We have found an answer to be false. All subsequent answers don't matter
*/

void
IW_Logexp_AND_Grouping::_mark_remaining_results_not_needed(int * zresults, int istart) const
{
  int istop = _istart + _n;
//cerr << "Marking results from " << istart << " to " << istop << " _n = " << _n << '\n';

  while (istart < istop)
  {
    zresults[istart] = IW_LOGEXP_RESULT_NOT_NEEDED;
    istart++;
  }

  return;
}

void
IW_Logexp_AND_Grouping::set_not_needed(int * zresults) const
{
  for (int i = _istart; i < _istart + _n; i++)
  {
    zresults[_istart + i] = IW_LOGEXP_RESULT_NOT_NEEDED;
  }

  return;
}


IW_Logexp_Low_Priority_and_Grouping::IW_Logexp_Low_Priority_and_Grouping()
{
  _n = 0;

  _istart = -1;

  _first_and_chunk = nullptr;

  _next = nullptr;

  return;
}

IW_Logexp_Low_Priority_and_Grouping::~IW_Logexp_Low_Priority_and_Grouping()
{
  if (nullptr != _first_and_chunk)
    delete _first_and_chunk;

  if (nullptr != _next)
    delete _next;

  return;
}

int
IW_Logexp_Low_Priority_and_Grouping::debug_print(std::ostream & os) const
{
  os << "Low priority and grouping starting at " << _istart << " covering " << _n << " results. Result ";

  int zresult;
  if (_result.value(zresult))
    os << zresult << '\n';
  else
    os << " unknown\n";

  if (_first_and_chunk)
    _first_and_chunk->debug_print(os);

  if (nullptr == _next)
    return os.good();

  return _next->debug_print(os);


}

void
IW_Logexp_Low_Priority_and_Grouping::reset()
{
  _result.unset();

  if (_first_and_chunk)
    _first_and_chunk->reset();

  if (_next)
    _next->reset();

  return;
}

//#define DEBUG_LOW_PRIORITY_AND_G_INITIALISE

/*
  Pretty confusing stuff here.

  istart is an index into the array of operators. If there are N operators,
  there will be N+1 results
*/

int
IW_Logexp_Low_Priority_and_Grouping::initialise(const resizable_array<int> & the_operators,
                                        int & istart)
{
  assert (0 == _n);
  assert (nullptr == _first_and_chunk);
  assert (nullptr == _next);

#if defined(DEBUG_LOW_PRIORITY_AND_G_INITIALISE)
  cerr << "Initialising low priority and grouping with " << the_operators.number_elements() << " operators, istart = " << istart << '\n';
#endif

  _istart = istart;

  _first_and_chunk = new IW_Logexp_AND_Grouping();

  IW_Logexp_AND_Grouping * current = _first_and_chunk;

  int nop = the_operators.number_elements();

  while (1)
  {
    if (! current->initialise(the_operators, istart))
    {
      cerr << "IW_Logexp_Low_Priority_and_Grouping::initialise: failed to initialise chunk\n";
      return 0;
    }

#ifdef DEBUG_LOW_PRIORITY_AND_G_INITIALISE
    cerr << "And group consumed " << current->number_results() << " results, istart = " << istart << '\n';
#endif

    if (istart > nop)
      break;
    
    if (IW_LOGEXP_LOW_PRIORITY_AND == the_operators[istart - 1])   // need to move on to the next grouping
      break;

    IW_Logexp_AND_Grouping * tmp = new IW_Logexp_AND_Grouping;

    current->set_next(tmp);

    current = tmp;
  }

  _n = istart - _istart;

#if defined(DEBUG_LOW_PRIORITY_AND_G_INITIALISE)
  cerr << "At end, istart = " << istart << ", nop = " << nop << ", n = " << _n << '\n';
#endif

  if (istart > nop)
    return 1;

  _next = new IW_Logexp_Low_Priority_and_Grouping();

  return _next->initialise(the_operators, istart);
}

//#define  DEBUG_LOW_P_AND_EVALUATE

int
IW_Logexp_Low_Priority_and_Grouping::evaluate(int * the_results,
                                              int & zresult)
{
#ifdef DEBUG_LOW_P_AND_EVALUATE
  cerr << "Checking lpand grouping starting at " << _istart;
  if (_result.value(zresult))    // already determined
    cerr << " already determined " << zresult << '\n';
  else
    cerr << " unknown\n";
#endif

  if (_result.value(zresult))    // already determined
    return 1;

  IW_Logexp_AND_Grouping * current = _first_and_chunk;
  assert (nullptr != current);

// Since we just have OR operators, we just look for the first grouping
// that is true

  while (nullptr != current)
  {
    int tmpresult;

#ifdef DEBUG_LOW_P_AND_EVALUATE
    cerr << "Trying to evaluate AND group at " << current->istart() << '\n';
#endif

    if (! current->evaluate(the_results, tmpresult))    // cannot be evaluated yet
      return 0;

    if (tmpresult)    // could be evaluated and is true
    {
      _result.set(1);    // we only need one component to be true
      zresult = 1;
      
      return 1;
    }

    current = current->next();
  }

// none of our components are true, we are done

  _result.set(0);
  zresult = 0;

  return 1;
}

IWString
IW_Logical_Expression::ToString() const {
  IWString result;
  if (_operator.empty()) {
    return result;
  }
  
  for (int i =  0; i < _unary_operator.number_elements(); ++i) {
    if (i > 0) {
      switch (_operator[i - 1]) {
        case IW_LOGEXP_AND:
          result << kAnd;
          break;
        case IW_LOGEXP_OR:
          result << kOr;
          break;
        case IW_LOGEXP_XOR:
          result << kXor;
          break;
        case IW_LOGEXP_LOW_PRIORITY_AND:
          result << kLPAnd;
          break;
        default:
          cerr << "IW_LOGEXP_LOW_PRIORITY_AND::ToString:unrecognised operator " << _operator[i] << '\n';
          result << '?';
      }
    }
    if (_unary_operator[i] == 0) {
      result << '!';
    }
    result << '.';
  }

  return result;
}

// parse `s` into a sequence of tokens.
// Each token in `s` must be either
//     a dot '.' or a negated dot '!.'
//     an operator
static constexpr int kDot = 10;
static constexpr int kNegatedDot = 11;
int
DotsAndOperators(const IWString& s, resizable_array<int>& result) {

  result.resize_keep_storage(0);

  for (int i = 0; i < s.length(); ++i) {
    const char c = s[i];
    char next_c;
    if (i == s.length() - 1) {
      next_c = ' ';
    } else {
      next_c = s[i + 1];
    }
    if (c == '!'  && next_c == '.') {
      result << kNegatedDot;
      ++i;
    } else if (c == '.') {
      result << kDot;
    } else if (c == kAnd) {
      result << IW_LOGEXP_AND;
    } else if (c == kOr) {
      result << IW_LOGEXP_OR;
    } else if (c == kXor) {
      result << IW_LOGEXP_XOR;
    } else if (c == kLPAnd) {
      result << IW_LOGEXP_LOW_PRIORITY_AND;
    } else {
      cerr << "DotsAndOperators:unrecognised character '" << c << " pos " << i << " in '" << s << "'\n";
      return 0;
    }
  }

  return result.number_elements();
}

static int
ContainsAlternatingDotsAndOperators(const resizable_array<int>& tokens) {
  // Must be an odd number of tokens.
  if (tokens.number_elements() % 2 != 1) {
    cerr << "ContainsAlternatingDotsAndOperators:incorrect token count " << tokens.size() << '\n';
    return 0;
  }

  for (int i = 0; i < tokens.number_elements(); ++i) {
    int ti = tokens[i];
    if (i % 2 == 0) {
      if (ti == kDot || ti == kNegatedDot) {
      } else {
        cerr << "Character " << i << " not dot, got " << ti << '\n';
        return 0;
      }
    } else {
      if (ti == IW_LOGEXP_AND || ti == IW_LOGEXP_OR || ti == IW_LOGEXP_XOR || ti == IW_LOGEXP_LOW_PRIORITY_AND) {
      } else {
        cerr << "Character " << i << " not operator, got " << ti << '\n';
        return 0;
      }
    }
  }

  return 1;
}

int
IW_Logical_Expression::BuildFromString(const IWString& s) {
  cerr << "Building from '" << s << "'\n";
  this->resize_keep_storage(0);
  _operator.resize_keep_storage(0);
  _unary_operator.resize_keep_storage(0);

  resizable_array<int> tokens;
  if (! DotsAndOperators(s, tokens)) {
    cerr << "IW_Logical_Expression::BuildFromString:cannot tokenise '" << s << "'\n";
    return 0;
  }

  if (! ContainsAlternatingDotsAndOperators(tokens)) {
    cerr << "IW_Logical_Expression::BuildFromString:invalid syntax '" << s << "'\n";
    return 0;
  }

  for (int i = 0; i < tokens.number_elements(); ++i) {
    if (tokens[i] == kDot) {
      _unary_operator << 1;
      this->add(IW_LOGEXP_UNKNOWN_RESULT);
    } else if (tokens[i] == kNegatedDot) {
      _unary_operator << 0;
      this->add(IW_LOGEXP_UNKNOWN_RESULT);
    }
    else {
      _operator << tokens[i];
    }
  }

  if (! ok()) {
    cerr << "IW_Logical_Expression::BuildFromString:invalid state '" << s << "'\n";
    debug_print(cerr);
    return 0;
  }

  _first_low_priority_and_grouping = nullptr;

  return 1;
}

bool
IW_Logical_Expression::operator==(const IW_Logical_Expression& rhs) const {
  if (_operator.size() != rhs._operator.size()) {
    return false;
  }
  const int nop = _operator.number_elements();
  for (int i = 0; i < nop; ++i) {
    if (_operator[i] != rhs._operator[i]) {
      return false;
    }
  }

  const int nunary = _unary_operator.number_elements();
  for (int i = 0; i < nunary; ++i) {
    if (_unary_operator[i] != rhs._unary_operator[i]) {
      return false;
    }
  }

  return true;
}
