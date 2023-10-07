#ifndef FOUNDATIONAL_IWMISC_MATCHER_H
#define FOUNDATIONAL_IWMISC_MATCHER_H

#include <algorithm>
#include <iostream>

#include "Foundational/iwstring/iwstring.h"

// A class that is at the heart of substructure matching.

namespace iwmatcher {

// The indices into the _mask array.
constexpr int kSet = 0;
constexpr int kMin = 1;
constexpr int kMax = 2;

template <typename T>
class Matcher {
  private:
    // Set to 1 if 
    // 0 _single_value is valid
    // 1 _min_val is valid;
    // 2 _min_val is valid;
    // 3 not used.
    // By making a union, returning _set will be non zero if any of the
    // components of _mask are set.
    union MaskAndInt {
      unsigned char _mask[4];
      unsigned int _set;
    };

    MaskAndInt _mask_and_int;

    // If we are matching just a single value, it is held here.
    T _single_value;

    // If there is a min or max, those are held here.
    T _min_val;
    T _max_val;

    // If more than one value is being matched, hold those
    // here. Note that the first value will be in _single_value
    // so this only gets populated if we have > 1 feature to match.
    T* _other_values;
    int _number_other_values;

    // Private functions
    void _reset_other_values();

  public:
    Matcher();
    ~Matcher();

    int ok() const;
    int debug_print(std::ostream &) const;

    // For compatibility with Min_Max_Specifier which was public resizable_array
    int number_elements() const;
    unsigned int size() const { return number_elements();}
    // This is only defined for r = 0;
    int resize(int r);
    // return the smallest explicit value if present.
    int min_explicit_value(T& v) const;
    // Enable indexing
    T operator[](int ndx) const;

    int is_set() const { return _mask_and_int._set;}
    int match_any() const { return ! is_set();}
    void set_match_any(int s);

    void reset();

    int initialise(const const_IWSubstring & buffer);

    int add(T);    // we override this so we can maintain _is_set
    int add_if_not_already_present(T);  // must override.

    // Note that these do no checking. Call ok() afterwards if needed.
    int set_min(T v);
    int set_max(T v);

    int adjust_to_accommodate(const T v);

    int min(T & v) const;
    int max(T & v) const;

    int matches(T) const;

    int write_msi(std::ostream &os, const char * name, int indentation = 0) const;

    int write_compact_representation(std::ostream &) const;  // something like '3,4', '>1', '<6'

    Matcher<T> & operator=(const Matcher<T> &);
    bool operator==(const Matcher<T>& rhs) const;

    // Return a vector containing any numeric values we match.
    resizable_array<T> ValuesMatched() const;

    // Returns true if any of the explicit values match `condition`.
    template <typename C> bool AnyValue(C condition) const;
    // For each explicit value, update the value with `fn`.
    template <typename C> void UpdateValues(C fn);
};

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(MATCHER_IMPLEMENTATION)

template <typename T>
Matcher<T>::Matcher() {
  _mask_and_int._set = 0;

  _other_values = nullptr;
  _number_other_values = 0;
}

template <typename T>
Matcher<T>::~Matcher() {
  if (_other_values != nullptr) {
    delete [] _other_values;
  }
}

template <typename T>
int
Matcher<T>::min(T& v) const {
  if (_mask_and_int._mask[kMin] == 0) {
    return 0;
  }
  v = _min_val;
  return 1;
}

template <typename T>
int
Matcher<T>::max(T& v) const {
  if (_mask_and_int._mask[kMax] == 0) {
    return 0;
  }
  v = _max_val;
  return 1;
}

template <typename T>
int
Matcher<T>::matches(T v) const {
  // Check the most common case first.
  if (_mask_and_int._mask[kSet]) {
    if (_single_value == v) {
      return 1;
    }
  }

  int passed_conditions = 0;
  if (_mask_and_int._mask[kMin] != 0) {
    if (v < _min_val) {
      return 0;
    }
    ++passed_conditions;
  }
  if (_mask_and_int._mask[kMax] != 0) {
    if (v > _max_val) {
      return 0;
    }
    ++passed_conditions;
  }

  if (_number_other_values == 0) {
    if (_mask_and_int._set == 0) {  // Nothing specified, match anything.
      return 1;
    }
    return passed_conditions;
  }

  // Having two values should be common, so a quick check before starting a loop.
  if (v == _other_values[0]) {
    return 1;
  }

  for (int i = 1; i < _number_other_values; ++i) {
    if (v == _other_values[i]) {
      return 1;
    }
  }

  return 0;
}

template <typename T>
int
Matcher<T>::debug_print(std::ostream& output) const {
  output << "Matcher";
  if (_mask_and_int._mask[kMin]) {
    output << " min " << _min_val;
  }
  if (_mask_and_int._mask[kMax]) {
    output << " max " << _max_val;
  }
  if (_mask_and_int._mask[kSet]) {
    output << " value " << _single_value;
  }
  if (_number_other_values > 0) {
    output << " values";
    for (int i = 0; i < _number_other_values; ++i) {
      output << ' ' << _other_values[i];
    }
  }

  output << '\n';

  return output.good();
}

template <class T>
std::ostream &
operator << (std::ostream & os, const Matcher<T> & matcher) 
{
  os << "Matcher: ";

  T tmp;
  if (matcher.min(tmp))
    os << "min " << tmp;
  if (matcher.max(tmp))
    os << " max " << tmp;
  
  if (matcher.match_any()) {
    os << " matches any value";
    return os;
  }

  const resizable_array<T> matched = matcher.ValuesMatched();
  os << " matches these values";
  for (int m : matched) {
    os << ' ' << m;
  }

  return os;
}

template <typename T>
void
Matcher<T>::_reset_other_values() {
  delete [] _other_values;
  _other_values = nullptr;
  _number_other_values = 0;
}

template <typename T>
int
Matcher<T>::ok() const {
  if (_mask_and_int._mask[kMin] && _mask_and_int._mask[kMax]) {
    if (_min_val > _max_val) {
      //std::cerr << "maxmin invert\n";
      return 0;
    }
  }

  // If there are multiple values, the first must be in _single_value.
  if (_mask_and_int._mask[kSet] == 0 && _other_values != nullptr) {
    //std::cerr << "single value not set, but others are\n";
    return 0;
  }

  if (_mask_and_int._mask[kSet]) {
    if (_mask_and_int._mask[kMin] && _single_value < _min_val) {
      //std::cerr << "Min set but single value inconsistent\n";
      return 0;
    }
    if (_mask_and_int._mask[kMax]  && _single_value > _max_val) {
      //std::cerr << "Max set but single value inconsistent\n";
      return 0;
    }
  }

  for (int i = 0; i < _number_other_values; ++i) {
    if (_mask_and_int._mask[kMin] && _other_values[i] < _min_val) {
      return 0;
    }
    if (_mask_and_int._mask[kMax] && _other_values[i] > _max_val) {
      return 0;
    }
  }

  return 1;
}

template <typename T>
void
Matcher<T>::reset() {
  _mask_and_int._set = 0;

  if (_other_values != nullptr) {
    _reset_other_values();
  }
}

template <typename T>
void
Matcher<T>::set_match_any(int s) {
  if (s == 1) {
    _mask_and_int._set = 0;
    _reset_other_values();
    return;
  }
  // Not sure what to do if s is 0??
  std::cerr << "Matcher::set_match_any:what to do with 0\n";
  return;
}

template <typename T>
int
Matcher<T>::set_min(T v) {
  _min_val = v;
  _mask_and_int._mask[kMin] = 1;
  return 1;
}

template <typename T>
int
Matcher<T>::set_max(T v) {
  _max_val = v;
  _mask_and_int._mask[kMax] = 1;
  return 1;
}

template <typename T>
int
Matcher<T>::add(T v) {
  if (_mask_and_int._mask[kSet] == 0) {
    _single_value = v;
    _mask_and_int._mask[kSet] = 1;
    return 1;
  }

  if (_other_values == nullptr) {
    _other_values = new T[1];
    _number_other_values = 1;
    _other_values[0] = v;
    return 1;
  }

  T * new_other = new T[_number_other_values + 1];
  std::copy_n(_other_values, _number_other_values, new_other);
  new_other[_number_other_values] = v;
  ++_number_other_values;
  delete [] _other_values;
  _other_values = new_other;
  return 1;
}

template <typename T>
int
Matcher<T>::add_if_not_already_present(T value) {
  if (_mask_and_int._mask[kSet] == 0) {
    return add(value);
  }
  if (value == _single_value) {
    return 0;
  }

  for (int i = 0; i < _number_other_values; ++i) {
    if (value == _other_values[i]) {
      return 0;
    }
  }

  return add(value);
}

/*
  THis is copied from minmaxspc.h

  An initialising string can look like

  1      a digit
  2,3,4  a sequence of digits
  >1     a range
  <2     a range
*/

template <typename T>
int
Matcher<T>::initialise(const const_IWSubstring & buffer)
{
  reset();

  if (buffer.starts_with("<="))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(2);

    T v;
    if (! mybuffer.numeric_value(v))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_max(v);

    return 1;
  }

  if (buffer.starts_with('<'))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(1);

    int m;
    if (! mybuffer.numeric_value(m))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_max(m - 1);    // works for ints, but not floats

    return 1;
  }

  if (buffer.starts_with(">="))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(2);

    T v;
    if (! mybuffer.numeric_value(v))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_min(v);

    return 1;
  }

  if (buffer.starts_with('>'))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(1);

    int m;
    if (! mybuffer.numeric_value(m))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_min(m + 1);      // works with integers but not floats

    return 1;
  }

  int i = 0;
  const_IWSubstring token;

  char separator;
  if (buffer.contains(','))
    separator = ',';
  else
    separator = ' ';

  while (buffer.nextword(token, i, separator))
  {
    if (token.starts_with('='))
      token.remove_leading_chars(1);

    int j = token.index('-');
    if (j > 0)    // looks like range specification
    {
      const_IWSubstring r1, r2;
      if (! token.split(r1, '-', r2) || 0 == r1.length() || 0 == r2.length())
      {
        std::cerr << "Min_Max_Specifier::initialise:invalid range specification '" << token << "'\n";
        return 0;
      }

      T vr1, vr2;
      if (! r1.numeric_value(vr1) || ! r2.numeric_value(vr2) || vr1 > vr2)
      {
        std::cerr << "Min_Max_Specifier::initialise:invalid range specification '" << token << "'\n";
        return 0;
      }
      while (vr1 <= vr2)    // theoretically T could be float or double type...
      {
        add(vr1);
        vr1 += static_cast<T>(1);
      }
    }
    else
    {
      T m;

      if (! token.numeric_value(m))
      {
        std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
        return 0;
      }
      add(m);
    }
  }

  return _mask_and_int._set;
}

template <typename T>
int
Matcher<T>::write_compact_representation(std::ostream & os) const
{
  if (_mask_and_int._mask[kMin]) {
    os << ">=" << _min_val;
    return os.good();
  }

  if (_mask_and_int._mask[kMax]) {
    os << "<=" << _max_val;
    return os.good();
  }

  if (! _mask_and_int._mask[kSet]) {
    return os.good();
  }

  os << _single_value;
  for (int i = 0; i < _number_other_values; i++)
  {
    os << ',';
    os << _other_values[i];
  }

  return os.good();
}


template <typename T>
Matcher<T> &
Matcher<T>::operator = (const Matcher<T> & rhs)
{
  assert(ok());
  assert(rhs.ok());

  _mask_and_int = rhs._mask_and_int;
  if (_mask_and_int._mask[kSet]) {
    _single_value = rhs._single_value;
  }
  if (_mask_and_int._mask[kMin]) {
    _min_val = rhs._min_val;
  }
  if (_mask_and_int._mask[kMax]) {
    _max_val = rhs._max_val;
  }

  if (rhs._other_values == nullptr) {
    if (rhs._other_values == nullptr) {
      return *this;
    }

    _reset_other_values();
    return *this;
  }

  // There are values in rhs.
  // See if we can just copy the values over.
  if (_number_other_values >= rhs._number_other_values) {
    std::copy_n(rhs._other_values, rhs._number_other_values, _other_values);
    _number_other_values = rhs._number_other_values;
    return *this;
  }

  delete [] _other_values;
  _other_values = new T[rhs._number_other_values];
  std::copy_n(rhs._other_values, rhs._number_other_values, _other_values);
  _number_other_values = rhs._number_other_values;
  return *this;
}

template <typename T>
bool
Matcher<T>::operator==(const Matcher<T>& rhs) const {
  if (_mask_and_int != rhs._mask_and_int) {
    return 0;
  }
  if (_mask_and_int._set[kSet] && _single_value != rhs._single_value) {
    return 0;
  }
  if (_mask_and_int._set[kMin] && _min_val != rhs._min_val) {
    return 0;
  }
  if (_mask_and_int._set[kMax] && _max_val != rhs._max_val) {
    return 0;
  }

  // If both null, we are done.
  if (_other_values == nullptr && rhs._other_values == nullptr) {
    return 1;
  }

  if (_number_other_values != rhs._number_other_values) {
    return 0;
  }

  return std::equal(_other_values, _other_values + _number_other_values, rhs._other_values);
}

template <typename T>
resizable_array<T>
Matcher<T>::ValuesMatched() const {
  resizable_array<T> result;

  if (! _mask_and_int._mask[kSet]) {
    return result;
  }

  result << _single_value;
  if (_number_other_values == 0) {
    return result;
  }

  result.add(_other_values, _number_other_values);

  return result;
}

template <typename T>
int
Matcher<T>::number_elements() const {
  if (_mask_and_int._mask[kSet] == 0) {
    return 0;
  }

  return 1 + _number_other_values;
}

template <typename T>
int
Matcher<T>::resize(int r) {
  if (r != 0) {
    std::cerr << "Matcher<T>::resize:only supported for r = 0 " << r << '\n';
    return 0;
  }

  if (_mask_and_int._mask[kSet] == 0) {
    return 1;
  }

  _mask_and_int._mask[kSet] = 0;
  _reset_other_values();

  return 1;
}

template <typename T>
int
Matcher<T>::min_explicit_value(T& v) const {
  if (_mask_and_int._mask[kSet] == 0) {
    return 0;
  }

  T result = _single_value;
  for (int i = 0; i < _number_other_values; ++i) {
    if (_other_values[i] < result) {
      result = _other_values[i];
    }
  }
  return result;
}

template <typename T>
T
Matcher<T>::operator[](int ndx) const {
  if (_mask_and_int._mask[kSet] == 0) {
    std::cerr << "Matcher::operator[] no values\n";
    ::abort();
  }

  if (ndx == 0) {
    return _single_value;
  }

  if (_number_other_values < ndx) {
    std::cerr << "Matcher::operator[] out of range\n";
    ::abort();
  }

  return _other_values[ndx - 1];
}

template <typename T>
int
Matcher<T>::write_msi(std::ostream &os, const char * name, int indentation) const {
  if (! is_set()) {
    return 1;
  }

  IWString ind;
  if (indentation) {
    ind.extend(indentation, ' ');
  }

  if (_mask_and_int._mask[kMin]) {
    os << ind << "(A I min_" << name << " " << _min_val << ")\n";
  }
  if (_mask_and_int._mask[kMax]) {
    os << ind << "(A I max_" << name << " " << _max_val << ")\n";
  }

  if (_mask_and_int._mask[kSet] == 0) {
    return os.good();
  }

  os << ind << "(A I " << name << ' ';
  if (_other_values == nullptr) {
    os << _single_value << ")\n";
    return os.good();
  }

  // Have multiple values.
  os << '(' << _single_value;
  for (int i = 0; i < _number_other_values; ++i) {
    os << ' ' << _other_values[i];
  }
  os << "))\n";

  return os.good();
}

template <typename T>
int
Matcher<T>::adjust_to_accommodate(const T v)
{
  assert (ok());

  if (matches(v))    // nothing to do!
    return 0;

// If neither min nor max are specified, we will need to add this to the archive.

  int check_archive = 1;

// Does not match, something is wrong.

  if (_mask_and_int._mask[kMin] && v < _min_val) {
    _min_val = v;
    check_archive = 0;
  }

  if (_mask_and_int._mask[kMax] && v > _max_val) {
    _max_val = v;
    check_archive = 0;
  }

  if (! check_archive)
    return 1;

  add(v);

  return 1;
}

template <typename T> template <typename C>
bool
Matcher<T>::AnyValue(C condition) const {
  if (_mask_and_int._mask[kSet] == 0) {
    return false;
  }

  if (condition(_single_value)) {
    return true;
  }

  for (int i = 0; i < _number_other_values; ++i) {
    if (condition(_other_values[i])) {
      return true;
    }
  }

  return false;
}

template <typename T> template <typename C>
void
Matcher<T>::UpdateValues(C fn) {
  if (! _mask_and_int._mask[kSet]) {
    return;
  }

  _single_value = fn(_single_value);
  for (int i = 0; i < _number_other_values; ++i) {
    _other_values[i] = fn(_other_values[i]);
  }

  return;
}


#endif  // IW_IMPLEMENTATIONS_EXPOSED || MATCHER_IMPLEMENTATION

}  // namespace iwmatcher

#endif // FOUNDATIONAL_IWMISC_MATCHER_H
