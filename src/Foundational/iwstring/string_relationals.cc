#include <stdlib.h>

#include "iwstring.h"

static int
common_less_than (const char * s1, 
                  int n1,
                  const char * s2,
                  int n2)
{
  int chars_to_compare = n1;
  if (n2 < n1)
    chars_to_compare = n2;

  int rc = ::strncmp (s1, s2, chars_to_compare);

  if (0 != rc)    // they differ
    return rc < 0;

// At this stage, both are the same up to the number of chars we have compared

  if (n1 == n2)   // both the same size, one is not less than the other
    return 0;

  if (n1 < n2)    // by definition, shorter string is less than
    return 1;

  return 0;
}

int
IWString::operator < (const IWString & rhs) const
{
  return common_less_than (_things, _number_elements, rhs._things, rhs._number_elements);
}

int
IWString::operator < (const const_IWSubstring & rhs) const
{
  return common_less_than (_things, _number_elements, rhs.rawchars (), rhs.length ());
}

int
const_IWSubstring::operator < (const IWString & rhs) const
{
  return common_less_than (_data, _nchars, rhs.rawchars (), rhs.length ());
}

int
const_IWSubstring::operator < (const const_IWSubstring & rhs) const
{
  return common_less_than (_data, _nchars, rhs._data, rhs._nchars);
}

static int
common_greater_than (const char * s1, 
                     int n1,
                     const char * s2,
                     int n2)
{
  int chars_to_compare = n1;
  if (n2 < n1)
    chars_to_compare = n2;

  int rc = ::strncmp (s1, s2, chars_to_compare);

  if (0 != rc)    // they differ
    return rc > 0;

// At this stage, both are the same up to the number of chars we have compared

  if (n1 == n2)   // both the same size, one is not less than the other
    return 0;

  if (n1 > n2)    // by definition, longer string is greater than
    return 1;

  return 0;
}

int
IWString::operator > (const IWString & rhs) const
{
  return common_greater_than (_things, _number_elements, rhs._things, rhs._number_elements);
}

int
IWString::operator > (const const_IWSubstring & rhs) const
{
  return common_greater_than (_things, _number_elements, rhs.rawchars (), rhs.length ());
}

int
const_IWSubstring::operator > (const IWString & rhs) const
{
  return common_greater_than (_data, _nchars, rhs.rawchars (), rhs.length ());
}

int
const_IWSubstring::operator > (const const_IWSubstring & rhs) const
{
  return common_greater_than (_data, _nchars, rhs._data, rhs._nchars);
}
