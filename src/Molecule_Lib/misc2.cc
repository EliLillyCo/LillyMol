#include <iostream>
#include <memory>

#ifdef UNIX
#include <unistd.h>
#endif

#include "Foundational/data_source/string_data_source.h"
#include "Foundational/iwmisc/iwconfig.h"

#include "misc2.h"

using std::cerr;

/*
  This file contains miscelaneous functions which cannot go in misc1.cc
  because of misc1 instantiates templates.
*/

int
iw_rename(const char* old_name, const char* new_name) {
  assert(nullptr != old_name);
  assert(nullptr != new_name);

  assert(strlen(old_name));
  assert(strlen(new_name));

  return rename(old_name, new_name);
}

int
iw_getpid() {
  return IW_GETPID();
}

/*
  An int comparitor for use with qsort.
*/

int
int_comparitor_larger(const int* pi1, const int* pi2) {
  if (*pi1 < *pi2) {
    return -1;
  } else if (*pi1 == *pi2) {
    return 0;
  } else {
    return 1;
  }
}

int
uint64_comparitor_smaller(const uint64_t* pi1, const uint64_t* pi2) {
  if (*pi1 < *pi2) {
    return 1;
  } else if (*pi1 == *pi2) {
    return 0;
  } else {
    return -1;
  }
}

void
iwxor(const int* i1, int* i2, int n) {
  for (int i = 0; i < n; i++) {
    i2[i] = (i2[i] ^ i1[i]);
  }

  return;
}

#define DEBUG_SHORTENING

/*
  This function is used in the ring finding stuff based on bonds.

  We have two sets of bonds. See if the first, I1, can be XOR'd with
  I2 to shorten I2. If so, do the XOR, and return 1
*/

/*int
try_shortening (const int * i1,
                int * i2,
                int n)
{
// First count bits

  int bits_in_i2 = 0;
  int bits_in_xor = 0;
  int bits_in_common = 0;
  for (int i = 0; i < n; i++)
  {
#ifdef DEBUG_SHORTENING
    cerr << "i = " << i << " i1 = " << i1[i] << " i2 = " << i2[i] << '\n';
#endif

    if (i2[i])
    {
      bits_in_i2++;
      if (0 == i1[i])
        bits_in_xor++;
      else
        bits_in_common++;
    }
    else if (i1[i])
    {
      bits_in_xor++;
    }
  }

#ifdef DEBUG_SHORTENING
  cerr << "i2 contains " << bits_in_i2 << " bits, xor has " << bits_in_xor << '\n';
#endif

  if (0 == bits_in_common)
    return 0;

  if (bits_in_xor >= bits_in_i2)
    return 0;

  xor (i1, i2, n);     // do the XOR operation

  return 1;
}*/

/*
  A frequent operation when parsing smiles/smarts is to extract one or
  more digits.  We return the number of characters we fully parse.
*/

int
fetch_numeric(const const_IWSubstring& zstring, int& result, int max_chars) {
  if (max_chars > zstring.length()) {
    max_chars = zstring.length();
  } else if (0 == max_chars) {
    max_chars = zstring.length();
  }

  return fetch_numeric_char(zstring.rawchars(), result, max_chars);
}

int
fetch_numeric_char(const char* zstring, int& result, int max_chars) {
  int consumed = 0;
  result = 0;

  for (int i = 0; i < max_chars; i++) {
    int tmp = zstring[i] - '0';
    if (tmp < 0 || tmp > 9) {
      return consumed;
    }

    result = result * 10 + tmp;
    consumed++;
  }

  return consumed;
}

std::optional<int>
FetchNumeric(const const_IWSubstring& s, int& ndx) {
  int result = 0;
  int consumed = 0;

  for (; ndx < s.length(); ++ndx, ++consumed) {
    const int tmp = s[ndx] - '0';
    if (tmp < 0 || tmp > 9) {
      break;
    }
    result = result * 10 + tmp;
  }

  if (consumed == 0) {
    return std::nullopt;
  }

  return result;
}

void
iwabort() {
  abort();
}

/*
  How many combinations possible choosing K items from N
  Same is N! / K! (N-K)!

  Need to evaluate the expression somewhat carefully to avoid overflow

  We determine which is larger, K or (N-K)
*/

#ifdef IW_COMBINATORIAL_COMBINATIONS_NEEDED
#include "iwfactorial.h"

static IW_Factorial<uint64_t> iwfactorial(60);

uint64_t
iw_combinatorial_combinations(int n, int k) {
  assert(n > 0);
  assert(k > 0);
  assert(n >= k);

  // Some simple boundary cases first

  if (k == n) {
    return 1;
  }

  if (1 == k) {
    return n;
  }

  int larger_denominator, smaller_denominator;

  if (k > n - k) {
    larger_denominator = k;
    smaller_denominator = n - k;
  } else {
    larger_denominator = n - k;
    smaller_denominator = k;
  }

  // Work out N! / LARGER_DENOMINATOR!

  uint64_t rc = static_cast<uint64_t>(n);

  int j = n - 1;
  while (j > larger_denominator) {
    uint64_t tmp = rc * static_cast<uint64_t>(j);

    if (tmp < rc)  // must have overflowed
    {
      cerr << "Overflow in iw_combinatorial_combinations " << n << "," << k << '\n';
      return 0;
    }

    rc = tmp;

    j--;
  }

  rc = rc / iwfactorial[smaller_denominator];

// cerr << "Dividing by " << iwfactorial[smaller_denominator] << " to yield " << rc << '\n';

// #define DEBUG_IW_COMBINATORIAL_COMBINATIONS
#ifdef DEBUG_IW_COMBINATORIAL_COMBINATIONS

  // Note this doesn't work if things overflow the factorial object

  cerr << n << " choose " << k << " is " << rc << '\n';
  uint64_t bf = (iwfactorial[n] / (iwfactorial[k] * iwfactorial[n - k]));
  cerr << "Brute force " << bf << '\n';
  assert(bf == rc);

#endif

  return rc;
}

#endif

/*
  The structure reading functions need the ability to skip over the rest
  of a bad connection table
*/

template <typename T>
int
skip_to_string(T& input, const char* target, int report_discard) {
  if (input.at_eof()) {
    return 0;
  }

  int records_discarded = 0;
  int start_record = input.lines_read();

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with(target)) {
      if (records_discarded && report_discard) {
        cerr << "skip to string: '" << target << "' discarded " << records_discarded
             << " records starting at " << start_record << '\n';
      }
      return 1;
    }

    records_discarded++;
  }

  if (report_discard) {
    cerr << "skip to string: reached EOF, no '" << target << "' found\n";
  }

  return 0;
}

template int skip_to_string<iwstring_data_source>(iwstring_data_source&, char const*,
                                                  int);
template int skip_to_string<String_Data_Source>(String_Data_Source&, char const*, int);

int
identify_plus_positions(const const_IWSubstring& buffer, resizable_array<int>& pos,
                        const char needle) {
  pos.resize_keep_storage(0);

  const int n = buffer.length();

  int in_square_bracket = 0;
  for (int i = 0; i < n; ++i) {
    if ('[' == buffer[i]) {
      in_square_bracket = 1;
    } else if (in_square_bracket) {
      if (']' == buffer[i]) {
        in_square_bracket = 0;
      }
    } else if (needle == buffer[i]) {
      pos.add(i);
    }
  }

  return pos.number_elements();
}

int
SplitOnPlusses(const const_IWSubstring& buffer,
               resizable_array_p<const_IWSubstring>& parts, const char split) {
  parts.resize_keep_storage(0);

  const int n = buffer.length();
  if (n <= 0) {
    return 0;
  }

  int cstart = 0;
  int in_square_bracket = 0;
  for (int i = 0; i < n; ++i) {
    if ('[' == buffer[i]) {
      in_square_bracket = 1;
    } else if (in_square_bracket) {
      if (']' == buffer[i]) {
        in_square_bracket = 0;
      }
    } else if (split == buffer[i]) {
      const_IWSubstring* componentPtr = new const_IWSubstring;

      if (i > 0 && i > cstart) {
        buffer.from_to(cstart, i - 1, *componentPtr);
      }

      parts.add(componentPtr);
      cstart = i + 1;
    }
  }

  // Add whatever came after the last plus.

  const_IWSubstring* lastComponentPtr = new const_IWSubstring;

  if (n > cstart) {
    buffer.from_to(cstart, n - 1, *lastComponentPtr);
  }

  parts.add(lastComponentPtr);

  return parts.number_elements();
}

namespace misc2 {

int
MatchingOpenCloseChar(const const_IWSubstring& s, int i, const char open,
                      const char close) {
  assert(s[i] == open);

  int level = 0;
  for (; i < s.length(); ++i) {
    if (s[i] == close) {
      level--;
      if (level == 0) {
        return i;
      }
    } else if (s[i] == open) {
      level++;
    }
  }

  return -1;  // Did not find matching char.
}

}  // namespace misc2
