#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <iomanip>
#include <assert.h>

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwstring.h"
#include "iwbits.h"
#include "iwbits_support.h"

using namespace std;

void
IW_Bits_Base::_default_values()
{
  _nbits  = 0;

  _bits = NULL;

  _whole_bytes = 0;
  _extra_bits = 0;
}

IW_Bits_Base::IW_Bits_Base()
{
  _default_values();

  return;
}

IW_Bits_Base::IW_Bits_Base(int initial_nbits)
{
  _default_values();

  assert (0 == _nbits);

  allocate_space_for_bits(initial_nbits);

  return;
}

void
IW_Bits_Base::_copy_bits(const IW_Bits_Base & rhs)
{
  int ncopy = rhs._whole_bytes;
  if (rhs._extra_bits)
    ncopy++;

  assert (_whole_bytes >= rhs._whole_bytes);    

  memcpy(_bits, rhs._bits, ncopy);

  return;
}

IW_Bits_Base::IW_Bits_Base(const IW_Bits_Base & rhs)
{
  assert (rhs.ok());

  _default_values();

  allocate_space_for_bits(rhs._nbits);

  _copy_bits(rhs);

  return;
}

IW_Bits_Base &
IW_Bits_Base::operator =(const IW_Bits_Base & rhs)
{
  assert (rhs.ok());

  if (_nbits != rhs._nbits)
  {
    if (NULL != _bits)
      delete [] _bits;

    _default_values();

    allocate_space_for_bits(rhs._nbits);
  }

  assert (NULL != _bits);

  _copy_bits(rhs);

  return *this;
}

//#define IWB_CHECK_ALREADY_DELETED

IW_Bits_Base::~IW_Bits_Base()
{
#ifdef IWB_CHECK_ALREADY_DELETED
  if (-807 == _whole_bytes)
    cerr << "Deleting an already deleted IW_Bits_Base\n";
  _whole_bytes = -807;
#endif

  if (NULL != _bits)
  {
    delete [] _bits;
    _bits = NULL;
  }

#ifdef IWB_CHECK_ALREADY_DELETED
  _whole_bytes = _nbits = -807;
#endif

  return;
}

int
IW_Bits_Base::ok() const
{
  if (0 == _nbits && NULL == _bits && 
      0 == _whole_bytes && 0 == _extra_bits)
    return 1;

// At this stage, at least some things are non zero.

  if (_extra_bits >= IW_BITS_PER_BYTE)
    return 0;

  if (_nbits < 0)
    return 0;

  return 1;
}

int
IW_Bits_Base::debug_print(std::ostream & os) const
{
  os << "IW_Bits_Base details ";

  os << "nbits = " << _nbits << endl;

  printon(os);

  os << endl;

  return 1;
}

int
IW_Bits_Base::is_empty() const
{
  return (0 == _nbits && NULL == _bits);
}

/*
  We always deal with multiples of 8 bits
*/

int
IW_Bits_Base::allocate_space_for_bits(int nb)
{
  assert (nb > 0);

  if (NULL == _bits)
    return _allocate_space_for_bits(nb);
  else if (nb != _nbits)
    return _increase_size_for_bits(nb);

  return 1;
}

/*
  Always allocate in multiples of a word
*/

int
IW_Bits_Base::_allocate_space_for_bits(int nb)
{
  int whole_words = nb / IW_BITS_PER_WORD;

  if (0 == nb % IW_BITS_PER_WORD)
  {
    _extra_bits = 0;
  }
  else
  {
    whole_words++;
    _extra_bits  = nb % IW_BITS_PER_BYTE;
  }

  _whole_bytes = nb / IW_BITS_PER_BYTE;

  _bits = new unsigned char[whole_words * IW_BYTES_PER_WORD];

  if (NULL == _bits)
  {
    cerr << "Memory failure in IW_Bits_Base::allocate_space_for_bits\n";
    return 0;
  }

  _nbits = nb;

  memset(_bits, 0, whole_words * IW_BYTES_PER_WORD);

  return 1;
}

int
IW_Bits_Base::_increase_size_for_bits(int nb)
{
  int whole_words = nb / IW_BITS_PER_WORD;
  if (0 != nb % IW_BITS_PER_WORD)
    whole_words++;

  int new_whole_bytes = nb / IW_BITS_PER_BYTE;

  int new_extra_bits  = nb % IW_BITS_PER_BYTE;

  int bytes_needed = new_whole_bytes;

  if (new_extra_bits)
    bytes_needed++;

  unsigned char * new_bits = new unsigned char[whole_words * IW_BYTES_PER_WORD];

// Copy the old data

  int istop = _whole_bytes;
  if (_extra_bits)
    istop++;

  for (int i = 0; i < istop; i++)
  {
    new_bits[i] = _bits[i];
  }

  int last_byte = whole_words * IW_BYTES_PER_WORD;
  for (int i = istop; i < last_byte; i++)
  {
    new_bits[i] = 0;
  }

  delete [] _bits;   // done with this now

  _nbits = nb;
  _bits = new_bits;
  _whole_bytes = new_whole_bytes;
  _extra_bits = new_extra_bits;

  return 1;
}
/*
  This next code is taken from bitcount.c in libg++ 2.7.0a
*/

/* bit_count[I] is number of '1' bits in I. */
static const unsigned char
four_bit_count[16] = { 
    0, 1, 1, 2,
    1, 2, 2, 3,
    1, 2, 2, 3,
    2, 3, 3, 4};

static const unsigned char
eight_bit_count[256] = {
    0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7,
    5, 6, 6, 7, 6, 7, 7, 8
};

int
count_bits_set(const unsigned char * s, int nbytes)
{
  int rc = 0;

  for (int i = 0; i < nbytes; i++)
  {
    rc += eight_bit_count[s[i]];
  }

  return rc;
}

static unsigned int all_bits_32 = 0xffffffff;
static unsigned int all_bits_8  = 0xff;

static inline int
_BS_count_word(unsigned int word)
{
  int count = 0;
  while (word > 0)
  {
    count += four_bit_count[word & 15];
    word >>= 4;
  }
  return count;
}

/*
  Note that this function is making use of the fact that any bits beyond
  _nbits must be 0. If this ever becomes a problem, change it.
*/

int
IW_Bits_Base::nset() const
{
  assert (ok());

  int istop = _whole_bytes;
  if (_extra_bits)
    istop++;

  int rc = 0;
  const unsigned char * c = (const unsigned char *) _bits;

  for (int i = 0; i < istop; i++)
  {
    rc += eight_bit_count[*c++];
  }

  return rc;
}

bool
IW_Bits_Base::any_bits_set() const
{
  if (NULL == _bits)
    return 0;

  int istop = _whole_bytes;
  if (_extra_bits)
    istop++;

  for (int i = 0; i < istop; i++)
  {
    if (0 != _bits[i])
      return 1;
  }

  return 0;    // every byte is zero, nothing set
}


template <typename T>
int
IW_Bits_Base::construct_from_array_of_ints(const T * ii,
                                            int nb)
{
  assert (nb > 0);

  if (_nbits)      // zero out any pre-existing information
    clear();

  allocate_space_for_bits(nb);

  for (int i = 0; i < _whole_bytes; i++)
  {
    unsigned int u = _bits[i];

    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (*ii)
        u |= one_bit_8[j];
      ii++;
    }
    _bits[i] = u;
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (*ii)
        _bits[_whole_bytes] |= one_bit_8[i];
      ii++;
    }
  }

  return 1;
}

template int IW_Bits_Base::construct_from_array_of_ints(const int *, int);
template int IW_Bits_Base::construct_from_array_of_ints(const unsigned int *, int);

#ifdef OLD_VERSION_ADSKJHAKSJD

int
IW_Bits_Base::construct_from_array_of_ints(const int * ii,
                                            int nb)
{
  assert (nb > 0);

  if (_nbits)      // zero out any pre-existing information
    clear();

  allocate_space_for_bits(nb);

  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (*ii)
        _bits[i] |= one_bit_8[j];
      ii++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (*ii)
        _bits[_whole_bytes] |= one_bit_8[i];
      ii++;
    }
  }

  return 1;
}

#endif

int
IW_Bits_Base::construct_from_array_of_bits(const unsigned char * raw_bits,
                                            const int nb)
{
  assert (nb > 0);

  allocate_space_for_bits(nb);

  memcpy(_bits, raw_bits, _nbits / IW_BITS_PER_BYTE);

  return 1;
}

/*
  construct from something like '0001110101'
*/

int
IW_Bits_Base::construct_from_ascii_01_representation(const char * s, int n)
{
  allocate_space_for_bits(n);

  int eight_character_segments = n / 8;

  if (eight_character_segments)
  {
    for (int i = 0; i < eight_character_segments; i++)
    {
      unsigned char zbyte = 0;

      for (int j = 0; j < 8; j++)
      {
        if ('1' == s[j])
          zbyte |= one_bit_8[j];
        else if ('0' != s[j])
        {
          cerr << "IW_Bits_Base::construct_from_ascii_01_representation: cannot interpret bit " << (i * 8 + j) << " '" << s[j] << "'\n";
          return 0;
        }
      }

      s+= 8;
      _bits[i] = zbyte;
    }

    n = n - (eight_character_segments * 8);
  }

  for (int i = 0; i < n; i++)     // won't execute of 0 == n
  {
    if ('1' == s[i])
      _bits[_whole_bytes] |= one_bit_8[i];
  }

  return 1;
}

/*
  Construct from something like '0 1 0 0 1 1 0 1 1'.
  n is the length of the string.
  Note that there will always be an odd number of characters

  Note that we don't check that every 2nd character is a space.
*/

int
IW_Bits_Base::construct_from_ascii_01_representation_with_spaces(const char * s, int n)
{
  if (n == n / 2 * 2)
  {
    cerr << "IW_Bits_Base::construct_from_ascii_01_representation_with_spaces: must be an odd number of characters, " << n << " is invalid\n";
    return 0;
  }

  if (' ' == *s)
  {
    cerr << "IW_Bits_Base::construct_from_ascii_01_representation_with_spaces: cannot start with space\n";
    return 0;
  }

  if (' ' != s[1])
  {
    cerr << "IW_Bits_Base::construct_from_ascii_01_representation_with_spaces: 2nd char not a space\n";
    return 0;
  }

  int nb = n / 2 + 1; // two chars per bit, odd number of characters: 3 chars means 2 bits

  allocate_space_for_bits(nb);  

  int eight_character_segments = nb / IW_BITS_PER_BYTE;

//cerr << " from ascii,  n = " << n << " nb = " << nb << ", " << eight_character_segments << " 8 character segments\n";

  if (eight_character_segments)
  {
    for (int i = 0; i < eight_character_segments; i++)
    {
      unsigned char zbyte = 0;

      for (int j = 0; j < IW_BITS_PER_BYTE * 2; j += 2)
      {
        if ('1' == s[j])
          zbyte |= one_bit_8[j/2];
        else if ('0' != s[j])
        {
          cerr << "IW_Bits_Base::construct_from_ascii_01_representation_with_spaces: cannot interpret bit " << (i * 8 + j / 2) << " '" << s[j] << "'\n";
          return 0;
        }
      }

      s+= IW_BITS_PER_BYTE * 2;
      _bits[i] = zbyte;
    }

    n = n - (eight_character_segments * IW_BITS_PER_BYTE * 2);
  }

  for (int i = 0; i < n; i += 2)     // won't execute of 0 == n
  {
    if ('1' == s[i])
      _bits[_whole_bytes] |= one_bit_8[i / 2];
  }

  return 1;
}

/*int
IW_Bits_Base::construct_from_ascii_01_representation(const char * s, int n)
{
  allocate_space_for_bits(n);

  for (int i = 0; i < n; i++)
  {
    if ('1' == s[i])
      set(i);
    else if ('0' == s[i])
      ;
    else
    {
      cerr << "IW_Bits_Base::construct_from_ascii_01_representation: cannot interpret bit " << i << " '" << s[i] << "'\n";
      return 0;
    }
  }

  return 1;
}*/

/*
  We are parsing a number or range. Consume digits from the front
*/

static int
get_number(const_IWSubstring & buffer,
            int & zresult)
{
  zresult = 0;
   
  int chars_consumed = 0;

  while (buffer.length())
  {
    int j = buffer[0] - '0';

    if (j > 9 || j < 0)
      return chars_consumed;

    zresult = zresult * 10 + j;
    chars_consumed++;
    buffer.remove_leading_chars(1);
  }

  return chars_consumed;
}

/*
  The sparse representation looks like

  1,5-7,88;nbits
*/

int
IW_Bits_Base::construct_from_sparse_representation(const const_IWSubstring & bit_representation)
{
  const_IWSubstring s1, s2;

  if (! bit_representation.split(s1, ';', s2) || 0 == s2.length())
  {
    cerr << "IW_Bits_Base::construct_from_sparse_representation: no nbits '" << bit_representation << "'\n";
    return 0;
  }

  s2.strip_leading_blanks();

  int nb;
  if (! s2.numeric_value(nb) || nb < 0)
  {
    cerr << "IW_Bits_Base::construct_from_sparse_representation: invalid nbits value '" << bit_representation << "'\n";
    return 0;
  }

  allocate_space_for_bits(nb);

// Tokens must be either a number or a range

  int i = 0;
  const_IWSubstring token;
  while (s1.nextword(token, i, ','))
  {
    int n1;
    int chars_consumed = get_number(token, n1);
    if (0 == chars_consumed)
    {
      cerr << "IW_Bits_Base::construct_from_sparse_representation: invalid specifier '" << token << "'\n";
      return 0;
    }

    if (n1 < 0 || n1 >= nb)
    {
      cerr << "IW_Bits_Base::construct_from_sparse_representation: invalid bit number " << n1 << " nb = " << nb << endl;
      return 0;
    }

    if (0 == token.length())     // was just a single number
    {
      set(n1);
      continue;
    }

    if ('-' != token[0])     // must be a range
    {
      cerr << "IW_Bits_Base::construct_from_sparse_representation: invalid range '" << token << "'\n";
      return 0;
    }

    token.remove_leading_chars(1);

    int n2;
    chars_consumed = get_number(token, n2);

    if (0 == chars_consumed || token.length())
    {
      cerr << "IW_Bits_Base::construct_from_sparse_representation: invalid range '" << token << "'\n";
      return 0;
    }

    if (n2 < n1 || n2 >= nb)
    {
      cerr << "IW_Bits_Base::construct_from_sparse_representation: invalid range " << n1 << ", " << n2 << " nb = " << nb << endl;
      return 0;
    }

    set_all_bits(n1, n2, 1);
  }

  return 1;
}

/*int
IW_Bits_Base::append_sparse_form(IWString & buffer) const
{
  buffer.make_room_for_extra_items(nb);    // can never be that many

  for (int i = 0; i < nb; i++)
  {
  }

  return 1;
}*/

#include "dy_fingerprint.h"

/*
  The code to determine the number of bytes is lifted from du_ascii2bin
*/

int
IW_Bits_Base::construct_from_daylight_ascii_bit_rep(const char * ascii,
                      const int nchars)
{
  unsigned int nbytes = (nchars - 1) / 4;
  nbytes *= 3;
  int i = ascii[nchars - 1] - '0';
  nbytes -= (3 - i);

  int number_bits = nbytes * IW_BITS_PER_BYTE;

  if (! allocate_space_for_bits(number_bits))
  {
    cerr << "IW_Bits_Base::construct_from_daylight_ascii_bit_rep: cannot allocate " << number_bits << " bits\n";
    return 0;
  }

  if (0 == number_bits)
    return 1;

#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced(stderr);
#endif

  unsigned int bytes_read;
  if (! (du_ascii2bin(ascii, nchars, _bits, bytes_read)))
  {
    cerr << "IW_Bits_Base::construct_from_daylight_ascii_bit_rep: du_ascii2bin failed\n";
    cerr << "nchars = " << nchars << endl;
    return 0;
  }
#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced(stderr);
#endif

  return 1;
}

int
IW_Bits_Base::construct_from_daylight_ascii_representation(const const_IWSubstring & ascii)
{
  return construct_from_daylight_ascii_bit_rep(ascii.rawchars(), ascii.length());
}

static int
determine_nbits_nset(const const_IWSubstring & fp,
                      int i, int & nbits, int & nset)
{
  assert (';' == fp[i]);

  i++;

  const_IWSubstring token;

  if (! fp.nextword(token, i, ';') || ! fp.nextword(token, i, ';'))
  {
    cerr << "determine_nbits_nset::bad data '" << fp << "'\n";
    return 0;
  }

  if (! fp.nextword(token, i, ';') || ! token.numeric_value(nbits) || nbits <= 0)
  {
    cerr << "determine_nbits_nset::bad nbits '" << fp << "'\n";
    return 0;
  }

  if (! fp.nextword(token, i, ';') || ! token.numeric_value(nset) || nset < 0 || nset > nbits)
  {
    cerr << "determine_nbits_nset::bad nset '" << fp << "'\n";
    return 0;
  }

  return 1;
}

/*
  At some stage, try to remove this because valgrind will complain about
  reading from unitialised memory
*/

static int
determine_nbits_nset (const char * fp, int & nbits, int & nset)
{
  if (2 != IW_SSCANF(fp, ";%*d;%*d;%d;%d", &nbits, &nset))
  {
    cerr << "determine_nbits_nset: cannot determine nbits, nset '" << fp << "'\n";
    return 0;
  }

  if (nbits <= 0 || nset < 0 || nset > nbits)
  {
    cerr << "determine_nbits_nset: bad nbits (" << nbits << ") nset (" << nset << ") value '";
    while (*fp != '>')
    {
      cerr << *fp;
      fp++;
    }
    
    cerr << "', nset = " << nset << ", nbits = " << nbits << endl;
    return 0;
  }

  return 1;
}

/*
  ASCII is points to the info from the TDT. IT does NOT include
  the dataitem name
  The info about size and such starts at TDT + NCHARS
*/

int
IW_Bits_Base::_construct_from_tdt_record(const char * ascii,
                                         int nchars,
                                         int check_nset)
{
  int expected_nbits = -1;   // initialise negative so it won't be checked
  int expected_nset;
  if (check_nset)
    determine_nbits_nset(ascii + nchars, expected_nbits, expected_nset);

  if (! construct_from_daylight_ascii_bit_rep(ascii, nchars))
  {
    cerr << "IW_Bits_Base::construct_from_tdt_record: inner call failed\n";
    return 0;
  }

  if (! check_nset)
    return 1;

  if (expected_nset != nset())
  {
    cerr << "IW_Bits_Base::_construct_from_tdt_record: check nset failed\n";
    cerr << "Expected " << expected_nset << " got " << nset() << " bits set\n";
    return 0;
  }

  return 1;
}

int
IW_Bits_Base::construct_from_tdt_record(const IWString & tdt_record,
                                         int check_nset)
{
//assert (tdt_record.ends_with('>'));    // new version may include newline

  int openangle = tdt_record.index('<');

  int semicolon = tdt_record.index(';');

  return _construct_from_tdt_record(tdt_record.rawchars() + openangle + 1, 
                  semicolon - openangle - 1, check_nset);
}

int
IW_Bits_Base::construct_from_tdt_record(const const_IWSubstring & tdt_record,
                                         int check_nset)
{
//assert (tdt_record.ends_with('>'));

  int openangle = tdt_record.index('<');

  int semicolon = tdt_record.index(';');

  return _construct_from_tdt_record(tdt_record.rawchars() + openangle + 1,
              semicolon - openangle - 1, check_nset);
}

/*
  This variant is used by IW_DY_Fingerprint which needs to have nset returned
*/

int
IW_Bits_Base::construct_from_tdt_record_nset(const const_IWSubstring & tdt_record,
                      int & nset)
{
//assert (tdt_record.ends_with('>'));

  int openangle = tdt_record.index('<');

  assert (openangle > 0);

  int semicolon = tdt_record.index(';');

  const char * s = tdt_record.rawchars();

  int file_nbits;
//if (2 != IW_SSCANF(s + semicolon, ";%*d;%*d;%d;%d", &file_nbits, &nset))
  if (! determine_nbits_nset(tdt_record, semicolon, file_nbits, nset))
  {
    cerr << "IW_Bits_Base::construct_from_tdt_record_nset: cannot parse nbits/nset\n";
    cerr << tdt_record << endl;
    return 0;
  }

  assert (file_nbits >= nset);

  allocate_space_for_bits(file_nbits);

  if (! construct_from_daylight_ascii_bit_rep(s + openangle + 1,
                       semicolon - openangle - 1))
  {
    cerr << "IW_Bits_Base::construct_from_tdt_record_nset: inner call failed\n";
    return 0;
  }

#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced(stderr);
#endif

  return 1;
}

//extern "C" int bits_in_common(const unsigned int *, const unsigned int *, int);

/*
  Because we have carefully allocated words worth of storage, we can
  process the bit arrays as unsigned ints
*/

int
IW_Bits_Base::bits_in_common(const IW_Bits_Base & f2) const
{
  assert (ok());
  assert (f2.ok());
  assert (_nbits == f2._nbits);

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  return ::bits_in_common((const unsigned int *) _bits, (const unsigned int *) f2._bits, number_words);
}

/*
  Process in word chunks for efficiency
*/

void
IW_Bits_Base::iwand(const IW_Bits_Base & f2) 
{
  assert (ok());
  assert (f2.ok());
  assert (_nbits == f2._nbits);

  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) f2._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  for (int i = 0; i < number_words; i++)
  {
    w1[i] = (w1[i] & w2[i]);
  }

  return;
}

/*
  Somewhat special purpose function to AND another bit set with
  this one, and indicate whether or not the operation changed THIS

  Needs whole word optimisation
*/

void
IW_Bits_Base::iwand(const IW_Bits_Base & f2, int & changed) 
{
  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) f2._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  changed = 0;

  for (int i = 0; i < number_words; i++)
  {
    unsigned int tmp = w1[i] & w2[i];
    if (tmp != w1[i])
    {
      w1[i] = tmp;
      changed++;
    }
  }

  return;
}

void
IW_Bits_Base::iwor(const IW_Bits_Base & f2) 
{
  assert (ok());
  assert (f2.ok());
  assert (_nbits == f2._nbits);

  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) f2._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  for (int i = 0; i < number_words; i++)
  {
    w1[i] |= w2[i];
  }

  return;
}

/*
  Somewhat special purpose function to OR another bit set with
  this one, and indicate whether or not the operation changed THIS
*/

void
IW_Bits_Base::iwor(const IW_Bits_Base & f2, int & changed) 
{
  assert (ok());
  assert (f2.ok());
  assert (_nbits == f2._nbits);

  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) f2._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  changed = 0;
  for (int i = 0; i < number_words; i++)
  {
    unsigned int tmp = w1[i] | w2[i];
    if (tmp != _bits[i])
    {
      w1[i] = tmp;
      changed++;
    }
  }

  return;
}

/*
*/

void
IW_Bits_Base::iwxor(const IW_Bits_Base & f2)
{
  assert (ok());
  assert (f2.ok());
  assert (f2._nbits == _nbits);

  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) f2._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  for (int i = 0; i < number_words; i++)
  {
    w1[i] ^= w2[i];
  }

  return;
}

/*
  Needed for ring finding code. We do this by an XOR
  followed by an AND
*/

void
IW_Bits_Base::unset_bits_in_rhs(const IW_Bits_Base & f2)
{
  assert (ok());
  assert (f2.ok());
  assert (f2._nbits == _nbits);

  int istop = _whole_bytes;
  if (_extra_bits)
    istop++;

  for (int i = 0; i < istop; i++)
  {
    _bits[i] = _bits[i] & (_bits[i] ^ f2._bits[i]);
  }

  return;

  unsigned int * w1 = reinterpret_cast<unsigned int *>(_bits);
  const unsigned int * w2 = (const unsigned int *) f2._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;
  if (0 != _nbits % IW_BITS_PER_WORD)
    number_words++;

  for (int i = 0; i < number_words; i++)
  {
    w1[i] = w1[i] & (w1[i] ^ w2[i]);
  }

  return;
}

/*
  We can have fast 0,1 output if we pre-compute some of the output
*/

static const_IWSubstring bxtable []  = { "0000",
                  "0001",
                  "0010",
                  "0011",
                  "0100",
                  "0101",
                  "0110",
                  "0111",
                  "1000",
                  "1001",
                  "1010",
                  "1011",
                  "1100",
                  "1101",
                  "1110",
                  "1111" };

#define LEFT_NIBBLE  0xf0
#define RIGHT_NIBBLE 0x0f

static void
append_ascii_01(IWString & s, unsigned char c)
{
  unsigned int nibble = (c & LEFT_NIBBLE) >> 4;

#ifdef TEST_NIBBLE
  cerr << "For character '" << int(c) << "' left nibble is " << nibble;
#endif

  s += bxtable[nibble];

  nibble = (c & RIGHT_NIBBLE);

#ifdef TEST_NIBBLE
  cerr << " right nibble " << nibble << endl;
#endif

  s += bxtable[nibble];

  return;
}

int
IW_Bits_Base::append_ascii_01_representation(IWString & s) const
{
  if (s.elements_allocated() - s.number_elements() < _nbits)
    s.resize(s.number_elements() + _nbits + 1);

  const unsigned char * b = (const unsigned char *) _bits;
  for (int i = 0; i < _whole_bytes; i++)
  {
    append_ascii_01(s, b[i]);
  }

  if (_extra_bits)
  {
    const unsigned char b = _bits[_whole_bytes];

    for (int i = 0; i < _extra_bits; i++)
    {
      if (one_bit_8[i] & b)
        s += '1';
      else
        s += '0';
    }
  }

  return 1;
}


int 
print_byte(std::ostream & os, unsigned char bits, int nbits,
            const char t, const char f,
            int include_space)
{
  assert (nbits >= 0 && nbits <= IW_BITS_PER_BYTE);

  unsigned char first_bit_only = ((unsigned int) 1) << 7;

  for (int i = 0; i < nbits; i++)
  {
    if (first_bit_only == (bits & first_bit_only))
      os << t;
    else
      os << f;
    bits <<= 1;
    if (include_space)
      os << ' ';
  }

  return 1;
}


int
print_bits(std::ostream & os, const void * bits, int nbits,
            const char t, const char f,
            int include_space)
{
  const unsigned char * s = (const unsigned char *) bits;
  for (int i = 0; i < nbits; i += IW_BITS_PER_BYTE)
  {
    int bits_to_print = nbits - i;
    if (bits_to_print > IW_BITS_PER_BYTE)
      bits_to_print = IW_BITS_PER_BYTE;
    if (! print_byte(os, *s, bits_to_print, t, f, include_space))
      return 0;
    s++;
  }

  return os.good();
}

int
IW_Bits_Base::printon(std::ostream & os, const char t, const char f,
                       int include_space) const
{
//cerr << "printon writing " << _whole_bytes << " bytes\n";
  for (int i = 0; i < _whole_bytes; i++)
  {
    print_byte(os, _bits[i], IW_BITS_PER_BYTE, t, f, include_space);
  }

  if (_extra_bits)
    return print_bits(os, _bits + _whole_bytes, _extra_bits, t, f, include_space);
  else
    return os.good();
}


void
word_string_form(IWString & buffer, const unsigned int zword,
                  int bits_to_print,
                  const char t, const char f,
                  int include_space)
{
  for (int i = 0; i < bits_to_print; i++)
  {
    if (include_space)
      buffer += ' ';
    if (one_bit_32[i] & zword)
      buffer += t;
    else
      buffer += f;
  }
}

void
byte_string_form(IWString & buffer, const unsigned char zbyte,
                  int bits_to_print,
                  const char t, const char f,
                  int include_space)
{
  for (int i = 0; i < bits_to_print; i++)
  {
    if (include_space)
      buffer += ' ';
    if (one_bit_8[i] & zbyte)
      buffer += t;
    else
      buffer += f;
  }

  return;
}

/*
  Append the string form of the fingerprint to BUFFER
  Designed to avoid C++ I/O. Use instead of printon with your own buffer
*/

int
IW_Bits_Base::append_string_form(IWString & buffer, const char t, const char f,
                       int include_space) const
{
  buffer.resize(buffer.nchars() + _nbits + _nbits * (include_space != 0));

  for (int i = 0; i < _whole_bytes; i++)
  {
    byte_string_form(buffer, _bits[i], IW_BITS_PER_BYTE, t, f, include_space);
  }

  if (_extra_bits)
    byte_string_form(buffer, _bits[_whole_bytes], _extra_bits, t, f, include_space);

  return 1;
}

static int
_first_bit(unsigned int w)
{
  for (int i = 0; i < 32; i++)
  {
    if (one_bit_32[i] == (one_bit_32[i] & w))
      return i;
  }

  return -1;
}

static int
_first_bit(unsigned char c)
{
  for (int i = 0; i < 8; i++)
  {
    if (one_bit_8[i] == (one_bit_8[i] & c))
      return i;
  }

  return -1;
}

int
IW_Bits_Base::set_vector(int * vec) const
{
  assert (ok());
  assert (_nbits);

  int * iptr = vec;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
        *iptr = 1;
      else
        *iptr = 0;
      iptr++;
#ifdef DEBUG_SET_VECTOR
      if (_bits[i] & one_bit_8[j])
        cerr << "bit " << ((i * 8) + j) << " set\n";
      else
        cerr << "bit " << ((i * 8) + j) << " not set\n";
#endif
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
        *iptr = 1;
      else
        *iptr = 0;

      iptr++;
    }
  }

  return 1;
}

int
IW_Bits_Base::set_vector(int * vec, int istop) const
{
  assert (ok());
  assert (_nbits);

  int rc = 0;
  for (int i = 0; i < _whole_bytes && rc < istop; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE && rc < istop; j++)
    {
      if (_bits[i] & one_bit_8[j])
        vec[rc] = 1;
      else
        vec[rc] = 0;
      rc++;
    }
  }

  if (_extra_bits && rc < istop)
  {
    for (int i = 0; i < _extra_bits && rc < istop; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
        vec[rc] = 1;
      else
        vec[rc] = 0;

      rc++;
    }
  }

  return 1;
}

/*
  When doing population studies, we need to know which bits are turned on for
  each fingerprint.
*/

int
IW_Bits_Base::increment_vector(int * population) const
{
  assert (ok());
  assert (_nbits);

  int * iptr = population;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
        (*iptr)++;
      iptr++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
        (*iptr)++;
      iptr++;
    }
  }

  return 1;
}

int
IW_Bits_Base::increment_vector(int * population,
                                int increment) const
{
  assert (ok());
  assert (_nbits);

  int * iptr = population;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
        (*iptr) += increment;
      iptr++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
        (*iptr) += increment;
      iptr++;
    }
  }

  return 1;
}
/*
  Is an individual bit set?
  First determine in which word the bit is located.
*/

int
IW_Bits_Base::is_set(int ibit) const
{
  assert (ok());
  assert (ibit >= 0 && ibit <= _nbits);

  int ibyte = ibit / IW_BITS_PER_BYTE;

  int bit_offset = ibit % IW_BITS_PER_BYTE;

//cerr << "Bit " << ibit << " is in byte " << ibyte << " offset = " << bit_offset << endl;

//return one_bit_8[bit_offset] == (_bits[ibyte] & one_bit_8[bit_offset]);
  return (_bits[ibyte] & one_bit_8[bit_offset]);
}

void
IW_Bits_Base::set(int ibit, int new_value)
{
  assert (ok());
  assert (ibit >= 0);

  if (ibit >= _nbits)
  {
    cerr << "IW_Bits_Base::set: illegal bit " << ibit << " nbits = " << _nbits << endl;
    assert (NULL == "Death");
  }

  int ibyte = ibit / IW_BITS_PER_BYTE;

//cerr << "Setting bit " << ibit << " to " << new_value << " ibyte " << ibyte << " new value " << hex << int(_bits[ibyte]) << dec << endl;

  int bit_offset = ibit % IW_BITS_PER_BYTE;

  if (new_value)          // set the bit
    _bits[ibyte] |= one_bit_8[bit_offset];
  else
  {
    unsigned char tmp = one_bit_8[bit_offset];

    _bits[ibyte] &= ~tmp;
  }

  return;
}

void
IW_Bits_Base::clear()
{
  assert (ok());

  int nb = _whole_bytes;
  if (_extra_bits)
    nb++;

  memset(_bits, 0, nb);

  return;
}

float
IW_Bits_Base::compute_weight(const float * weight) const
{
  assert (ok());

  float rc = 0.0;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
        rc += *weight;

      weight++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
        rc += *weight;
      
      weight++;
    }
  }

  return rc;
}

int
IW_Bits_Base::compute_weight(const int * weight) const
{
  assert (ok());

  int rc = 0;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
        rc += *weight;

      weight++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
        rc += *weight;
      
      weight++;
    }
  }

  return rc;
}

/*
  Same as the basic compute weight, except that we add MISSING_BIT_VALUE
  whenever we find a bit set here, and a zero value in WEIGHT
*/

int
IW_Bits_Base::compute_weight_inc_missing(const int * weight,
                                          int missing_bit_value) const
{
  assert (ok());

  int rc = 0;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
      {
        if (*weight)
          rc += *weight;
        else
          rc += missing_bit_value;
      }

      weight++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
      {
        if (*weight)
          rc += *weight;
        else
          rc += missing_bit_value;
      }
      
      weight++;
    }
  }

  return rc;
}

float
IW_Bits_Base::compute_weight_inc_missing(const float * weight,
                                          int missing_bit_value) const
{
  assert (ok());

  float rc = 0;
  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (_bits[i] & one_bit_8[j])
      {
        if (*weight)
          rc += *weight;
        else
          rc += missing_bit_value;
      }

      weight++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (_bits[_whole_bytes] & one_bit_8[i])
      {
        if (*weight)
          rc += *weight;
        else
          rc += missing_bit_value;
      }
      
      weight++;
    }
  }

  return rc;
}

/*
  This function depends on bits off the end being zero.

  using memcmp slows things down.
*/

int
IW_Bits_Base::operator == (const IW_Bits_Base & rhs) const
{
  assert (ok());
  assert (rhs.ok());

  if (_nbits != rhs._nbits)
    return 0;

#ifdef EQ_USES_MEMCMP
  int whole_bytes = _whole_bytes;
  if (_extra_bits)
    whole_bytes++;

  return 0 == memcmp(_bits, rhs._bits, whole_bytes);
#endif

  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) rhs._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;

  for (int i = 0; i < number_words; i++)
  {
    if (w1[i] != w2[i])
      return 0;
  }

  int extra_bytes = (_whole_bytes) % IW_BYTES_PER_WORD;
  if (_extra_bits)
    extra_bytes++;

  if (0 == extra_bytes)
    return 1;

  return 0 == memcmp(_bits + number_words * IW_BYTES_PER_WORD, rhs._bits + number_words * IW_BYTES_PER_WORD, extra_bytes);
}

int
operator !=(const IW_Bits_Base & lhs, const IW_Bits_Base & rhs)
{
  return ! (lhs == rhs);
}

IW_Bits_Base
operator & (const IW_Bits_Base & lhs, const IW_Bits_Base & rhs)
{
  assert (lhs.ok());
  assert (rhs.ok());

  IW_Bits_Base rc = lhs;
  rc.iwand(rhs);

  return IW_Bits_Base(rc);
}


int
IW_Bits_Base::first_bit() const
{
  unsigned int * w1 = (unsigned int *) _bits;

  int number_words = _nbits / IW_BITS_PER_WORD;

  for (int i = 0; i < number_words; i++)
  {
    if (0 == w1[i])
      continue;

    return IW_BITS_PER_WORD * i + _first_bit(w1[i]);
  }
  
  int extra_bytes = _whole_bytes % IW_BYTES_PER_WORD;

  for (int i = 0; i < extra_bytes; i++)
  {
    unsigned char c = _bits[number_words * IW_BYTES_PER_WORD + i];

    if (0 == c)
      continue;

    return number_words * IW_BITS_PER_WORD + i * IW_BITS_PER_BYTE + _first_bit(c);
  }

  if (_extra_bits)
  {
    unsigned char c = _bits[_whole_bytes];
    if (0 == c)
      return -1;

    return _whole_bytes * IW_BITS_PER_BYTE + _first_bit(c);
  }

  return -1;

  return -1;
}

char *
IW_Bits_Base::daylight_ascii_representation(int & nchars) const
{
  assert (0 == (0x00000003 & _nbits));    // must be a multiple of 4 bits

  return du_bin2ascii(&nchars, _nbits / IW_BITS_PER_BYTE, (char *) _bits);
}

int
IW_Bits_Base::daylight_ascii_representation(IWString & result) const
{
  assert (0 == (0x00000003 & _nbits));    // must be a multiple of 4 bits

  int nchars;

  char * ascii = du_bin2ascii(&nchars, _nbits / IW_BITS_PER_BYTE, (char *) _bits);

  result.set_and_assume_ownership(ascii, nchars);

  return 1;
}

int
IW_Bits_Base::write_daylight_ascii_representation(std::ostream & output, const IWString & tag) const
{
  output << tag;
  if (! tag.ends_with('<'))
    output << '<';

  write_daylight_ascii_representation(output);

  output << ">\n";

  return output.good();
}

int
IW_Bits_Base::write_daylight_ascii_representation(std::ostream & output) const
{
  IWString zascii;

  zascii.resize(_nbits / 6 + 48);

  if (! daylight_ascii_representation(zascii))
  {
    cerr << "IW_Bits_Base::write_daylight_ascii_representation: cannot form ascii representation\n";
    return 0;

  }

  int ns = nset();

  zascii << ';' << _nbits << ';' << ns << ';' << _nbits << ';' << ns << ";1";

  output << zascii;

  return output.good();
}

int
IW_Bits_Base::daylight_ascii_representation_including_nset_info(IWString & zascii) const
{
  daylight_ascii_representation(zascii);

  int ns = nset();

  zascii << ';' << _nbits << ';' << ns << ';' << _nbits << ';' << ns << ";1";

  return zascii.length();
}

//#define DEBUG_OPERATOR_PLUSEQ

/*
  To avoid a whole bunch of messy bit manipulations we restruct operator +=
  to bits in which there are IW_BITS_PER_BYTE bits only.
*/

void
IW_Bits_Base::operator += (const IW_Bits_Base & rhs)
{
  if (0 != _nbits % IW_BITS_PER_BYTE)
  {
    cerr << "IW_Bits_Base::operator += incorrect bit count(s) " << _nbits << endl;
    assert (NULL == "Cannot do this!");
    return;
  }

#ifdef DEBUG_OPERATOR_PLUSEQ
  cerr << "Combine " << _nbits << " bits in " << _whole_bytes << " bytes, nset = " << nset() << endl;
  cerr << "   with " << rhs._nbits << " bits in " << rhs._whole_bytes << " bytes, nset = " << rhs.nset() << endl;
#endif

  int total_bits = _nbits + rhs._nbits;

  int wb = _whole_bytes;   // save the value, will be changed by allocate_space_...

  allocate_space_for_bits(total_bits);

#ifdef DEBUG_OPERATOR_PLUSEQ
  cerr << "After allocating space " << nset() << " in " << _whole_bytes << " bytes\n";
#endif

// Now append the bits from rhs. How many words to copy?

  int n2 = rhs._whole_bytes;
  if (rhs._extra_bits)
    n2++;

#ifdef DEBUG_OPERATOR_PLUSEQ
  cerr << "Appending " << n2 << " bytes from rhs starting " << wb << endl;
#endif

  for (int i = 0; i < n2; i++)
  {
    _bits[wb + i] = rhs._bits[i];
  }

#ifdef DEBUG_OPERATOR_PLUSEQ
  cerr << "After join, nbits = " << _nbits << " nset = " << nset() << endl;
  cerr << "Whole words = " << _whole_words << " extra = " << _extra_bits << endl;
#endif

  return;
}

#ifdef __i386__

static const unsigned int first_bits_32[32] = {
  0x00000080,
  0x000000c0,
  0x000000e0,
  0x000000f0,
  0x000000f8,
  0x000000fc,
  0x000000fe,
  0x000000ff,
  0x000080ff,
  0x0000c0ff,
  0x0000e0ff,
  0x0000f0ff,
  0x0000f8ff,
  0x0000fcff,
  0x0000feff,
  0x0000ffff,
  0x0080ffff,
  0x00c0ffff,
  0x00e0ffff,
  0x00f0ffff,
  0x00f8ffff,
  0x00fcffff,
  0x00feffff,
  0x00ffffff,
  0x80ffffff,
  0xc0ffffff,
  0xe0ffffff,
  0xf0ffffff,
  0xf8ffffff,
  0xfcffffff,
  0xfeffffff,
  0xffffffff
};

static const unsigned int last_bits_32[32] = {
  0x01000000,
  0x03000000,
  0x07000000,
  0x0f000000,
  0x1f000000,
  0x3f000000,
  0x7f000000,
  0xff000000,
  0xff010000,
  0xff030000,
  0xff070000,
  0xff0f0000,
  0xff1f0000,
  0xff3f0000,
  0xff7f0000,
  0xffff0000,
  0xffff0100,
  0xffff0300,
  0xffff0700,
  0xffff0f00,
  0xffff1f00,
  0xffff3f00,
  0xffff7f00,
  0xffffff00,
  0xffffff01,
  0xffffff03,
  0xffffff07,
  0xffffff0f,
  0xffffff1f,
  0xffffff3f,
  0xffffff7f,
  0xffffffff
};

#else

static const unsigned int first_bits_32[32] = {
  0x80000000,
  0xC0000000,
  0xE0000000,
  0xF0000000,
  0xF8000000,
  0xFC000000,
  0xFE000000,
  0xFF000000,
  0xFF800000,
  0xFFC00000,
  0xFFE00000,
  0xFFF00000,
  0xFFF80000,
  0xFFFC0000,
  0xFFFE0000,
  0xFFFF0000,
  0xFFFF8000,
  0xFFFFC000,
  0xFFFFE000,
  0xFFFFF000,
  0xFFFFF800,
  0xFFFFFC00,
  0xFFFFFE00,
  0xFFFFFF00,
  0xFFFFFF80,
  0xFFFFFFC0,
  0xFFFFFFE0,
  0xFFFFFFF0,
  0xFFFFFFF8,
  0xFFFFFFFC,
  0xFFFFFFFE,
  0xFFFFFFFF
};

static const unsigned int last_bits_32[32] = {
  0x00000001,
  0x00000003,
  0x00000007,
  0x0000000F,
  0x0000001F,
  0x0000003F,
  0x0000007F,
  0x000000FF,
  0x000001FF,
  0x000003FF,
  0x000007FF,
  0x00000FFF,
  0x00001FFF,
  0x00003FFF,
  0x00007FFF,
  0x0000FFFF,
  0x0001FFFF,
  0x0003FFFF,
  0x0007FFFF,
  0x000FFFFF,
  0x001FFFFF,
  0x003FFFFF,
  0x007FFFFF,
  0x00FFFFFF,
  0x01FFFFFF,
  0x03FFFFFF,
  0x07FFFFFF,
  0x0FFFFFFF,
  0x1FFFFFFF,
  0x3FFFFFFF,
  0x7FFFFFFF,
  0xFFFFFFFF
};

#endif

static const unsigned char first_bits_8[8] = {
  0x80,
  0xc0,
  0xe0,
  0xf0,
  0xf8,
  0xfc,
  0xfe,
  0xff
};

static const unsigned char last_bits_8[8] = {
  0x01,
  0x03,
  0x07,
  0x0F,
  0x1F,
  0x3F,
  0x7F,
  0xFF
};

/*
  Careful to keep any extra bits clear
*/

void
IW_Bits_Base::set_all()
{
  assert (ok());

  unsigned int * w1 = (unsigned int *) _bits;

  int number_words = _nbits / IW_BITS_PER_WORD;

  for (int i = 0; i < number_words; i++)
  {
    w1[i] = all_bits_32;
  }

  int extra_bytes = _whole_bytes % IW_BYTES_PER_WORD;

  for (int i = 0; i < extra_bytes; i++)
  {
    _bits[number_words * IW_BYTES_PER_WORD + i] = static_cast<unsigned char>(0xff);
  }

  if (_extra_bits)
    _bits[number_words * IW_BYTES_PER_WORD + extra_bytes] = first_bits_8[_extra_bits];

  return;
}

//#define DEBUG_SET_ALL_BITS

/*
  Set all the bits in a given range
*/

void
IW_Bits_Base::set_all_bits(int istart, int istop, int v)
{
#ifdef SHOW_PATTERNS
  static int first = 1;
  if (first)
  {
    for (int i = 0; i < IW_BITS_PER_BYTE; i++)
    {
      unsigned int f = first_bits_8[i];
      cerr << i << " first ";
      int nset = 0;
      for (int j = 0; j < IW_BITS_PER_BYTE; j++)
      {
        if (f & one_bit_8[j])
        {
          cerr << '1';
          nset++;
        }
        else 
          cerr << '0';

      }
      cerr << " nset = " << nset << endl;
      cerr <<   "last ";
      unsigned int l = last_bits_8[i];
      nset = 0;
      for (int j = 0; j < IW_BITS_PER_BYTE; j++)
      {
        if (l & one_bit_8[j])
        {
          nset++;
          cerr << '1';
        }
        else 
          cerr << '0';
      }
      cerr << " nset = " << nset << endl;
    }
    first = 0;
  }
#endif
  assert (istart >= 0 && istart <= _nbits);
  assert (istop >= 0 && istop <= _nbits && istop >= istart);

#ifdef DEBUG_SET_ALL_BITS
  cerr << "Setting all bits between " << istart << " and " << istop << " to " << v << endl;
  cerr << "initially " << hex << _bits[0] << ',' << _bits[1] << dec << endl;
#endif

  if (istart == istop)
  {
    set(istart, v);
    return;
  }

// Need special care with any pieces of bytes which need to be set

  int byte_start = istart / IW_BITS_PER_BYTE;
  int bit_start  = istart % IW_BITS_PER_BYTE;

  int byte_stop = istop / IW_BITS_PER_BYTE;
  int bit_stop  = istop % IW_BITS_PER_BYTE;

// Special care is needed if there is just one word being processed
// Need to create masks with the appropriate bits

  if (byte_start == byte_stop)
  {
    unsigned char mask = last_bits_8[7 - bit_start] & first_bits_8[bit_stop];
    if (0 == v)      // turning the bits off
      _bits[byte_start] = _bits[byte_start] ^ mask;
    else
      _bits[byte_start] = _bits[byte_start] | mask;

#ifdef DEBUG_SET_ALL_BITS
    cerr << "Single byte result " << hex << _bits[byte_start] << dec << endl;
#endif

    return;
  }
  
  if (bit_start)   // leading bits need to be handled
  {
    if (v)       // turn on all bits from bstart onwards
      _bits[byte_start] = _bits[byte_start] | last_bits_8[7 - bit_start];
    else         // turn off all bits from bstart onwards
      _bits[byte_start] = _bits[byte_start] & first_bits_8[bit_start - 1];

    byte_start++;
  }

// Think of the case of istop = 15; In that case, byte_stop is 1, and
// bit_stop is 7. But in order to set all bits 8-15, we need to
// process _bits[1], so increment byte_stop

  if (7 == bit_stop)
    byte_stop++;
  else
  {
    if (v)       // turn on  all bits from 0 to bstart
      _bits[byte_stop] = _bits[byte_stop] | first_bits_8[bit_stop];
    else         // turn off all bits from 0 to bit_start
      _bits[byte_stop] = _bits[byte_stop] & last_bits_8[7 - bit_stop - 1];
  }

#ifdef DEBUG_SET_ALL_BITS
  cerr << "byte_start = " << byte_start << " bit_start = " << bit_start << endl;
  cerr << "byte_stop = " << byte_stop << " bit_stop = " << bit_stop << endl;
#endif

// only fill in whole words if there are whole words to be filled in

  if (byte_stop > byte_start)
  {
    unsigned char mask;
    if (v)
      mask = all_bits_8;
    else
      mask = 0;

    for (int i = byte_start; i < byte_stop; i++)
    {
      _bits[i] = mask;
    }
  }

#ifdef DEBUG_SET_ALL_BITS
  cerr << "Result is " << hex << _bits[0] << ',' << _bits[1] << dec << endl;
#endif

  return;
}

/*
  Someone wants our bits re-arranged. We take a cross reference array and
  put the result in B2

  xref[1] = 2

  means that if bit 1 is set, set bit 2 in B2

  if any xref entry is negative, that means don't move the bit
*/

int
IW_Bits_Base::rearrange(const int * xref, IW_Bits_Base & b2) const
{
  assert (ok());
  assert (0 == _extra_bits);      // current code would overflow the XREF array - fix if ever this becomes a problem

//b2.clear();    WE ASSUME THAT B2 IS CLEARED - perhaps it isn't

  if (b2.nbits() < _nbits)
    b2.allocate_space_for_bits(_nbits);

  if (0 == _nbits)
    return 1;

  int istop = _whole_bytes;
  if (_extra_bits)
    istop++;

  for (int i = 0; i < istop; i++)
  {
    if (0 == _bits[i])    // no bits set this word
      continue;

    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      int zbit = (i * IW_BITS_PER_BYTE) + j;

      if (xref[zbit] < 0)    // don't move this bit
        continue;

      if (xref[zbit] >= _nbits)
      {
        cerr << "IW_Bits_Base::rearrange: cross reference entry out of range\n";
        cerr << "i = " << i << " xref = " << xref[i] << " nbits = " << _nbits << endl;
        return 0;
      }

      if (_bits[i] & one_bit_8[j])    // the J'th bit in the word is set
        b2.set(xref[zbit]);
    }
  }

  return 1;
}

void
IW_Bits_Base::flip()
{
  unsigned int * w1 = (unsigned int *) _bits;

  int number_words = _whole_bytes / IW_BYTES_PER_WORD;

  for (int i = 0; i < number_words; i++)
  {
    w1[i] = ~ (w1[i]);
  }

  int extra_bytes = _whole_bytes % IW_BYTES_PER_WORD;

//cerr << "Bit vector with " << _nbits << " bits, extra bytes = " << extra_bytes << endl;

  for (int i = 0; i < extra_bytes; i++)
  {
    _bits[number_words * IW_BYTES_PER_WORD + i] = ~ (_bits[number_words * IW_BYTES_PER_WORD + i]);
  }

  if (_extra_bits)
  {
    unsigned char c = _bits[number_words * IW_BYTES_PER_WORD + extra_bytes];

    _bits[number_words * IW_BYTES_PER_WORD + extra_bytes] = c ^ first_bits_8[_extra_bits];
  }

  return;
}

int
IW_Bits_Base::fold(int nfold)
{
  assert (nfold > 0);
  assert (0 == _whole_bytes % 4);
  assert (0 == _extra_bits);

  unsigned int * w = (unsigned int *) _bits;

  int number_words = _nbits / IW_BITS_PER_WORD;

  int nw2 = number_words / 2;

  for (int i = 0; i < nw2; i++)
  {
    w[i] = (w[i]) | (w[nw2 + i]);
  }

  _nbits = _nbits / 2;
  _whole_bytes = _whole_bytes / 2;

  if (1 == nfold)
    return 1;

  return fold(nfold - 1);
}

/*
  Is bit object RHS fully contained within THIS
*/

int
IW_Bits_Base::is_subset(const IW_Bits_Base & rhs) const
{
  if (_nbits != rhs._nbits)
  {
    cerr << "IW_Bits_Base::is_subset: different sized bits cannot be subsets\n";
    abort();
  }

  unsigned int * w1 = (unsigned int *) _bits;
  unsigned int * w2 = (unsigned int *) rhs._bits;

  int number_words = _nbits / IW_BITS_PER_WORD;

  for (int i = 0; i < number_words; i++)
  {
    unsigned int tmp = (w1[i]) & (w2[i]);
    if (tmp != w2[i])
      return 0;
  }

  int extra_bytes = _whole_bytes % IW_BYTES_PER_WORD;

  for (int i = 0; i < extra_bytes; i++)
  {
    unsigned char c2 = rhs._bits[number_words * IW_BYTES_PER_WORD + i];

    unsigned char c = (_bits[number_words * IW_BYTES_PER_WORD + i]) & c2;

    if (c != c2)
      return 0;
  }

  if (_extra_bits)
  {
    unsigned char c2 = rhs._bits[_whole_bytes];
    unsigned char c = (_bits[_whole_bytes]) & c2;

    if (c != c2)
      return 0;
  }

  return 1;
}

int
IW_Bits_Base::shift(int s)
{
  if (s > 0)
    return _right_shift(s);
  else if (s < 0)
    return _left_shift(-s);

  return 1;
}

/*
  For now, this is a very inefficient implementation
*/

int
IW_Bits_Base::_right_shift(int s)
{
  assert (s > 0);
  int old_nbits = _nbits;

  if (! allocate_space_for_bits(_nbits + s))
  {
    return 0;
  }

  for (int i = 0; i < old_nbits; i++)
  {
    if (is_set(old_nbits - i - 1))
      set(_nbits - i - 1, 1);
    else
      set(_nbits - i - 1, 0);
  }

  for (int i = 0; i < s; i++)
  {
    set(i, 0);
  }

  return 1;
}

int
IW_Bits_Base::_left_shift(int s)
{
  assert (s > 0);

  cerr << "Not implemented\n";
  abort();
  return 1;
}

/*
  Dec 98, added for Bob Coner.
  Someone has a contiguous chunk of memory and needs an object in that
  chunk. We copy the object and then place the bits just after the
  object. The target's pointer is adjusted and we return the
  location where the next object could be placed.
*/

void *
IW_Bits_Base::copy_to_contiguous_storage(void * p) const
{
  assert (ok());

  if (NULL == _bits)
  {
    cerr << "IW_Bits_Base::copy_to_contiguous: cannot copy empty object\n";
    return NULL;
  }

  memcpy(p, this, sizeof(IW_Bits_Base));

// Now copy the bytes

  int bytes_to_copy = _whole_bytes;
  if (_extra_bits)
    bytes_to_copy++;

  memcpy(reinterpret_cast<unsigned char *>(p) + sizeof(IW_Bits_Base), _bits, bytes_to_copy);

  IW_Bits_Base * target = reinterpret_cast<IW_Bits_Base *>(p);

  target->_bits = (unsigned char *)((unsigned char *) target + sizeof(IW_Bits_Base));

  return (reinterpret_cast<unsigned char *>(p) + sizeof(IW_Bits_Base) + bytes_to_copy);
}

void *
IW_Bits_Base::copy_to_contiguous_storage_gpu(void * p) const
{
  assert (ok());

  if (NULL == _bits)
  {
    cerr << "IW_Bits_Base::copy_to_contiguous: cannot copy empty object\n";
    return NULL;
  }

  memcpy(p, &_nbits, sizeof(int));

// Now copy the bytes

  int bytes_to_copy = _whole_bytes;
  if (_extra_bits)
    bytes_to_copy++;

  memcpy(reinterpret_cast<unsigned char *>(p) + sizeof(int), _bits, bytes_to_copy);

  return(reinterpret_cast<unsigned char *>(p) + sizeof(int) + bytes_to_copy);
}

const void *
IW_Bits_Base::build_from_contiguous_storage(const void * p,
                                             int allocate_arrays)
{
  if (allocate_arrays && NULL != _bits)
  {
    delete [] _bits;
    _nbits = 0;
    _bits = NULL;
  }

  memcpy(this, p, sizeof(IW_Bits_Base));

  p = reinterpret_cast<const unsigned char *>(p) + sizeof(IW_Bits_Base);

  int bytes_containing_bits = _whole_bytes;
  if (_extra_bits)
    bytes_containing_bits++;

  if (allocate_arrays)
  {
    allocate_space_for_bits(_nbits);

    memcpy(_bits, p, bytes_containing_bits);
  }
  else
    _bits = (unsigned char *)p;    // loss of const here, dangerous if destructor ever called

  p = reinterpret_cast<const unsigned char *>(p) + bytes_containing_bits;

  return p;
}

/*
  For every different byte pattern, what is the first bit that is set
*/

static int first_bit_set [] = {
      99,    /* this is wrong, 0 has no bits set, check for 0 first */
      7,  /* 1 */
      6,  /* 2 */
      6,  /* 3 */
      5,  /* 4 */
      5,  /* 5 */
      5,  /* 6 */
      5,  /* 7 */
      4,  /* 8 */
      4,  /* 9 */
      4,  /* 10 */
      4,  /* 11 */
      4,  /* 12 */
      4,  /* 13 */
      4,  /* 14 */
      4,  /* 15 */
      3,  /* 16 */
      3,  /* 17 */
      3,  /* 18 */
      3,  /* 19 */
      3,  /* 20 */
      3,  /* 21 */
      3,  /* 22 */
      3,  /* 23 */
      3,  /* 24 */
      3,  /* 25 */
      3,  /* 26 */
      3,  /* 27 */
      3,  /* 28 */
      3,  /* 29 */
      3,  /* 30 */
      3,  /* 31 */
      2,  /* 32 */
      2,  /* 33 */
      2,  /* 34 */
      2,  /* 35 */
      2,  /* 36 */
      2,  /* 37 */
      2,  /* 38 */
      2,  /* 39 */
      2,  /* 40 */
      2,  /* 41 */
      2,  /* 42 */
      2,  /* 43 */
      2,  /* 44 */
      2,  /* 45 */
      2,  /* 46 */
      2,  /* 47 */
      2,  /* 48 */
      2,  /* 49 */
      2,  /* 50 */
      2,  /* 51 */
      2,  /* 52 */
      2,  /* 53 */
      2,  /* 54 */
      2,  /* 55 */
      2,  /* 56 */
      2,  /* 57 */
      2,  /* 58 */
      2,  /* 59 */
      2,  /* 60 */
      2,  /* 61 */
      2,  /* 62 */
      2,  /* 63 */
      1,  /* 64 */
      1,  /* 65 */
      1,  /* 66 */
      1,  /* 67 */
      1,  /* 68 */
      1,  /* 69 */
      1,  /* 70 */
      1,  /* 71 */
      1,  /* 72 */
      1,  /* 73 */
      1,  /* 74 */
      1,  /* 75 */
      1,  /* 76 */
      1,  /* 77 */
      1,  /* 78 */
      1,  /* 79 */
      1,  /* 80 */
      1,  /* 81 */
      1,  /* 82 */
      1,  /* 83 */
      1,  /* 84 */
      1,  /* 85 */
      1,  /* 86 */
      1,  /* 87 */
      1,  /* 88 */
      1,  /* 89 */
      1,  /* 90 */
      1,  /* 91 */
      1,  /* 92 */
      1,  /* 93 */
      1,  /* 94 */
      1,  /* 95 */
      1,  /* 96 */
      1,  /* 97 */
      1,  /* 98 */
      1,  /* 99 */
      1,  /* 100 */
      1,  /* 101 */
      1,  /* 102 */
      1,  /* 103 */
      1,  /* 104 */
      1,  /* 105 */
      1,  /* 106 */
      1,  /* 107 */
      1,  /* 108 */
      1,  /* 109 */
      1,  /* 110 */
      1,  /* 111 */
      1,  /* 112 */
      1,  /* 113 */
      1,  /* 114 */
      1,  /* 115 */
      1,  /* 116 */
      1,  /* 117 */
      1,  /* 118 */
      1,  /* 119 */
      1,  /* 120 */
      1,  /* 121 */
      1,  /* 122 */
      1,  /* 123 */
      1,  /* 124 */
      1,  /* 125 */
      1,  /* 126 */
      1,  /* 127 */
      0,  /* 128 */
      0,  /* 129 */
      0,  /* 130 */
      0,  /* 131 */
      0,  /* 132 */
      0,  /* 133 */
      0,  /* 134 */
      0,  /* 135 */
      0,  /* 136 */
      0,  /* 137 */
      0,  /* 138 */
      0,  /* 139 */
      0,  /* 140 */
      0,  /* 141 */
      0,  /* 142 */
      0,  /* 143 */
      0,  /* 144 */
      0,  /* 145 */
      0,  /* 146 */
      0,  /* 147 */
      0,  /* 148 */
      0,  /* 149 */
      0,  /* 150 */
      0,  /* 151 */
      0,  /* 152 */
      0,  /* 153 */
      0,  /* 154 */
      0,  /* 155 */
      0,  /* 156 */
      0,  /* 157 */
      0,  /* 158 */
      0,  /* 159 */
      0,  /* 160 */
      0,  /* 161 */
      0,  /* 162 */
      0,  /* 163 */
      0,  /* 164 */
      0,  /* 165 */
      0,  /* 166 */
      0,  /* 167 */
      0,  /* 168 */
      0,  /* 169 */
      0,  /* 170 */
      0,  /* 171 */
      0,  /* 172 */
      0,  /* 173 */
      0,  /* 174 */
      0,  /* 175 */
      0,  /* 176 */
      0,  /* 177 */
      0,  /* 178 */
      0,  /* 179 */
      0,  /* 180 */
      0,  /* 181 */
      0,  /* 182 */
      0,  /* 183 */
      0,  /* 184 */
      0,  /* 185 */
      0,  /* 186 */
      0,  /* 187 */
      0,  /* 188 */
      0,  /* 189 */
      0,  /* 190 */
      0,  /* 191 */
      0,  /* 192 */
      0,  /* 193 */
      0,  /* 194 */
      0,  /* 195 */
      0,  /* 196 */
      0,  /* 197 */
      0,  /* 198 */
      0,  /* 199 */
      0,  /* 200 */
      0,  /* 201 */
      0,  /* 202 */
      0,  /* 203 */
      0,  /* 204 */
      0,  /* 205 */
      0,  /* 206 */
      0,  /* 207 */
      0,  /* 208 */
      0,  /* 209 */
      0,  /* 210 */
      0,  /* 211 */
      0,  /* 212 */
      0,  /* 213 */
      0,  /* 214 */
      0,  /* 215 */
      0,  /* 216 */
      0,  /* 217 */
      0,  /* 218 */
      0,  /* 219 */
      0,  /* 220 */
      0,  /* 221 */
      0,  /* 222 */
      0,  /* 223 */
      0,  /* 224 */
      0,  /* 225 */
      0,  /* 226 */
      0,  /* 227 */
      0,  /* 228 */
      0,  /* 229 */
      0,  /* 230 */
      0,  /* 231 */
      0,  /* 232 */
      0,  /* 233 */
      0,  /* 234 */
      0,  /* 235 */
      0,  /* 236 */
      0,  /* 237 */
      0,  /* 238 */
      0,  /* 239 */
      0,  /* 240 */
      0,  /* 241 */
      0,  /* 242 */
      0,  /* 243 */
      0,  /* 244 */
      0,  /* 245 */
      0,  /* 246 */
      0,  /* 247 */
      0,  /* 248 */
      0,  /* 249 */
      0,  /* 250 */
      0,  /* 251 */
      0,  /* 252 */
      0,  /* 253 */
      0,  /* 254 */
      0   /* 255 */
};

int
IW_Bits_Base::next_on_bit(int & ndx) const
{
  assert (ndx >= 0);

  if (ndx >= _nbits)
    return -1;

// Check the current byte first

  int b = ndx / IW_BITS_PER_BYTE;

  unsigned char zbyte = _bits[b];

  int istop = IW_BITS_PER_BYTE;
  if (b == _whole_bytes)
    istop = _extra_bits;

#ifdef DEBUG_NEXT_ON_BIT
  cerr << "next_on_bit::ndx " << ndx << " _whole_bytes " << _whole_bytes << " istop " << istop << " zbyte " << static_cast<int>(zbyte) << endl;
#endif

  for (int i = ndx % IW_BITS_PER_BYTE; i < istop; i++)
  {
    ndx++;
    if (one_bit_8[i] & zbyte)
      return ndx - 1;
  }

  for (int i = b + 1; i < _whole_bytes; i++)
  {
    if (0 != _bits[i])
    {
      int rc = first_bit_set[_bits[i]];

#ifdef DEBUG_NEXT_ON_BIT
      cerr << "Non zero bit, i = " << i << " rc = " << rc << " bit " << static_cast<int>(_bits[i]) << endl;
#endif

      ndx += rc + 1;
      return ndx - 1;
    }

    ndx += IW_BITS_PER_BYTE;
  }

#ifdef DEBUG_NEXT_ON_BIT
  cerr << " b = " << b << endl;
#endif

  if (b == _whole_bytes)
    return -1;

  if (_extra_bits)
  {
    zbyte = _bits[_whole_bytes];

    for (int i = 0; i < _extra_bits; i++)
    {
      ndx++;
      if (one_bit_8[i] & zbyte)
        return ndx - 1;
    }
  }

  ndx = _nbits + 1;

  return -1;
}
