#include <stdlib.h>
#include <memory>
#include <fstream>
#include <iostream>
#include <iomanip>
using std::cerr;
using std::endl;

#ifdef __WIN32__
  #include <winsock2.h>
  #include <limits.h>
  #include <iomanip>
  #pragma comment(lib, "Ws2_32.lib")
#else
  #include <netinet/in.h>
#endif

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwbits/dy_fingerprint.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/md5.h"
#include "Foundational/iwmisc/timsort.hpp"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iwstring.h"

#include "sparse_fp_creator.h"

Sparse_Fingerprint_Creator::Sparse_Fingerprint_Creator ()
{
  return;
};

#ifdef SEEMS_TO_NOT_WORK

int
Sparse_Fingerprint_Creator::operator == (const Sparse_Fingerprint_Creator & rhs) const
{
  if (_fp.size() != rhs._fp.size())
    return 0;

  FPHash::const_iterator i = _fp.begin();
  FPHash::const_iterator j = rhs._fp.begin();

  while (1)
  {
    if ((*i).first != (*j).first)
      return 0;

    if ((*i).second != (*j).second)
      return 0;

    ++i;
    ++j;

    if (i == _fp.end())
      return 1;

    if (j == rhs._fp.end())   // should not happen
      return 1;
  }

  return 1;
}
#endif

void
Sparse_Fingerprint_Creator::hit_bit(unsigned int b)
{
  auto f = _fp.find(b);
  if (f == _fp.end())
    _fp[b] = 1;
  else
    f->second++;
}

void
Sparse_Fingerprint_Creator::hit_bit(unsigned int b, int c)
{
  auto f = _fp.find(b);
  if (f == _fp.end())
    _fp[b] = c;
  else
    f->second += c;
}

int
Sparse_Fingerprint_Creator::operator == (const Sparse_Fingerprint_Creator & rhs) const
{
  if (_fp.size() != rhs._fp.size())
    return 0;

//cerr << "Size matches, checking " << _fp.size() << " bits\n";

//for (auto i = _fp.begin(); i != _fp.end(); ++i)
  for (FPHash::const_iterator i = _fp.begin(); i != _fp.end(); ++i)
  {
//  const auto f = rhs._fp.find((*i).first);
    const FPHash::const_iterator f = rhs._fp.find((*i).first);

    if (f == rhs._fp.end())
      return 0;

    if ((*i).second != (*f).second)
      return 0;
  }

  return 1;
}

void
Sparse_Fingerprint_Creator::copy_bits_to_unsigned_int_array(unsigned int * tmp,
                                     int & ndx) const
{
  ndx = 0;
  for (FPHash::const_iterator i = _fp.begin(); i != _fp.end(); i++)
  {
    unsigned int b = (*i).first;

    tmp[ndx] = b;
    ndx++;
  }

  if (static_cast<unsigned int>(ndx) != _fp.size())
    cerr << "Yipes, ndx " << ndx << " vs size " << _fp.size() << endl;
  assert (_fp.size() == static_cast<unsigned int>(ndx));

  return;
}

int
Sparse_Fingerprint_Creator::fill_count_array(const unsigned int * b,
                                        int * c, 
                                        int n) const
{
  int rc = n;

  for (int i = 0; i < n; i++)
  {
    FPHash::const_iterator f = _fp.find(b[i]);
    if (f == _fp.end())
    {
      cerr << "Sparse_Fingerprint_Creator::fill_count_array:no data for " << b[i] << endl;
      rc = 0;
      c[i] = 0;
    }
    else
      c[i] = (*f).second;
  }

  return rc;
}

static int
unsigned_int_comparitor(const void * p1, const void * p2)
{
  const unsigned int u1 = *(reinterpret_cast<const unsigned int *>(p1));
  const unsigned int u2 = *(reinterpret_cast<const unsigned int *>(p2));

  if (u1 < u2)
    return -1;

  if (u1 > u2)
    return 1;

  return 0;
}

class Unsigned_Int_Comparitor
{
  public:

    int operator() (unsigned int u1, unsigned int u2) const
      {
        if (u1 < u2)
          return -1;
        if (u1 > u2)
          return 1;
        return 0;
      }
};

static Unsigned_Int_Comparitor uic;

int
Sparse_Fingerprint_Creator::write_constant_width_fingerprint(unsigned int nb,
                                 const IWString & tag,
                                 std::ostream & os) const
{
  int * tmp = new_int(nb); std::unique_ptr<int[]> free_tmp(tmp);

  return _write_constant_width_fingerprint(nb, tmp, tag, os);
}

int
Sparse_Fingerprint_Creator::_write_constant_width_fingerprint(unsigned int nb,
                                 int * tmp,
                                 const IWString & tag,
                                 std::ostream & os) const
{
  for (FPHash::const_iterator i = _fp.begin(); i != _fp.end(); i++)
  {
    unsigned int b = (*i).first;

    b = b % nb;

    tmp[b]++;
  }

  IW_Bits_Base dyfp;

  (void) dyfp.construct_from_array_of_ints(reinterpret_cast<const int *>(tmp), nb);

  if (iw_little_endian())
  {
    cerr << "Sparse_Fingerprint_Creator::_write_constant_width_fingerprint: code is broken for Little endian systems, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    abort();
  }

  IWString dy_ascii;
  dyfp.daylight_ascii_representation(dy_ascii);

  os << tag << dy_ascii << ">\n";

  return os.good();
}

inline void
Sparse_Fingerprint_Creator::_convert_to_unsigned_char(unsigned int b,
                                     unsigned char & count) const
{
  FPHash::const_iterator f = _fp.find(b);

  int c = (*f).second;

  assert (c > 0);

  if (c > 255)
    count = 255;
  else
    count = static_cast<unsigned char>(c); 

  return;
}

//#define DEBUG_ENCODE_WITH_COUNTS

int
Sparse_Fingerprint_Creator::daylight_ascii_form_with_counts_encoded(IWString & dyascii) const
{
  unsigned int * s = new unsigned int[_fp.size()]; std::unique_ptr<unsigned int[]> free_s(s);

  int ndx = 0;
  copy_bits_to_unsigned_int_array(s, ndx);

  if (0 == ndx)
    return 1;

//qsort (s, ndx, sizeof(unsigned int), unsigned_int_comparitor);
//iwqsort(s, ndx, uic);
  gfx::timsort(s, s + ndx);

  int words_needed = words_needed_for_counted_form(ndx);

#ifdef DEBUG_ENCODE_WITH_COUNTS
  cerr << "To encode " << _fp.size() << " bits we need " << words_needed << " words\n";
#endif

  unsigned int * tmp = new unsigned int[words_needed]; std::unique_ptr<unsigned int[]> free_tmp(tmp);

  return _daylight_ascii_form_with_counts_encoded(s, tmp, dyascii);
}

/*
  We need to be very careful with byte swapping. We swap the bytes in those words
  that are holding numbers, but not the words that are holding the bytes that 
  correspond to the counts
*/

int
Sparse_Fingerprint_Creator::_daylight_ascii_form_with_counts_encoded(const unsigned int * s,
                                  unsigned int * tmp,
                                  IWString & dyascii) const
{
  union 
  {
    unsigned int counts;
    unsigned char c[4];
  } counts;

  int nb = static_cast<int>(_fp.size());

  int ndx = 0;
  for (int i = 0; i < nb; i += 4)
  {
    counts.counts = 0;

    int jstop = nb - i;
    if (jstop > 4)
      jstop = 4;

    for (int j = 0; j < jstop; j++)
    {
      unsigned int b = s[i + j];

      _convert_to_unsigned_char(b, counts.c[j]);

      b = htonl(b);    // byte swap if needed

      tmp[ndx] = b;
      ndx++;
    }

    tmp[ndx] = counts.counts;
    ndx++;
  }

#ifdef DEBUG_ENCODE_WITH_COUNTS
  cerr << "After building array, ndx = " << ndx << endl;
  for (int i = 0; i < ndx; i++)
  {
    cerr << tmp[i] << endl;
  }
  const unsigned char * uc = (const unsigned char *) tmp;
  for (int i = 0; i < ndx * IW_BYTES_PER_WORD; i++)
  {
    unsigned int c = uc[i];
    cerr << hex << c;
  }
  cerr << endl;
#endif

  return form_sparse_fingerprint(ndx, tmp, dyascii);
}

int
Sparse_Fingerprint_Creator::daylight_ascii_form_with_counts_encoded (const const_IWSubstring & tag,
                                                        IWString & dyascii) const
{
  IWString tmp;

  if (! daylight_ascii_form_with_counts_encoded(tmp))
    return 0;

  dyascii << tag;
  if (! tag.ends_with('<'))
    dyascii << '<';

  dyascii << tmp;

  dyascii << '>';

  return 1;
}

int
Sparse_Fingerprint_Creator::write_fingerprint (const IWString & tag,
                                               std::ostream & output) const
{
  IWString dyascii;

  (void) daylight_ascii_form_with_counts_encoded(dyascii);

  output << tag << dyascii << ">\n";

  return output.good();
}

int
Sparse_Fingerprint_Creator::write_fingerprint (const IWString & tag,
                                               IWString_and_File_Descriptor & output) const
{
  IWString dyascii;

  (void) daylight_ascii_form_with_counts_encoded(dyascii);

  output << tag << dyascii << ">\n";

  return output.good();
}

int
Sparse_Fingerprint_Creator::create_from_array_of_ints (const int * key, int n)
{
  assert (0 == _fp.size());

  for (int i = 0; i < n; i++)
  {
    if (key[i] > 0)
      _fp[i] += key[i];
  }

  return 1;
}

int
Sparse_Fingerprint_Creator::debug_print (std::ostream & os) const
{
  os << "Sparse fingerprint with " << _fp.size() << " bits set\n";
  for (const auto& iter : _fp)
  {
    os << iter.second << " hits to bit " << iter.first << endl;
  }

  return os.good();
}

template <typename O>
int
Sparse_Fingerprint_Creator::to_svml (O & output) const
{
  unsigned int * s = new unsigned int[_fp.size()]; std::unique_ptr<unsigned int[]> free_s(s);

  int ndx = 0;
  copy_bits_to_unsigned_int_array(s, ndx);

  iwqsort(s, ndx, uic);

  const auto n = _fp.size();

  for (size_t i = 0; i < n; ++i)
  {
    FPHash::const_iterator f = _fp.find(s[i]);

    const int c = (*f).second;

    output << ' ' << s[i] << ':' << c;
  }

  return 1;
}

template int Sparse_Fingerprint_Creator::to_svml(std::ofstream &) const;

int
Sparse_Fingerprint_Creator::flatten_to_01()
{
  int rc = 0;

  for (IW_Hash_Map<unsigned int, int>::iterator i = _fp.begin(); i != _fp.end(); i++)
  {
    if ((*i).second <= 1)
      continue;

    (*i).second = 1;
    rc++;
  }

  return rc;
}

int
Sparse_Fingerprint_Creator::increment_vector(int * v) const
{
  for (IW_Hash_Map<unsigned int, int>::const_iterator i = _fp.begin(); i != _fp.end(); i++)
  {
    v[(*i).first]++;
  }

  return static_cast<int>(_fp.size());
}

static IWDigits iwdigits;

int
Sparse_Fingerprint_Creator::write_in_svml_form(IWString & output) const
{
  if (0 == iwdigits.number_elements())
  {
    iwdigits.set_leading_string(":");
    iwdigits.initialise(256);
  }

  unsigned int * s = new unsigned int[_fp.size()]; std::unique_ptr<unsigned int[]> free_s(s);

  int ndx = 0;
  copy_bits_to_unsigned_int_array(s, ndx);

  if (0 == ndx)
    return 1;

  qsort(s, ndx, sizeof(unsigned int), unsigned_int_comparitor);

  for (int i = 0; i < ndx; i++)
  {
    unsigned int b = s[i];

    FPHash::const_iterator f = _fp.find(b);

    int c = (*f).second;

    if (i > 0)
      output << ' ';

    output << b;

//  cerr << "Writing bit " << b << " value " << c << endl;

    iwdigits.append_number(output, c);    
  }

  return 1;
}

template <typename O>
int 
Sparse_Fingerprint_Creator::write_as_feature_count(const char sep, O & output) const
{
  if (0 == iwdigits.number_elements())
  {
    iwdigits.set_leading_string(":");
    iwdigits.initialise(256);
  }

  unsigned int * s = new unsigned int[_fp.size()]; std::unique_ptr<unsigned int[]> free_s(s);

  int ndx = 0;
  copy_bits_to_unsigned_int_array(s, ndx);

  if (0 == ndx)
    return 1;

//qsort (s, ndx, sizeof(unsigned int), unsigned_int_comparitor);
//iwqsort(s, ndx, uic);
//gfx::timsort(s, s + ndx);

  qsort(s, ndx, sizeof(unsigned int), unsigned_int_comparitor);

  for (int i = 0; i < ndx; i++)
  {
    const unsigned int b = s[i];

    const FPHash::const_iterator f = _fp.find(b);

    int c = (*f).second;
//  cerr << "Writing bit " << b << " value " << c << endl;

    if (i > 0)
      output << sep;

    output << b;

    iwdigits.append_number(output, c);
  }

  return 1;
}

#ifdef MAY_NEED_THIS_SOMETIME
static void
do_hex_append(const unsigned char * digest,
              const int n,
              std::ostream & output)
{
  output << std::ios::hex;
  for (int i = 0; i < 16; ++i)
  {
    output << digest[i];
  }
  output << std::ios::dec;

  return;
}
#endif

template <typename O>
int 
Sparse_Fingerprint_Creator::write_as_md5_sum(O & output) const
{
  return unordered_map_to_md5(_fp, output);
}

IWString
Sparse_Fingerprint_Creator::FixedWidthFingerprint(int nbits) const {
  IW_Bits_Base bits(nbits);
  bits.clear();
  for (auto [bit, _] : _fp) {
    int b = bit % nbits;
    bits.set(b);
  }
  IWString result;
  bits.daylight_ascii_representation_including_nset_info(result);
  return result;
}

IWString
Sparse_Fingerprint_Creator::BitsWithoutCounts() const {
  unsigned int * bits = new unsigned int[_fp.size()]; std::unique_ptr<unsigned int[]> free_bits(bits);

  int ndx = 0;
  for (auto [bit, _] : _fp) {
    bits[ndx] = bit;
    ++ndx;
  }
  gfx::timsort(bits, bits + ndx);

  for (int i = 0; i < ndx; ++i) {
    bits[i] = htonl(bits[i]);
  }

  int allocated = 0;
  char * daylight = du_bin2ascii(&allocated, ndx * IW_BYTES_PER_WORD,
                                 reinterpret_cast<char *>(bits));

  IWString result;
  result.set_and_assume_ownership(daylight, allocated);
  return result;
}

#ifdef __GNUG__
//template class hashtable<pair<unsigned int const, int>, unsigned int, hash<unsigned int>, _Select1st<pair<unsigned int const, int> >, equal_to<unsigned int>, allocator<int> >;
template class IW_Hash_Map<unsigned int, int>;
template int Sparse_Fingerprint_Creator::write_as_feature_count(const char, IWString_and_File_Descriptor &) const;
template int Sparse_Fingerprint_Creator::write_as_md5_sum(IWString_and_File_Descriptor &) const;
#endif
