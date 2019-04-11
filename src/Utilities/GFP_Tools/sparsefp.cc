#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef _WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

#include <limits>
#include <algorithm>

#define IWQSORT_FO_IMPLEMENTATION 1

#include "cmdline.h"
#include "misc.h"
#include "iwqsort.h"
#include "sparse_fp_creator.h"

#include "sparsefp.h"
#include "dyfp.h"
#include "various_distance_metrics.h"

extern unsigned char bic_table[256*256];

static int warn_empty_data = 0;

void
set_sparsefp_warn_empty_data(int s)
{
  warn_empty_data = s;
}

static double continuous_tanimoto_exponent = 1.0;

void
set_continuous_tanimoto_exponent(double s)
{
  continuous_tanimoto_exponent = s;
}

int 
Sparse_Fingerprint::check_sorted () const
{
  unsigned int prev = _bit[0];
  int need_to_die = 0;

  for (int i = 1; i < _nbits; i++)
  {
    if (_bit[i] <= prev)
    {
      cerr << "Out of order, i = " << i << " prev " << prev << " now " << _bit[i] << endl;
      need_to_die = 1;
    }

    prev = _bit[i];
  }

  if (need_to_die)
  {
    cerr << "Sparse_Fingerprint::check_sorted: bits out of order\n";
    debug_print(cerr);
    abort();
  }

  return 1;
}

Sparse_Fingerprint::Sparse_Fingerprint ()
{
  _nbits = 0;
  _nset = 0;

  _bit = NULL;
  _count = NULL;
  _sum_squared = 0;
  _norm = 0.0;

  return;
}

Sparse_Fingerprint::Sparse_Fingerprint (const Sparse_Fingerprint & rhs)
{
  _nbits = rhs._nbits;
  _nset  = rhs._nset;
  _norm  = rhs._norm;
  _sum_squared = rhs._sum_squared;

  if (NULL == rhs._bit && NULL == rhs._count)
  {
    _bit = NULL;
    _count = NULL;

    return;
  }

  assert(_nbits > 0);

  _bit = new unsigned int[_nbits];
  _count = new int[_nbits];

  copy_vector(_bit, rhs._bit, _nbits);
  copy_vector(_count, rhs._count, _nbits);

  return;
}

Sparse_Fingerprint::~Sparse_Fingerprint ()
{
  assert (-15 != _nbits);

  if (NULL != _bit)
    delete [] _bit;

  if (NULL != _count)
    delete [] _count;

  _nbits = -15;

  return;
}

int
Sparse_Fingerprint::resize (int n)
{
  if (NULL != _bit)
    delete [] _bit;

  if (NULL != _count)
    delete [] _count;

  if (0 == n)
  {
    _nbits = 0;
    _bit = NULL;
    _nset = 0;
    _count = NULL;
    return 1;
  }

  _bit = new unsigned int[n];

  _count = new int[n];

  if (NULL == _count)
  {
    cerr << "Sparse_Fingerprint::resize: memory failure for " << n << " bits\n";
    _nset = _nbits = 0;
    return 0;
  }

  _nbits = n;
  _nset = 0;

  return 1;
}

Sparse_Fingerprint &
Sparse_Fingerprint::operator = (const Sparse_Fingerprint & rhs)
{
  if (0 == rhs._nbits)
  {
    resize(0);
    return *this;
  }

  if (rhs._nbits != _nbits)
    resize(rhs._nbits);

  copy_vector(_bit, rhs._bit, _nbits);
  copy_vector(_count, rhs._count, _nbits);

  _nset = rhs._nset;
  _norm = rhs._norm;

  return *this;
}

int
Sparse_Fingerprint::construct_from_tdt_record (const const_IWSubstring & buffer)
{
  const_IWSubstring daylight = buffer;

//assert (buffer.ends_with ('>'));   // new version may have newline

  daylight.remove_up_to_first('<');
  daylight.chop();

  if (0 == daylight.length())
  {
    if (warn_empty_data)
      cerr << "Sparse_Fingerprint::construct_from_tdt_record: empty dataitem\n";
    resize(0);
    return 1;
  }

  return construct_from_daylight_ascii_representation(daylight);
}

int
Sparse_Fingerprint::debug_print (std::ostream & os) const
{
  os << "Sparse fingerprint with " << _nbits << " bits\n";
  for (int i = 0; i < _nbits; i++)
  {
    os << ' ' << i << " bit " << _bit[i];
    if (NULL != _count)
      os << " hit " << _count[i] << " times";
    os << endl;
  }

  return os.good();
}

similarity_type_t
Sparse_Fingerprint::tanimoto(const Sparse_Fingerprint & rhs) const
{
  return _tanimoto_with_counts(rhs);
}

similarity_type_t
Sparse_Fingerprint::distance (const Sparse_Fingerprint & rhs) const
{
  return static_cast<similarity_type_t>(1.0) - tanimoto(rhs);
}

similarity_type_t
Sparse_Fingerprint::optimistic_distance (const Sparse_Fingerprint & rhs,
                                         const Tversky & tv) const
{
  similarity_type_t d = distance(rhs);

  similarity_type_t tv1 = tversky_distance(rhs, tv);
  similarity_type_t tv2 = rhs.tversky_distance(*this, tv);

  if (d < tv1 && d < tv2)
    return d;

  if (tv1 < d && tv1 < tv2)
    return tv1;

  return tv2;
}

similarity_type_t
Sparse_Fingerprint::tversky_distance (const Sparse_Fingerprint & rhs, 
                                      const Tversky & tv) const
{
  if (tv.treat_non_colliding_as_01())
    return tversky_distance01(rhs, tv);

  return static_cast<similarity_type_t>(1.0) - tversky(rhs, tv);
}

/*
  We don't take any advantage of the fact that we have counts
*/

similarity_type_t
Sparse_Fingerprint::tversky_distance01 (const Sparse_Fingerprint & rhs,
                                        const Tversky & tv) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (_nbits == rhs._nbits)
      return static_cast<similarity_type_t>(0.0);

    return static_cast<similarity_type_t>(1.0);
  }

  const unsigned int * b1;
  int n1;
  const unsigned int * b2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    b1 = _bit;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    b1 = rhs._bit;
    b2 = _bit;
  }
    
  int bits_in_common = 0;

  for (int i = 0, j = 0; i < n1; i++)
  {
    unsigned int b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      bits_in_common++;
      j++;
    }
  }

  int just_a = _nbits - bits_in_common;
  int just_b = rhs._nbits - bits_in_common;

  similarity_type_t rc = static_cast<double>(bits_in_common) / static_cast<double>(tv.a() * just_a + tv.b() * just_b + bits_in_common);

  return static_cast<similarity_type_t>(1.0) - rc;          // convert to distance
}

similarity_type_t
Sparse_Fingerprint::tversky (const Sparse_Fingerprint & rhs,
                             const Tversky & tv) const
{
  if (tv.treat_non_colliding_as_01())
    return static_cast<similarity_type_t>(1.0) - tversky_distance01(rhs, tv);

  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (_nbits == rhs._nbits)    // both 0
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

#ifdef DEBUG_SPARSE_TVERSKY
  cerr << "Lhs contains " << _nbits << " bits, last " << _bit[_nbits - 1] <<
          " rhs contains " << rhs._nbits << " bits, last " << rhs._bit[rhs._nbits - 1] << endl;
#endif
          
  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    b1 = _bit;
    c1 = _count;
    b2 = rhs._bit;
    c2 = rhs._count;
  }
  else
  {
    n1 = rhs._nbits;
    b1 = rhs._bit;
    c1 = rhs._count;
    b2 = _bit;
    c2 = _count;
  }

#ifdef DEBUG_SPARSE_TVERSKY
  cerr << "LHS\n";
  for (int i = 0; i < _nbits; i++)
  {
    cerr << " i = " << i << " bit " << _bit[i] << " count " << _count[i] << endl;
  }
  cerr << "RHS\n";
  for (int i = 0; i < rhs._nbits; i++)
  {
    cerr << " i = " << i << " bit " << rhs._bit[i] << " count " << rhs._count[i] << endl;
  }

  cerr << "Scanning " << n1 << " bits\n";

#endif
    
  int bits_in_common = 0;

  register unsigned int b;   // declare here for efficiency
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
//    cerr << "  skipping j = " << j << " bit " << b2[j] << " need " << b << endl;
      j++;
    }

    if (b2[j] == b)
    {
      if (c1[i] < c2[j])
        bits_in_common += c1[i];
      else
        bits_in_common += c2[j];

      j++;
    }
  }

  int just_a = _nset - bits_in_common;
  int just_b = rhs._nset - bits_in_common;

  similarity_type_t rc = static_cast<double>(bits_in_common) / static_cast<double>(tv.a() * just_a + tv.b() * just_b + bits_in_common);

  return rc;
}

std::ostream &
operator << (std::ostream & os, const Sparse_Fingerprint & sfp)
{
  sfp.debug_print(os);

  return os;
}

int
Sparse_Fingerprint::is_set (unsigned int b) const
{
#ifdef IS_SET_SLOW
  for (int i = 0; i < _nbits; i++)
  {
    if (b == _bit[i])
      return 1;
  }

  return 0;
#endif

#ifdef DEBUG_IS_SET
  cerr << "Looking for bit " << b << endl;
#endif

  if (0 == _nbits)
    return 0;

  if (b < _bit[0])
    return 0;

  if (b > _bit[_nbits - 1])
    return 0;

  if (b == _bit[0] || b == _bit[_nbits - 1])
    return 1;

  int left = 0;
  int right = _nbits - 1;
  int middle;
  while ( (middle = (left + right) / 2) > left)
  {
    unsigned int mb = _bit[middle];

#ifdef DEBUG_IS_SET
    cerr << "left " << left << " middle " << middle << " bit " << mb << " right " << right << endl;
#endif

    if (b < mb)
      right = middle;
    else if (b > mb)
      left = middle;
    else
      return 1;
  }

  return 0;
}

similarity_type_t
Sparse_Fingerprint::_tanimoto_with_counts(const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_SPARSE_TANIMOTO_COUNT
#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int bits_in_common = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      bits_in_common += bic_table[c1[i] * 256 + c2[j]];
/*    if (c1[i] <= c2[j])
        bits_in_common += c1[i];
      else 
        bits_in_common += c2[j];*/

      j++;
    }
  }

  similarity_type_t rc = static_cast<similarity_type_t>(bits_in_common) / static_cast<similarity_type_t>(_nset + rhs._nset - bits_in_common);

#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << _nset << " bits and " << rhs._nset << " bits, " << bits_in_common << " bits_in_common, similarity " << rc << endl;
#endif

  return rc;
}

/*
  Version that does not depend on bic_table
*/

similarity_type_t
Sparse_Fingerprint::tanimoto_with_unlimited_counts(const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_SPARSE_TANIMOTO_COUNT
#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int bits_in_common = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (*b2 < b)
    {
      b2++;
    }

    if (*b2 == b)
    {
      if (c1[i] <= c2[j])
        bits_in_common += c1[i];
      else 
        bits_in_common += c2[j];

      j++;
    }
  }

  similarity_type_t rc = static_cast<similarity_type_t>(bits_in_common) / static_cast<similarity_type_t>(_nset + rhs._nset - bits_in_common);

#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << _nset << " bits and " << rhs._nset << " bits, " << bits_in_common << " bits_in_common, similarity " << rc << endl;
#endif

  return rc;
}

similarity_type_t
Sparse_Fingerprint::tanimoto_binary(const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_SPARSE_TANIMOTO_BINARY
#ifdef DEBUG_SPARSE_TANIMOTO_BINARY
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  int n1;
  const unsigned int * b2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    b1 = _bit;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    b1 = rhs._bit;
    b2 = _bit;
  }

  int bits_in_common = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; ++i)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      ++j;
    }

    if (b2[j] == b)
    {
      bits_in_common++;

      ++j;
    }
  }

  similarity_type_t rc = static_cast<similarity_type_t>(bits_in_common) / static_cast<similarity_type_t>(_nbits + rhs._nbits - bits_in_common);

#ifdef DEBUG_SPARSE_TANIMOTO_BINARY
  cerr << _nset << " bits and " << rhs._nset << " bits, " << bits_in_common << " bits_in_common, similarity " << rc << endl;
#endif

  return rc;
}

#ifdef OLD_VERSION

similarity_type_t
Sparse_Fingerprint::_tanimoto_with_counts (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_SPARSE_TANIMOTO_COUNT
#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  int n1;
  const unsigned int * b2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int bits_in_common = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      bits_in_common += bic_table[c1[i] * 256 + c2[j]];
//    if (c1[i] > c2[j])
//      bits_in_common += c2[j];
//    else                          // less than or equal
//      bits_in_common += c1[i];

      j++;
    }
  }

  similarity_type_t rc = static_cast<similarity_type_t>(bits_in_common) / static_cast<similarity_type_t>(_nset + rhs._nset - bits_in_common);

#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << _nset << " bits and " << rhs._nset << " bits, " << bits_in_common << " bits_in_common, similarity " << rc << endl;
#endif

  return rc;
}

#endif

#ifdef OLDER_VERSION
similarity_type_t
Sparse_Fingerprint::_tanimoto_with_counts (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_SPARSE_TANIMOTO_COUNT
#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  const unsigned int * b1end = b1 + n1;

  int bits_in_common = 0;

  register unsigned int b;
  while (b1 < b1end)
  {
    b = *b1;
    while (*b2 < b)
    {
      b2++;
      c2++;
    }

    if (*b2 == b)
    {
      bits_in_common += bic_table[(*c1) * 256 + (*c2)];
//    if (c1[i] > c2[j])
//      bits_in_common += c2[j];
//    else                          // less than or equal
//      bits_in_common += c1[i];

      b2++;
      c2++;
    }

    b1++;
    c1++;
  }

  similarity_type_t rc = static_cast<similarity_type_t>(bits_in_common) / static_cast<similarity_type_t>(_nset + rhs._nset - bits_in_common);

#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << _nset << " bits and " << rhs._nset << " bits, " << bits_in_common << " bits_in_common, similarity " << rc << endl;
#endif

  return rc;
}

#endif

/*
  At some stage, change this to just call dot_product and divide...
*/

similarity_type_t
Sparse_Fingerprint::cosine_measure (const Sparse_Fingerprint & rhs) const
{
  assert (NULL != _count);

//#define DEBUG_SPARSE_TANIMOTO_COUNT
#ifdef DEBUG_SPARSE_TANIMOTO_COUNT
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int d1d2 = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      d1d2 += c1[i] * c2[j];

      j++;
    }
  }

  similarity_type_t rc = static_cast<double>(d1d2) / (_norm * rhs._norm);

  return rc;
}

int
Sparse_Fingerprint::dot_product (const Sparse_Fingerprint & rhs) const
{
  assert (NULL != _count);

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int d1d2 = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      d1d2 += c1[i] * c2[j];

      j++;
    }
  }

  return d1d2;
}

/*
  Figuring out how many bits we have encoded is tricky. Here's a table of number of
  bits and number of words needed
    1  2
    2  3
    3  4 
    4  5
    5  7
    6  8
    7  9
    8 10
    9 12
   10 13
   11 14
   12 15
   13 17
   14 18
   15 19
   16 20
   17 22
*/

int
Sparse_Fingerprint::construct_from_daylight_ascii_representation(const const_IWSubstring & daylight)
{
  if (0 == daylight.length())
    return 1;

  IW_Bits_Base fp;

  if (! fp.construct_from_daylight_ascii_representation(daylight))
  {
    cerr << "Sparse_Fingerprint::_counted_form_construct_from_daylight_ascii_representation: cannot parse\n";
    cerr << daylight << endl;
    return 0;
  }

  int nb = fp.nbits();

  if (0 != nb % IW_BITS_PER_WORD)
  {
    cerr << "Sparse_Fingerprint::construct_from_daylight_ascii_representation: must be a multiple of " << IW_BITS_PER_WORD << " bits, " << nb << " is not\n";
    return 0;
  }

  if (0 == nb)
  {
    cerr << "Sparse_Fingerprint::construct_from_daylight_ascii_representation:no bits present! '" << daylight << "'\n";
    return 0;
  }

  int nwords = nb / IW_BITS_PER_WORD;

  int number_bits = nwords / 5 * 4;     // each 4 "bits" have a word with their counts

  int remainder = nwords % 5;

  if (0 != remainder)
    number_bits += remainder - 1;

//cerr << "Reading " << nb << " bits, which is " << nwords << " words. nbits " << number_bits << endl;

  resize(number_bits);

  if (! iw_little_endian())
    return _counted_form_construct_from_array_of_bits(fp.bits());

// need to do a byte swap. We only do byte swaps on the bit numbers, not the counts
    
  unsigned int * b = reinterpret_cast<unsigned int *>(const_cast<unsigned char *>(fp.bits()));

  unsigned int last_word = b[nwords - 1];

  for (int i = 0; i < nwords; i++)
  {
    if (4 != i % 5)              // swap words 0 1 2 3   5 6 7 8  10 11 12 13   15 16 ...
      b[i] = ntohl(b[i]);
  }

  b[nwords - 1] = last_word;

  return _counted_form_construct_from_array_of_bits(fp.bits());   // has been swapped
}

int
Sparse_Fingerprint::_counted_form_construct_from_array_of_bits (const void * voidb)
{
  assert (_nbits > 0);

#ifdef __GNUG__
  const unsigned int * b = static_cast<const unsigned int *>(voidb);
#else
  const unsigned int * b = reinterpret_cast<const unsigned int *>(voidb);
#endif

  union foo
  {
    unsigned int zbit;
    unsigned char count[IW_BYTES_PER_WORD];
  };

#ifdef __GNUG__
  const foo * fooptr = static_cast<const foo *>(voidb);
#else
  const foo * fooptr = reinterpret_cast<const foo *>(voidb);
#endif

//cerr << "Reading " << _nbits << " bits\n";

// This constant will be used for checking whether or not we have 4 words + a count word

  int nbits_minus_5 = _nbits - IW_BYTES_PER_WORD - 1;

  for (int i = 0, j = 0; i < _nbits; i += IW_BYTES_PER_WORD, j += IW_BYTES_PER_WORD + 1)
  {
    if (i <= nbits_minus_5)    // hopefully the most common case, we have 5 words to process
    {
      memcpy(_bit + i, b + j, IW_BYTES_PER_WORD * IW_BYTES_PER_WORD);
      const foo * c = fooptr + j + 4;
      _count[i] = c->count[0];
      _count[i + 1] = c->count[1];
      _count[i + 2] = c->count[2];
      _count[i + 3] = c->count[3];

//    cerr << "Read bits " << _bit[i] << ',' << _bit[i + 1] << ',' << _bit[i + 2] << ',' << _bit[i + 3] << endl;

      continue;
    }

//  We are at the end, we don't have a full set of words to process

    int extra_words = _nbits - i;

    memcpy(_bit + i, b + j, extra_words * IW_BYTES_PER_WORD);
    const foo * c = fooptr + j + extra_words;

    for (int k = 0; k < extra_words; k++)
    {
      _count[i + k] = c->count[k];
    }

    break;
  }

  _nset = sum_vector(_count, _nbits);
  
  _sum_squared = 0;
  for (int i = 0; i < _nbits; i++)
  {
    _sum_squared += _count[i] * _count[i];
  }

  _norm = sqrt(static_cast<double>(_sum_squared));

  check_sorted();

//#define ECHO_COUNTED_FINGERPRINTS
#ifdef ECHO_COUNTED_FINGERPRINTS
  cerr << "Just read fingerprint\n";
  debug_print(cerr);
#endif
  return 1;
}

int
Sparse_Fingerprint::next_bit_set (int & istart,
                                  unsigned int & zbit,
                                  int & hits) const
{
  if (istart >= _nbits)
    return 0;

  zbit = _bit[istart];
  hits = _count[istart];

  istart++;

  return 1;
}

int
Sparse_Fingerprint::_build_bit(int ndx,
                               const const_IWSubstring & s)
{
  const char sep = ',';

  int i = 0;
  const_IWSubstring token;

  s.nextword(token, i, sep);

  unsigned int b;

  if (! token.numeric_value(b))
  {
    cerr << "Sparse_Fingerprint::_build_bit:invalid bit '" << token << "'\n";
    return 0;
  }

  _bit[ndx] = b;

  if (! s.nextword(token, i, sep))
  {
    _count[ndx] = 1;
    return 1;
  }

  unsigned int c;
  if (! token.numeric_value(c))
  {
    cerr << "Sparse_Fingerprint::_build_bit:invalid count '" << token << "'\n";
    return 0;
  }

  _count[ndx] = c;

  return 1;
}

/*
  Somewhat hard to implement because we need to know the number of bits set up front.
*/

int
Sparse_Fingerprint::construct_from_sparse_ascii_representation (const const_IWSubstring & fp)
{
  if (0 == fp.length())
  {
    _nbits = 0;
    return 1;
  }

  int n = fp.nwords();

  if (! resize(n))
    return 0;

  const_IWSubstring token;
  int i = 0;
  int ndx = 0;
  while (fp.nextword(token, i))
  {
    if (! _build_bit(ndx, token))
    {
      cerr << "Sparse_Fingerprint::construct_from_sparse_ascii_representation:invalid bit/count specification '" << token << "'\n";
      return 0;
    }

//  cerr << "Bit " << ndx << " set to " << _bit[ndx] << " count " << _count[ndx] << endl;

    _nset += _count[ndx];

    ndx++;
  }

  assert (ndx == _nbits);

  return 1;
}

/*
  This is incorrect. With sparse fingerprints, there really is no concept of
  number of bits not set by either entity.
  Also, doesn't look like I'm correctly counting the bits either - trailing
  bits from longer vector omitted.
  Don't use this.
*/

similarity_type_t
Sparse_Fingerprint::fvb_modified_tanimoto (const Sparse_Fingerprint & rhs) const
{
//cerr << "NBITS " << _nbits << " nset " << _nset << " and rhs nbits " << rhs._nbits << " set " << rhs._nset << endl;

  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  const unsigned int * b1;
  int n1;
  const int * c1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    b1 = _bit;
    c1 = _count;
    b2 = rhs._bit;
    c2 = rhs._count;
  }
  else
  {
    n1 = rhs._nbits;
    b1 = rhs._bit;
    c1 = rhs._count;
    b2 = _bit;
    c2 = _count;
  }
    
  int n11 = 0;
  int n00 = 0;

  for (int i = 0, j = 0; i < n1; i++)
  {
    unsigned int b = b1[i];
    while (b2[j] < b)
    {
      j++;
      n00 += c2[j];
    }

    if (b2[j] == b)
    {
      if (c1[i] > c2[j])
      {
        n11 += c2[j];
      }
      else                          // less than or equal
      {
        n11 += c1[i];
      }
      j++;
    }
  }

  return fligner_verducci_blower(_nset + rhs._nset, _nset, rhs._nset, n00, n11);
}

/*
*/

similarity_type_t
Sparse_Fingerprint::manhattan_distance (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_MANHATTAN_MEASURE
#ifdef DEBUG_MANHATTAN_MEASURE
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;
  int n2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    n2 = rhs._nbits;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    n2 = _nbits;
    c2 = _count;
    b2 = _bit;
  }

  unsigned int d = 0;

  register unsigned int b;
  int j = 0;
  for (int i = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      d += c2[j];
      j++;
    }

    if (b2[j] == b)   // bit set in both
    {
      if (c1[i] > c2[j])
        d += c1[i] - c2[j];
      else
        d += c2[j] - c1[i];

      j++;
    }
    else          // bit not set in 2nd FP
      d += c1[i];
  }

// Remember, the *2 array is longer than the *1 array. Gather any
// unused items from *2

  while (j < n2)
  {
    d += c2[j];
    j++;
  }

#ifdef DEBUG_MANHATTAN_MEASURE
  cerr << "d = " << d << " becomes " << (1.0 /static_cast<double>(1 + d) ) << endl;
#endif

  return static_cast<similarity_type_t>(1.0 / static_cast<double>(1 + d));
}

/*
  Had to make a decision about how to handle the case where a bit is missing in one of the
  molecules, but not in the other.
  Do we add 1.0 to the distance, or just ignore it. I get better looking distances when I
  ignore those bits...

*/

similarity_type_t
Sparse_Fingerprint::soergel_variant_similarity (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

//#define DEBUG_SOERGEL_MEASURE
#ifdef DEBUG_SOERGEL_MEASURE
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  double rc = 0.0;
  int p = 0;

  register unsigned int b;
  int j = 0;
  for (int i = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      rc += 1.0;
      p++;
      j++;
    }

    if (b2[j] == b)   // bit set in both
    {
      if (c1[i] > c2[j])
      {
        rc += static_cast<double>(c1[i] - c2[j]) / static_cast<double>(c1[i]);
        p++;
      }
      else if (c1[i] == c2[j])
        p++;
      else
      {
        rc += static_cast<double>(c2[j] - c1[i]) / static_cast<double>(c2[j]);
        p++;
      }

      j++;
    }
    else
    {
      rc += 1.0;
      p++;
    }
  }

#ifdef DEBUG_SOERGEL_MEASURE
  cerr << "p = 0, must be identical\n";
#endif

  if (0 == p)
    return 1.0f;

#ifdef DEBUG_SOERGEL_MEASURE
  cerr << "rc = " << rc << " across " << p << " bits, becomes " << (rc / static_cast<double>(p)) << endl;
#endif

  return 1.0f - static_cast<similarity_type_t>(rc / static_cast<double>(p));
}

//#define DEBUG_SOERGEL_MEASURE


similarity_type_t
Sparse_Fingerprint::soergel_similarity (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
  {
    if (0 == _nbits && 0 == rhs._nbits)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

#ifdef DEBUG_SOERGEL_MEASURE
  cerr << "Comparing non colliding counted fingerprint with " << _nbits << " bits " << _nset << " set\n";
  int bic = 0;
  for (int i = 0; i < _nbits; i++)
  {
    unsigned int b = _bit[i];
    cerr << ' ' << i << " bit " << b << " set " << _count[i] << " times";
    if (rhs.is_set(b))
    {
      bic++;
      cerr << " *";
    }
    cerr << endl;
  }
  cerr << bic << " bits in common, and       " << rhs;
#endif

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int numerator = 0;
  int denominator = 0;

  register unsigned int b;
  int j = 0;
  for (int i = 0; i < n1; i++)
  {
    b = b1[i];

#ifdef DEBUG_SOERGEL_MEASURE
    cerr << "B1 " << b << " (i = " << i << "), c = " << c1[i] << "\n";
#endif

    while (b2[j] < b)
    {
      numerator += c2[j];
      denominator += c2[j];
#ifdef DEBUG_SOERGEL_MEASURE
      cerr << " incrementing second vector, b = " << b2[j] << " count " << c2[j] << endl;
#endif
      j++;
    }

    if (b2[j] == b)   // bit set in both
    {
#ifdef DEBUG_SOERGEL_MEASURE
      cerr << " match for bit " << b << ", i = " << i << " c1 " << c1[i] << ", j = " << j << " c2 " << c2[j] << endl;
#endif

      if (c1[i] > c2[j])
      {
        numerator += c1[i] - c2[j];
        denominator += c1[i];
      }
      else
      {
        numerator += c2[j] - c1[i];
        denominator += c2[j];
      }

      j++;
    }
    else
    {
      numerator += c1[i];
      denominator += c1[i];
#ifdef DEBUG_SOERGEL_MEASURE
      cerr << " bit 1 " << b << " ahead, c = " << c1[i] << endl;
#endif
    }
#ifdef DEBUG_SOERGEL_MEASURE
    cerr << " finished i = " << i << ", numerator " << numerator << " denominator " << denominator << endl;
#endif
  }

  return static_cast<float>(1.0f) - static_cast<similarity_type_t>(numerator) / static_cast<similarity_type_t>(denominator);
}

#ifdef MIGHT_BE_SLIGHTLY_SLOWER

int
Sparse_Fingerprint::bits_in_common (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
    return 0;

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int bits_in_common = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      bits_in_common += bic_table[c1[i] * 256 + c2[j]];

      j++;
    }
  }

  return bits_in_common;
}
#endif

int
Sparse_Fingerprint::bits_in_common (const Sparse_Fingerprint & rhs) const
{
  if (0 == _nbits || 0 == rhs._nbits)
    return 0;

  assert (NULL != _count);
  assert (_nset > 0 && rhs._nset > 0);

  const unsigned int * b1;
  const int * c1;
  int n1;
  const unsigned int * b2;
  const int * c2;

// Make sure the last item in b2 is greater than the last item in b1

  if (_bit[_nbits - 1] < rhs._bit[rhs._nbits - 1])
  {
    n1 = _nbits;
    c1 = _count;
    b1 = _bit;
    c2 = rhs._count;
    b2 = rhs._bit;
  }
  else
  {
    n1 = rhs._nbits;
    c1 = rhs._count;
    b1 = rhs._bit;
    c2 = _count;
    b2 = _bit;
  }

  int bits_in_common = 0;

  register unsigned int b;
  for (int i = 0, j = 0; i < n1; i++)
  {
    b = b1[i];
    while (b2[j] < b)
    {
      j++;
    }

    if (b2[j] == b)
    {
      if (c1[i] < c2[j])
        bits_in_common += c1[i];
      else
        bits_in_common += c2[j];

      j++;
    }
  }

  return bits_in_common;
}

int
Sparse_Fingerprint::truncate_counts_at (int c)
{
  int rc = 0;

  for (int i = 0; i < _nbits; i++)
  {
    if (_count[i] > c)
    {
      _count[i] = c;
      rc++;
    }
  }

  return rc;
}

/*
  compute fp1 - fp2
*/

void
Sparse_Fingerprint::vector_difference (const Sparse_Fingerprint & fp1, 
                                       const Sparse_Fingerprint & fp2)
{
  resize(fp1._nbits + fp2._nbits);   // worst case, no bits in common

  _nbits = 0;

  int i1 = 0;
  int i2 = 0;

  while (i1 < fp1._nbits && i2 < fp2._nbits)
  {
    if (fp2._bit[i2] < fp1._bit[i1])
    {
      _bit[_nbits] = fp2._bit[i2];
      _count[_nbits] = - fp2._count[i2];
      i2++;
    }
    else if (fp1._bit[i1] == fp2._bit[i2])
    {
      _bit[_nbits] = fp1._bit[i1];
      _count[_nbits] = fp1._count[i1] - fp2._count[i2];
      i1++;
      i2++;
    }
    else
    {
      _bit[_nbits] = fp1._bit[i1];
      _count[_nbits] = fp1._count[i1];
      i1++;
    }

    _nbits++;
  }

  while (i1 < fp1._nbits)
  {
    _bit[_nbits] = fp1._bit[i1];
    _count[_nbits] = fp1._count[i1];
    i1++;
    _nbits++;
  }

  while (i2 < fp2._nbits)
  {
    _bit[_nbits] = fp2._bit[i2];
    _count[_nbits] = - fp2._count[i2];
    i2++;
    _nbits++;
  }

// we don't set _nset, not sure what it would mean

  _nset = 0;

  return;
}

double
Sparse_Fingerprint::cosine_coefficient (const Sparse_Fingerprint & rhs) const
{
  if (0 == _sum_squared || 0 == rhs._sum_squared)
    return 0.0;

  int sum_product = 0;

  int i1 = 0;
  int i2 = 0;

#ifdef DEBUG_COSINE_COEFFICIENT
  cerr << "Cosine between " << _nbits << " bits and " << rhs._nbits << endl;
#endif

  while (i1 < _nbits && i2 < rhs._nbits)
  {
#ifdef DEBUG_COSINE_COEFFICIENT
    cerr << " cmp " << i1 << ' ' << _bit[i1] << " and " << i2 << ' ' << rhs._bit[i2] << ", counts " << _count[i1] << " and " << rhs._count[i2] << endl;
#endif
    if (_bit[i1] < rhs._bit[i2])
      i1++;
    else if (_bit[i1] == rhs._bit[i2])
    {
      sum_product += _count[i1] * _count[i1];
      i1++;
      i2++;
    }
    else
      i2++;
  }

  if (0 == sum_product)
    return 0.0;

#ifdef DEBUG_COSINE_COEFFICIENT
  cerr << " sum_product " << sum_product << " sum2 " << _sum_squared << " sum2 " << rhs._sum_squared << endl;
#endif

  return (static_cast<double>(sum_product) / sqrt(static_cast<double>(_sum_squared) * static_cast<double>(rhs._sum_squared)) + 1.0) * 0.5;
}

similarity_type_t
Sparse_Fingerprint::continuous_tanimoto (const Sparse_Fingerprint & rhs) const
{
  if (0 == _sum_squared || 0 == rhs._sum_squared)
    return 0.0;

  int i1 = 0;
  int i2 = 0;

  int sum_product = 0;

  while (i1 < _nbits && i2 < rhs._nbits)
  {
    if (_bit[i1] < rhs._bit[i2])
      i1++;
    else if (_bit[i1] == rhs._bit[i2])
    {
      sum_product += _count[i1] * rhs._count[i2];
      i1++;
      i2++;
    }
    else
      i2++;
  }

//cerr << "continuous_tanimoto_exponent " << continuous_tanimoto_exponent << " sum_product " << sum_product << " sq " << _sum_squared << " rhs " << rhs._sum_squared << endl;
//cerr << static_cast<double>(sum_product) / static_cast<double>(_sum_squared + rhs._sum_squared) << endl;

// Because we deal with counts, we will never have negative values, so we do NOT do the scaling from the -1.33 to 1 range

  if (1.0 != continuous_tanimoto_exponent)
    return pow((static_cast<double>(sum_product) / static_cast<double>(_sum_squared + rhs._sum_squared - sum_product)), continuous_tanimoto_exponent);
  else
    return (static_cast<float>(sum_product) / static_cast<float>(_sum_squared + rhs._sum_squared - sum_product));
}

int
Sparse_Fingerprint::append_daylight_ascii_form_with_counts_encoded (IWString & s) const
{
  if (0 == _nbits)
    return 1;

  int bytes_needed = _nbits * 5;   // need 5 bytes per sparse bit
  if (0 != bytes_needed % 4)
    bytes_needed = (bytes_needed / 4 + 1) * 4;

  IW_Bits_Base b(bytes_needed * IW_BITS_PER_BYTE);

  unsigned int * y = (unsigned int *)(b.bits());  // dangerous C type cast

  union 
  {
    unsigned int counts;
    unsigned char c[4];
  } counts;

  int counts_ndx = 0;   // index into union above

  int b_ndx = 0;   // index into y array

  counts.counts = 0;

  for (int i = 0; i < _nbits; i++)
  {
    y[b_ndx] = htonl(_bit[i]);
    b_ndx++;

    counts.c[counts_ndx] = static_cast<unsigned char>(_count[i]);
    counts_ndx++;

    if (4 == counts_ndx)
    {
      y[b_ndx] = counts.counts;
      b_ndx++;
      counts_ndx = 0;
    }
  }

  if (counts_ndx > 0)
  {
    y[b_ndx] = counts.counts;
    b_ndx++;
  }

  IWString tmp;
  b.daylight_ascii_representation(tmp);

  s << tmp;
  
  return 1;
}

int
Sparse_Fingerprint::count_for_bit (unsigned int b) const
{
  for (int i = 0; i < _nbits; i++)
  {
    if (b == _bit[i])
      return _count[i];
  }

  return 0;
}

class UnsignedIntComparator
{
  private:
  public:
    int operator () (unsigned int, unsigned int) const;
};

int
UnsignedIntComparator::operator () (unsigned int i1, unsigned int i2) const
{
  if (i1 < i2)
    return -1;

  if (i1 > i2)
    return 1;

  return 0;
}

int
Sparse_Fingerprint::build_from_sparse_fingerprint_creator (Sparse_Fingerprint_Creator & sfc)
{
  int n = sfc.nbits();

  if (! resize(n))
    return 0;

  if (0 == n)
    return 1;

  int notused;
  sfc.copy_bits_to_unsigned_int_array(_bit, notused);

  assert (notused == _nbits);

  UnsignedIntComparator uic;
  iwqsort(_bit, _nbits, uic);

  sfc.fill_count_array(_bit, _count, n);

  return 1;
}

#ifdef __GNUG__
template void iwqsort<unsigned, UnsignedIntComparator>(unsigned*, int, UnsignedIntComparator&);
template void iwqsort<unsigned, UnsignedIntComparator>(unsigned*, int, UnsignedIntComparator&, void*);
template void compare_two_items<unsigned, UnsignedIntComparator>(unsigned*, UnsignedIntComparator&, void*);
template void swap_elements<unsigned>(unsigned&, unsigned&, void*);
template void move_in_from_right<unsigned, UnsignedIntComparator>(unsigned*, int&, int&, UnsignedIntComparator&);
template void move_in_from_left<unsigned, UnsignedIntComparator>(unsigned*, int&, int&, int, UnsignedIntComparator&, void*);
#endif

void *
Sparse_Fingerprint::copy_to_contiguous_storage (void * p) const
{
  memcpy(p, this, sizeof(Sparse_Fingerprint));

  p = reinterpret_cast<Sparse_Fingerprint *>(p) + 1;

  memcpy(p, _bit, _nbits * sizeof(int));

  p = reinterpret_cast<int *>(p) + _nbits;

  memcpy(p, _count, _nbits * sizeof(int));

  p = reinterpret_cast<int *>(p) + _nbits;

  return p;
}

void *
Sparse_Fingerprint::copy_to_contiguous_storage_gpu (void * p) const
{
  memcpy(p, &_nbits, sizeof(int));

  p = reinterpret_cast<int *>(p) + 1;

  memcpy(p, &_nset, sizeof(int));

  p = reinterpret_cast<int *>(p) + 1;

  memcpy(p, _bit, _nbits * sizeof(int));

  p = reinterpret_cast<int *>(p) + _nbits;
 
  unsigned int terminate=std::numeric_limits<unsigned int>::max();
    
  memcpy(p, &terminate, sizeof(unsigned int));

  p = reinterpret_cast<unsigned int *>(p) + 1;

  memcpy(p, _count, _nbits * sizeof(int));

  p = reinterpret_cast<int *>(p) + _nbits;

  return p;
}

const void *
Sparse_Fingerprint::build_from_contiguous_storage (const void * p,
                                        int allocate_arrays)
{
  if (allocate_arrays && NULL != _bit)
  {
    delete [] _bit;
    delete [] _count;
  }

  memcpy(this, p, sizeof(Sparse_Fingerprint));

  p = reinterpret_cast<const Sparse_Fingerprint *>(p) + 1;

  if (allocate_arrays)
  {
    resize(_nbits);

    memcpy(_bit, p, _nbits * sizeof(int));
  }
  else
    _bit = (unsigned int *) p;

  p = reinterpret_cast<const int *>(p) + _nbits;

  if (allocate_arrays)
    memcpy(_count, p, _nbits * sizeof(int));
  else
    _count = (int *) p;

  p = reinterpret_cast<const int *>(p) + _nbits;

  return p;
}

int
Sparse_Fingerprint::remove_bit (unsigned int b)
{
  auto f = std::lower_bound(_bit, _bit + _nbits, b);

  if (f == (_bit + _nbits) || *f != b)    // bit not present to be removed
    return 0;

  _nset -= _count[f - _bit];

  _nbits--;

  for (auto i = f - _bit; i < _nbits; ++i)
  {
    _bit[i] = _bit[i+1];
    _count[i] = _count[i+1];
  }


  return 1;
}

int
Sparse_Fingerprint::set_count (unsigned int b, int c)
{
  auto f = std::lower_bound(_bit, _bit + _nbits, b);

  if (f == (_bit + _nbits) || *f != b)    // bit not present to be removed
    return 0;

  const auto ndx = f - _bit;

  _count[ndx] = c;

  return 1;
}
