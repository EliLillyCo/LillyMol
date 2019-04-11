#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwstring.h"

#include "dyfp.h"
#include "tversky.h"
#include "various_distance_metrics.h"

void
IWDYFP::_default_values()
{
  _nset   = IWDYFP_NSET_NOT_COMPUTED;
  _whole_words = 0;
}

IWDYFP::IWDYFP()
{
  _default_values();

  return;
}

IWDYFP::IWDYFP(int initial_nbits) : IW_Bits_Base(initial_nbits)
{
  _default_values();

  assert (0 == (0x00000007 & _nbits));    // must be a multiple of 8 bits

  _determine_whole_words();

  return;
}

IWDYFP &
IWDYFP::operator = (const IWDYFP & rhs)
{
  IW_Bits_Base::operator = (rhs);

  _nset = rhs._nset;
  _whole_words = rhs._whole_words;

  return *this;
}

IWDYFP::~IWDYFP()
{
  if (-807 == _nset)
    cerr << "Deleting an already deleted IWDYFP\n";
  _nset = -807;
}

int
IWDYFP::ok() const
{
  if (! IW_Bits_Base::ok())
    return 0;

  if (IWDYFP_NSET_NOT_COMPUTED == _nset && is_empty())
    return 1;

// At this stage, at least some things are non zero.

  if (_nset >= 0)
    ;
  else if (IWDYFP_NSET_NOT_COMPUTED == _nset)
    ;
  else
    return 0;

  return 1;
}

int
IWDYFP::debug_print(std::ostream & os) const
{
  os << "Details on IWDYFP, nset";
  if (IWDYFP_NSET_NOT_COMPUTED == _nset)
    os << " not computed ";
  else
    os << '=' << _nset;

  return IW_Bits_Base::debug_print(os);
}

int
IWDYFP::allocate_space_for_bits(int n)
{
  int rc = IW_Bits_Base::allocate_space_for_bits(n);

  if (0 == rc)
    return 0;

  _determine_whole_words();

  return rc;
}

/*
  If the last whole word is not all being used, we need to mask out
  any unused area.

  Sept 04. I'm not 100% sure this is still needed.
*/

static unsigned int * mask = NULL;

/*
  A little complicated because we need to do this in an endian
  independent fashion
*/

static void
establish_mask(int xtra,
                unsigned int & mask)
{
  unsigned char * c = reinterpret_cast<unsigned char *>(&mask);
  c[0] = one_bit_8[0];

  for (int i = 0; i < xtra; i++)
  {
    int j = i / 8;    // which byte are we doing

    c[j] = (c[j] | one_bit_8[i % 8]);

//  cerr << "Mask updated to " << mask << endl;

    mask = mask | (mask >> 1);
  }

//cerr << "for " << xtra << " extra bits, mask set to " << mask << endl;

  return;
}

static void
establish_mask()
{
  mask = new unsigned int[IW_BITS_PER_WORD];

  for (int i = 0; i < IW_BITS_PER_WORD; i++)
  {
    establish_mask(i, mask[i]);
  }

  return;
}

void
IWDYFP::_determine_whole_words()
{
  int xtra = _nbits % IW_BITS_PER_WORD;

  if (0 == xtra)
    _whole_words = _nbits / IW_BITS_PER_WORD;
  else
  {
    _whole_words = _nbits / IW_BITS_PER_WORD + 1;

    if (NULL == mask)
      establish_mask();

    unsigned int * tmp = reinterpret_cast<unsigned int *>(_bits) + _whole_words - 1;

    *tmp = ( (*tmp) & mask[xtra]);
  }

  return;
}

/*
  Ideally the number of bits set would be <= nbits(), but since
  we may have been created from a sparse fingerprint, we could have
  the situation where _nset is > _nbits
*/

int
IWDYFP::set_nset(int ns)
{
  assert (ns >= 0);

  _nset = ns;

  return 1;
}

int
IWDYFP::construct_from_tdt_record(const IWString & tdt_record)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (! IW_Bits_Base::construct_from_tdt_record_nset(tdt_record, _nset))
    return 0;

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_tdt_record(const const_IWSubstring & tdt_record)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (! IW_Bits_Base::construct_from_tdt_record_nset(tdt_record, _nset))
    return 0;

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_daylight_ascii_bit_rep(const const_IWSubstring & s)
{
  if (! IW_Bits_Base::construct_from_daylight_ascii_bit_rep(s.rawchars(), s.length()))
  {
    cerr << "IWDYFP::construct_from_daylight_ascii_bit_rep:cannot parse Daylight ascii form\n";
    return 0;
  }

  nset();

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_daylight_ascii_representation(const_IWSubstring const & buffer)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;

  const_IWSubstring mybuffer(buffer);

  int semicolon = buffer.find(';');

  if (semicolon < 0)
  {
    cerr << "IWDYFP::construct_from_daylight_ascii_representation: must contain semicolon\n";
    return 0;
  }

  if (! IW_Bits_Base::construct_from_daylight_ascii_bit_rep(buffer.rawchars(), semicolon))
    return 0;

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_sparse_representation(const const_IWSubstring & buffer)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (! IW_Bits_Base::construct_from_sparse_representation(buffer))
    return 0;

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_ascii_01_representation(const char * buffer, int nchars)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;

  int rc;
  if (' ' == buffer[1])    // didn't check that n > 0
    rc = IW_Bits_Base::construct_from_ascii_01_representation_with_spaces(buffer, nchars);
  else
    rc = IW_Bits_Base::construct_from_ascii_01_representation(buffer, nchars);

  if (0 == rc)
    return 0;

  nset();

  _determine_whole_words();

//cerr << "Fingerprint contains " << _nbits << " bits\n";

  return 1;
}

/*
  This is wrong if we have extra bits in _nset
*/

int
IWDYFP::write_daylight_ascii_representation(std::ostream & os,
                                const const_IWSubstring & dataitem_name)
{
  (void) nset();    // force computation of _nset if needed

  IWString ascii;
  if (! IW_Bits_Base::daylight_ascii_representation(ascii) || 0 == ascii.length())
  {
    cerr << "IWDYFP::write_daylight_ascii_representation: failure\n";
    return 0;
  }

  os << dataitem_name;
  if (! dataitem_name.ends_with('<'))
    os << '<';
  os << ascii << ';' << _nbits << ';' << _nset << ';' <<
        _nbits << ';' << _nset << ";1>\n";

  return os.good();
}

int
IWDYFP::daylight_ascii_representation(IWString & result)
{
  (void) nset();    // force computation of _nset if needed

  if (! IW_Bits_Base::daylight_ascii_representation(result) || 0 == result.length())
  {
    cerr << "IW_Bits_Base::daylight_ascii_representation: failure\n";
    return 0;
  }

  result.resize(result.length() + 25);

  result << ';' << _nbits << ';' << _nset << ';' << _nbits << ';' << _nset << ";1";

  return 1;
}

int
IWDYFP::daylight_ascii_tdt(IWString & result, const const_IWSubstring & tag)
{
  result.resize(_nbits / 2);    // significant overestimate

  result = tag;

  if (! result.ends_with('<'))
    result += '<';

  (void) nset();

  IWString ascii;
  (void) IW_Bits_Base::daylight_ascii_representation(ascii);

  result += ascii;

  result << ';' << _nbits << ';' << _nset << ';' << _nbits << ';' << _nset << ";1>";

  return 1;
}

/*
  And two fingerprints, and return the number of bits in common.
  The fingerprints are processed word at a time. Note that this is
  highly dependent on 4 bytes per unsigned int.
*/

extern "C" int bits_in_common(const unsigned int *, const unsigned int *, int);
extern "C" int bits_in_common8(const unsigned int *, const unsigned int *, int);

int
IWDYFP::_bits_in_common(const IWDYFP & f2) const
{
  if (f2._whole_bytes != _whole_bytes)
    cerr << "Dying, " << _whole_bytes << " vs " << f2._whole_bytes << endl;

  assert (f2._whole_bytes == _whole_bytes);

  return ::bits_in_common((const unsigned int *) _bits,
                          (const unsigned int *) f2._bits,
                           _whole_words);
}

/*
  Somewhat special purpose function to and another fingerprint with
  this one, and indicate whether or not the operation changed THIS
*/

void
IWDYFP::iwand(const IWDYFP & f2, int & changed) 
{
  IW_Bits_Base::iwand(f2, changed);

  if (changed)
    _nset = IWDYFP_NSET_NOT_COMPUTED;

  return;
}

void
IWDYFP::iwor(const IWDYFP & f2) 
{
  IW_Bits_Base::iwor(f2);

  _nset = IWDYFP_NSET_NOT_COMPUTED;

  return;
}

/*
  Somewhat special purpose function to OR another fingerprint with
  this one, and indicate whether or not the operation changed THIS
*/

void
IWDYFP::iwor(const IWDYFP & f2, int & changed) 
{
  IW_Bits_Base::iwor(f2, changed);

  if (changed)
    _nset = IWDYFP_NSET_NOT_COMPUTED;

  return;
}

/*
*/

void
IWDYFP::iwxor(const IWDYFP & f2)
{
  IW_Bits_Base::iwxor(f2);

  _nset = IWDYFP_NSET_NOT_COMPUTED;

  return;
}

/*
  f2 must contain >= the number of bits
*/

//#define DEBUG_IMPLEMENTATION_TANIMOTO

similarity_type_t
IWDYFP::_tanimoto_multiplier(IWDYFP & f2)
{
  int multiplier = 1;
  if (f2._whole_bytes > _whole_bytes)
     multiplier = f2._whole_bytes / _whole_bytes;

  int nc = _bits_in_common(f2);

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
  cerr << "nb " << _nbits << " nset " << _nset << " and " << f2._nset << " nc = " << nc << " _whole_words " << _whole_words << ", bytes " << _whole_bytes << endl;
  cerr << "computed nset " << compute_nset() << " and " << f2.compute_nset() << endl;
#endif

  if (0 == nc)
  {
    if (0 == _nset && 0 == f2._nset)      // otherwise we have identical molecules that will have non-zero distances
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
#endif

  return similarity_type_t(nc) /
         similarity_type_t(multiplier * _nset + f2._nset - nc);
}

similarity_type_t
IWDYFP::_tanimoto(IWDYFP & f2)
{
  assert (_whole_bytes == f2._whole_bytes);

  const int nc = _bits_in_common(f2);

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
  cerr << "nb " << _nbits << " nset " << _nset << " and " << f2._nset << " nc = " << nc << " _whole_words " << _whole_words << ", bytes " << _whole_bytes << endl;
  cerr << "computed nset " << compute_nset() << " and " << f2.compute_nset() << endl;
#endif

  if (0 == nc)
  {
    if (0 == _nset && 0 == f2._nset)      // otherwise we have identical molecules that will have non-zero distances
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
#endif

  return similarity_type_t(nc) /
         similarity_type_t(_nset + f2._nset - nc);
}

similarity_type_t
IWDYFP::tanimoto(IWDYFP & f2)
{
  if (_whole_bytes == f2._whole_bytes)    // the most common case
    return _tanimoto(f2);

  if (_whole_bytes > f2._whole_bytes)
    return _tanimoto_multiplier(f2);
  else
    return f2._tanimoto_multiplier(*this);
}

/*
  Tversky similarity 
*/

similarity_type_t
IWDYFP::tversky(IWDYFP & rhs,
                 const Tversky & tv)
{
  assert (_whole_bytes == rhs._whole_bytes);

  if (_nset < 0)
    (void) nset();

  if (rhs._nset < 0)
    (void) rhs.nset();

  if (0 == _nset || 0 == rhs._nset)
  {
    if (0 == _nset && 0 == rhs._nset)
      return static_cast<similarity_type_t>(1.0);

    if (! tv.nset_sensitive_zero_bit_similarity())
      return static_cast<similarity_type_t>(0.0);

    if (_nset)
      return static_cast<similarity_type_t>(1.0) / static_cast<similarity_type_t>(_nset + 1);
    else
      return static_cast<similarity_type_t>(1.0) / static_cast<similarity_type_t>(rhs._nset + 1);
  }

  int c = ::bits_in_common((const unsigned int *) _bits,
                           (const unsigned int *) rhs._bits,
                            _whole_words);

// Use Bradshaw's notation

  int a = _nset - c;
  int b = rhs._nset - c;

//cerr << "a = " << a << " b = " << b << " c = " << c << endl;
  return (float(c) / (tv.a() * float(a) + tv.b() * float(b) + float(c)));
}

/*
  The idea of a distance metric where you compute the Tversky
  both ways, and the Tanimoto, and then take the shortest distance
*/

similarity_type_t
IWDYFP::optimistic_distance(IWDYFP & f2,
                             const Tversky & tv)
{
  assert (_whole_bytes == f2._whole_bytes);

  if (_nset < 0)
    (void) nset();

  if (f2._nset < 0)
    (void) f2.nset();

  if (0 == _nset || 0 == f2._nset)
  {
    if (_nset == f2._nset)     // no bits set in either
      return static_cast<similarity_type_t>(0.0);

    if (! tv.nset_sensitive_zero_bit_similarity())
      return static_cast<similarity_type_t>(1.0);

    if (_nset)
      return static_cast<similarity_type_t>(_nset) / static_cast<similarity_type_t>(_nset + 1);
    else
      return static_cast<similarity_type_t>(f2._nset) / static_cast<similarity_type_t>(f2._nset + 1);
  }

  int c = ::bits_in_common((const unsigned int *) _bits,
                           (const unsigned int *) f2._bits,
                            _whole_words);

  if (0 == c)
    return static_cast<similarity_type_t>(1.0);

// Use Bradshaw's notation

  float a = static_cast<float>(_nset - c);
  float b = static_cast<float>(f2._nset - c);

//cerr << "a = " << a << " b = " << b << " c = " << c << endl;

  similarity_type_t tv1 = float(c) / (tv.a() * a + tv.b() * b + float(c));
  similarity_type_t tv2 = float(c) / (tv.b() * a + tv.a() * b + float(c));
  similarity_type_t tnm = float(c) / static_cast<float>(_nset + f2._nset - c);

// These are similarities right now, so return 1.0 - the largest similarity

//cerr << "Tv1 " << tv1 << ", tv2 " << tv2 << " tn " << tnm << endl;
  if (tv1 >= tv2 && tv1 >= tnm)
    return static_cast<similarity_type_t>(1.0) - tv1;

  if (tv2 >= tv1 && tv2 >= tnm)
    return static_cast<similarity_type_t>(1.0) - tv2;

  return static_cast<similarity_type_t>(1.0) - tnm;
}

/*
  Sometimes it is more convenient to specify the tversky coefficients
  up front, and then call the function with no arguments
*/

static tversky_coeff_t global_tversky_alpha = 1.0;
static tversky_coeff_t global_tversky_beta = 1.0;

int
set_global_tversky_alpha(tversky_coeff_t na)
{
  global_tversky_alpha = na;

  return 1;
}

int
set_global_tversky_beta(tversky_coeff_t nb)
{
  global_tversky_beta = nb;

  return 1;
}

int
set_global_tversky_coefficients(tversky_coeff_t na, tversky_coeff_t nb)
{
  global_tversky_alpha = na;
  global_tversky_beta  = nb;

  return 1;
}

similarity_type_t
IWDYFP::tversky(IWDYFP & f2)
{
  assert (_whole_bytes == f2._whole_bytes);

  if (_nset < 0)
    (void) nset();

  if (f2._nset < 0)
    (void) f2.nset();

  if (0 == _nset || 0 == f2._nset)
  {
    if (0 == _nset && 0 == f2._nset)
      return static_cast<similarity_type_t>(1.0);

    return static_cast<similarity_type_t>(0.0);
  }

  int c = ::bits_in_common((const unsigned int *) _bits,
                           (const unsigned int *) f2._bits,
                            _whole_words);

// Use Bradshaw's notation

  int a = _nset - c;
  int b = f2._nset - c;

  return (float(c) / (global_tversky_alpha * float(a) +
                        global_tversky_beta  * float(b) + float(c)));
}

similarity_type_t
IWDYFP::fraction_matched(IWDYFP & f2)
{
  (void) nbits();

  int nc;
  if (f2.nbits() > _nbits)
    nc = _bits_in_common(f2);
  else
    nc = f2._bits_in_common(*this);

  return similarity_type_t(nc) / similarity_type_t(_nbits);
}

similarity_type_t
tanimoto(IWDYFP & fp1, IWDYFP & fp2)
{
  return fp1.tanimoto(fp2);
}

similarity_type_t
fraction_matched(IWDYFP & fp1, IWDYFP & fp2)
{
  return fp1.fraction_matched(fp2);
}

similarity_type_t
tversky(IWDYFP & fp1, IWDYFP & fp2)
{
  return fp1.tversky(fp2);
}

/*
  The only complication is the need to maintain the _nset value
*/

void
IWDYFP::operator += (const IWDYFP & rhs)
{
  IW_Bits_Base::operator += (rhs);

  if (IWDYFP_NSET_NOT_COMPUTED == _nset)
    ;
  else if (IWDYFP_NSET_NOT_COMPUTED == rhs._nset)
    _nset = IWDYFP_NSET_NOT_COMPUTED;
  else
    _nset += rhs._nset;

  return;
}

int
IWDYFP::fold(int nfold)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;

  _whole_words = _whole_words / 2;

  return IW_Bits_Base::fold(nfold);
}

/*
  We overload this method so we can compute _nset at the same time
*/

int
IWDYFP::construct_from_array_of_ints(const int * ii,
                                            int nb)
{
  assert (nb > 0);

  if (_nbits)      // zero out any pre-existing information
    clear();

  allocate_space_for_bits(nb);
  _determine_whole_words();

  _nset = 0;

  for (int i = 0; i < _whole_bytes; i++)
  {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++)
    {
      if (*ii)
      {
        _bits[i] |= one_bit_8[j];
        _nset++;
      }
      ii++;
    }
  }

  if (_extra_bits)
  {
    for (int i = 0; i < _extra_bits; i++)
    {
      if (*ii)
      {
        _bits[_whole_bytes - 1] |= one_bit_8[i];
        _nset++;
      }
      ii++;
        
    }
  }

  return 1;
}

int
IWDYFP::construct_from_descriptor_tdt_record(const IWString & buffer)
{
  assert (buffer.ends_with('>'));

  const_IWSubstring tmp = buffer;

  tmp.remove_up_to_first('<');
  tmp.chop();

  return construct_from_descriptor_record(tmp);
}

int
/*
  Sometimes we may have descriptors that are to be interpreted as bits.
  Anything that isn't '0' is considered as an on bit
*/

IWDYFP::construct_from_descriptor_record(const const_IWSubstring & buffer)
{
  if (_nbits)
    clear();

  _nset = 0;

  int nb = buffer.nwords();

  allocate_space_for_bits(nb);
  _determine_whole_words();

  int j = 0;
  const_IWSubstring token;
  for (int i = 0; i < _nbits; i++)
  {
    if (! buffer.nextword(token, j))
    {
      cerr << "IW_DY_Fingerprint::construct_from_descriptor_record: not enough tokens\n";
      cerr << buffer << endl;
      return 0;
    }

    if ('0' != token)
    {
      set(i);
      _nset++;
    }
  }

  return 1;
}

int
IWDYFP::construct_from_hex(const const_IWSubstring & zhex)
{
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (! IW_Bits_Base::construct_from_hex(zhex))
    return 0;
 
  _determine_whole_words();

  nset();

  return 1;
}

similarity_type_t
IWDYFP::fvb_modified_tanimoto(IWDYFP & rhs)
{
  int ns1 = nset();
  int ns2 = rhs.nset();

  int n11;    // bits set in both

  if (0 == ns1 || 0 == ns2)
    n11 = 0;
  else
    n11 = _bits_in_common(rhs);

  int n00 = _nbits - ns1 - ns2 + n11;    // bits set in neither

  return fligner_verducci_blower(_nbits, _nset, rhs._nset, n00, n11);
}

similarity_type_t
IWDYFP::russel_rao(IWDYFP & rhs)
{
  if (0 == nset() || 0 == rhs.nset())
    return static_cast<similarity_type_t>(0.0);

  int n11 = _bits_in_common(rhs);

  return static_cast<similarity_type_t>(n11) / static_cast<similarity_type_t>(_nbits);
}

similarity_type_t
IWDYFP::forbes_similarity(IWDYFP & rhs)
{
  int ns1 = nset();
  int ns2 = rhs.nset();

  if (0 == ns1 || 0 == ns2)
    return static_cast<similarity_type_t>(0.0);

  int a = _bits_in_common(rhs);
  if (0 == a)
    return static_cast<similarity_type_t>(0.0);

//int b = ns1 - a;    // just set in just LHS
//int c = ns2 - a;    // just set in just RHS

  return static_cast<similarity_type_t>(a) / static_cast<similarity_type_t>(ns1 * ns2);
}

similarity_type_t
IWDYFP::simple_matching(IWDYFP & rhs)
{
  int ns1 = nset();
  int ns2 = rhs.nset();

  int a = _bits_in_common(rhs);

  int d = _nbits - ns1 - ns2 + a;    // bits set in neither

  return static_cast<similarity_type_t>(a + d) / static_cast<similarity_type_t>(_nbits);
}

similarity_type_t
IWDYFP::sorensendice(IWDYFP & rhs)
{
  const int bic = _bits_in_common(rhs);

  return 1.0f - static_cast<similarity_type_t>(bic) / static_cast<similarity_type_t>(nset() + rhs.nset());
}

similarity_type_t
IWDYFP::overlap(IWDYFP & rhs)
{
  const int bic = _bits_in_common(rhs);
  if (0 == bic)
    return static_cast<similarity_type_t>(1.0);

  const int ns1 = nset();
  const int ns2 = rhs.nset();

  if (ns1 < ns2)
    return static_cast<similarity_type_t>(1.0) - static_cast<similarity_type_t>(bic) / static_cast<similarity_type_t>(ns1);
  else
    return static_cast<similarity_type_t>(1.0) - static_cast<similarity_type_t>(bic) / static_cast<similarity_type_t>(ns2);
}

void *
IWDYFP::copy_to_contiguous_storage(void * p) const
{
//const unsigned char * initp = reinterpret_cast<const unsigned char *>(p);
  p = IW_Bits_Base::copy_to_contiguous_storage(p);

  size_t sbase = sizeof(IW_Bits_Base);

  memcpy(p, reinterpret_cast<const unsigned char *>(this) + sbase, sizeof(IWDYFP) - sbase);

  p = reinterpret_cast<unsigned char *>(p) + sizeof(IWDYFP) - sbase;

  return p;
}

void *
IWDYFP::copy_to_contiguous_storage_gpu(void * p) const
{
//const unsigned char * initp = reinterpret_cast<const unsigned char *>(p);
  p = IW_Bits_Base::copy_to_contiguous_storage_gpu(p);


  memcpy(p, &_nset, sizeof(int));

  p = reinterpret_cast<unsigned char *>(p) + sizeof(int);

  return p;
}

const void *
IWDYFP::build_from_contiguous_storage(const void * p,
                                       int allocate_arrays)
{
  p = IW_Bits_Base::build_from_contiguous_storage(p, allocate_arrays);

  size_t sbase = sizeof(IW_Bits_Base);

  memcpy(reinterpret_cast<unsigned char *>(this) + sbase, p, sizeof(IWDYFP) - sbase);

  return reinterpret_cast<const unsigned char *>(p) + sizeof(IWDYFP) - sbase;
}
