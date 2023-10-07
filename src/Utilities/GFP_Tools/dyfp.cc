#include "dyfp.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Foundational/iwstring/iwstring.h"
#include "tversky.h"
#include "various_distance_metrics.h"

using std::cerr;

float
Fraction(int numerator, int denomiator) {
  return static_cast<float>(numerator) / static_cast<float>(denomiator);
}

void
IWDYFP::_default_values() {
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  _whole_words = 0;
}

IWDYFP::IWDYFP() {
  _default_values();

  return;
}

IWDYFP::IWDYFP(int initial_nbits) : IW_Bits_Base(initial_nbits) {
  _default_values();

  assert(0 == (0x00000007 & _nbits));  // must be a multiple of 8 bits

  _determine_whole_words();

  return;
}

IWDYFP &
IWDYFP::operator=(const IWDYFP &rhs) {
  IW_Bits_Base::operator=(rhs);

  _nset = rhs._nset;
  _whole_words = rhs._whole_words;

  return *this;
}

IWDYFP::~IWDYFP() {
  if (-807 == _nset) {
    cerr << "Deleting an already deleted IWDYFP\n";
  }
  _nset = -807;
}

int
IWDYFP::ok() const {
  if (!IW_Bits_Base::ok()) {
    return 0;
  }

  if (IWDYFP_NSET_NOT_COMPUTED == _nset && is_empty()) {
    return 1;
  }

  // At this stage, at least some things are non zero.

  if (_nset >= 0) {
    ;
  } else if (IWDYFP_NSET_NOT_COMPUTED == _nset) {
    ;
  } else {
    return 0;
  }

  return 1;
}

int
IWDYFP::debug_print(std::ostream &os) const {
  os << "Details on IWDYFP, nset";
  if (IWDYFP_NSET_NOT_COMPUTED == _nset) {
    os << " not computed ";
  } else {
    os << '=' << _nset;
  }

  return IW_Bits_Base::debug_print(os);
}

int
IWDYFP::allocate_space_for_bits(int n) {
  int rc = IW_Bits_Base::allocate_space_for_bits(n);

  if (0 == rc) {
    return 0;
  }

  _determine_whole_words();

  return rc;
}

/*
  If the last whole word is not all being used, we need to mask out
  any unused area.

  Sept 04. I'm not 100% sure this is still needed.
*/

static unsigned int *mask = NULL;

/*
  A little complicated because we need to do this in an endian
  independent fashion
*/

static void
establish_mask(int xtra, unsigned int &mask) {
  unsigned char *c = reinterpret_cast<unsigned char *>(&mask);
  c[0] = one_bit_8[0];

  for (int i = 0; i < xtra; i++) {
    int j = i / 8;  // which byte are we doing

    c[j] = (c[j] | one_bit_8[i % 8]);

    //  cerr << "Mask updated to " << mask << '\n';

    mask = mask | (mask >> 1);
  }

  // cerr << "for " << xtra << " extra bits, mask set to " << mask << '\n';

  return;
}

static void
establish_mask() {
  mask = new unsigned int[IW_BITS_PER_WORD];

  for (int i = 0; i < IW_BITS_PER_WORD; i++) {
    establish_mask(i, mask[i]);
  }

  return;
}

void
IWDYFP::_determine_whole_words() {
  int xtra = _nbits % IW_BITS_PER_WORD;

  if (0 == xtra) {
    _whole_words = _nbits / IW_BITS_PER_WORD;
  } else {
    _whole_words = _nbits / IW_BITS_PER_WORD + 1;

    if (NULL == mask) {
      establish_mask();
    }

    unsigned int *tmp = reinterpret_cast<unsigned int *>(_bits) + _whole_words - 1;

    *tmp = ((*tmp) & mask[xtra]);
  }

  return;
}

/*
  Ideally the number of bits set would be <= nbits(), but since
  we may have been created from a sparse fingerprint, we could have
  the situation where _nset is > _nbits
*/

int
IWDYFP::set_nset(int ns) {
  assert(ns >= 0);

  _nset = ns;

  return 1;
}

int
IWDYFP::construct_from_tdt_record(const IWString &tdt_record) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (!IW_Bits_Base::construct_from_tdt_record_nset(tdt_record, _nset)) {
    return 0;
  }

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_tdt_record(const const_IWSubstring &tdt_record) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (!IW_Bits_Base::construct_from_tdt_record_nset(tdt_record, _nset)) {
    return 0;
  }

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_daylight_ascii_bit_rep(const const_IWSubstring &s) {
  if (!IW_Bits_Base::construct_from_daylight_ascii_bit_rep(s.rawchars(), s.length())) {
    cerr << "IWDYFP::construct_from_daylight_ascii_bit_rep:cannot parse Daylight ascii "
            "form\n";
    return 0;
  }

  nset();

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_daylight_ascii_representation(const_IWSubstring const &buffer) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;

  const_IWSubstring mybuffer(buffer);

  int semicolon = buffer.find(';');

  if (semicolon < 0) {
    cerr << "IWDYFP::construct_from_daylight_ascii_representation: must contain "
            "semicolon\n";
    return 0;
  }

  if (!IW_Bits_Base::construct_from_daylight_ascii_bit_rep(buffer.rawchars(),
                                                           semicolon)) {
    return 0;
  }

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_sparse_representation(const const_IWSubstring &buffer) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (!IW_Bits_Base::construct_from_sparse_representation(buffer)) {
    return 0;
  }

  _determine_whole_words();

  return 1;
}

int
IWDYFP::construct_from_ascii_01_representation(const char *buffer, int nchars) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;

  int rc;
  if (' ' == buffer[1]) {  // didn't check that n > 0
    rc = IW_Bits_Base::construct_from_ascii_01_representation_with_spaces(buffer, nchars);
  } else {
    rc = IW_Bits_Base::construct_from_ascii_01_representation(buffer, nchars);
  }

  if (0 == rc) {
    return 0;
  }

  nset();

  _determine_whole_words();

  // cerr << "Fingerprint contains " << _nbits << " bits\n";

  return 1;
}

/*
  This is wrong if we have extra bits in _nset
*/

int
IWDYFP::write_daylight_ascii_representation(std::ostream &os,
                                            const const_IWSubstring &dataitem_name) {
  (void)nset();  // force computation of _nset if needed

  IWString ascii;
  if (!IW_Bits_Base::daylight_ascii_representation(ascii) || 0 == ascii.length()) {
    cerr << "IWDYFP::write_daylight_ascii_representation: failure\n";
    return 0;
  }

  os << dataitem_name;
  if (!dataitem_name.ends_with('<')) {
    os << '<';
  }
  os << ascii << ';' << _nbits << ';' << _nset << ';' << _nbits << ';' << _nset
     << ";1>\n";

  return os.good();
}

int
IWDYFP::daylight_ascii_representation(IWString &result) {
  (void)nset();  // force computation of _nset if needed

  if (!IW_Bits_Base::daylight_ascii_representation(result) || 0 == result.length()) {
    cerr << "IW_Bits_Base::daylight_ascii_representation: failure\n";
    return 0;
  }

  result.resize(result.length() + 25);

  result << ';' << _nbits << ';' << _nset << ';' << _nbits << ';' << _nset << ";1";

  return 1;
}

int
IWDYFP::daylight_ascii_tdt(IWString &result, const const_IWSubstring &tag) {
  result.resize(_nbits / 2);  // significant overestimate

  result = tag;

  if (!result.ends_with('<')) {
    result += '<';
  }

  (void)nset();

  IWString ascii;
  (void)IW_Bits_Base::daylight_ascii_representation(ascii);

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
IWDYFP::_bits_in_common(const IWDYFP &rhs) const {
  if (rhs._whole_bytes != _whole_bytes) {
    cerr << "Dying, " << _whole_bytes << " vs " << rhs._whole_bytes << '\n';
  }

  assert(rhs._whole_bytes == _whole_bytes);

  return ::bits_in_common((const unsigned int *)_bits, (const unsigned int *)rhs._bits,
                          _whole_words);
}

/*
  Somewhat special purpose function to and another fingerprint with
  this one, and indicate whether or not the operation changed THIS
*/

void
IWDYFP::iwand(const IWDYFP &rhs, int &changed) {
  IW_Bits_Base::iwand(rhs, changed);

  if (changed) {
    _nset = IWDYFP_NSET_NOT_COMPUTED;
  }

  return;
}

void
IWDYFP::iwor(const IWDYFP &rhs) {
  IW_Bits_Base::iwor(rhs);

  _nset = IWDYFP_NSET_NOT_COMPUTED;

  return;
}

/*
  Somewhat special purpose function to OR another fingerprint with
  this one, and indicate whether or not the operation changed THIS
*/

void
IWDYFP::iwor(const IWDYFP &rhs, int &changed) {
  IW_Bits_Base::iwor(rhs, changed);

  if (changed) {
    _nset = IWDYFP_NSET_NOT_COMPUTED;
  }

  return;
}

/*
 */

void
IWDYFP::iwxor(const IWDYFP &rhs) {
  IW_Bits_Base::iwxor(rhs);

  _nset = IWDYFP_NSET_NOT_COMPUTED;

  return;
}

/*
  rhs must contain >= the number of bits
*/

// #define DEBUG_IMPLEMENTATION_TANIMOTO

similarity_type_t
IWDYFP::_tanimoto_multiplier(IWDYFP &rhs) {
  int multiplier = 1;
  if (rhs._whole_bytes > _whole_bytes) {
    multiplier = rhs._whole_bytes / _whole_bytes;
  }

  int nc = _bits_in_common(rhs);

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
  cerr << "nb " << _nbits << " nset " << _nset << " and " << rhs._nset << " nc = " << nc
       << " _whole_words " << _whole_words << ", bytes " << _whole_bytes << '\n';
  cerr << "computed nset " << compute_nset() << " and " << rhs.compute_nset() << '\n';
#endif

  if (0 == nc) {
    if (0 == _nset && 0 == rhs._nset) {  // otherwise we have identical molecules that
                                         // will have non-zero distances
      return static_cast<similarity_type_t>(1.0);
    }

    return static_cast<similarity_type_t>(0.0);
  }

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
#endif

  return similarity_type_t(nc) / similarity_type_t(multiplier * _nset + rhs._nset - nc);
}

similarity_type_t
IWDYFP::_tanimoto(IWDYFP &rhs) {
  assert(_whole_bytes == rhs._whole_bytes);

  const int nc = _bits_in_common(rhs);

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
  cerr << "nb " << _nbits << " nset " << _nset << " and " << rhs._nset << " nc = " << nc
       << " _whole_words " << _whole_words << ", bytes " << _whole_bytes << '\n';
  cerr << "computed nset " << compute_nset() << " and " << rhs.compute_nset() << '\n';
#endif

  if (0 == nc) {
    if (0 == _nset && 0 == rhs._nset) {  // otherwise we have identical molecules that
                                         // will have non-zero distances
      return static_cast<similarity_type_t>(1.0);
    }

    return static_cast<similarity_type_t>(0.0);
  }

#ifdef DEBUG_IMPLEMENTATION_TANIMOTO
#endif

  return similarity_type_t(nc) / similarity_type_t(_nset + rhs._nset - nc);
}

similarity_type_t
IWDYFP::tanimoto(IWDYFP &rhs) {
  if (_whole_bytes == rhs._whole_bytes) {  // the most common case
    return _tanimoto(rhs);
  }

  if (_whole_bytes > rhs._whole_bytes) {
    return _tanimoto_multiplier(rhs);
  } else {
    return rhs._tanimoto_multiplier(*this);
  }
}

/*
  Tversky similarity
*/

similarity_type_t
IWDYFP::tversky(IWDYFP &rhs, const Tversky &tv) {
  assert(_whole_bytes == rhs._whole_bytes);

  if (_nset < 0) {
    (void)nset();
  }

  if (rhs._nset < 0) {
    (void)rhs.nset();
  }

  if (0 == _nset || 0 == rhs._nset) {
    if (0 == _nset && 0 == rhs._nset) {
      return static_cast<similarity_type_t>(1.0);
    }

    if (!tv.nset_sensitive_zero_bit_similarity()) {
      return static_cast<similarity_type_t>(0.0);
    }

    if (_nset) {
      return static_cast<similarity_type_t>(1.0) /
             static_cast<similarity_type_t>(_nset + 1);
    } else {
      return static_cast<similarity_type_t>(1.0) /
             static_cast<similarity_type_t>(rhs._nset + 1);
    }
  }

  int c = ::bits_in_common((const unsigned int *)_bits, (const unsigned int *)rhs._bits,
                           _whole_words);

  // Use Bradshaw's notation

  int a = _nset - c;
  int b = rhs._nset - c;

  // cerr << "a = " << a << " b = " << b << " c = " << c << '\n';
  return (float(c) / (tv.a() * float(a) + tv.b() * float(b) + float(c)));
}

/*
  The idea of a distance metric where you compute the Tversky
  both ways, and the Tanimoto, and then take the shortest distance
*/

similarity_type_t
IWDYFP::optimistic_distance(IWDYFP &rhs, const Tversky &tv) {
  assert(_whole_bytes == rhs._whole_bytes);

  if (_nset < 0) {
    (void)nset();
  }

  if (rhs._nset < 0) {
    (void)rhs.nset();
  }

  if (0 == _nset || 0 == rhs._nset) {
    if (_nset == rhs._nset) {  // no bits set in either
      return static_cast<similarity_type_t>(0.0);
    }

    if (!tv.nset_sensitive_zero_bit_similarity()) {
      return static_cast<similarity_type_t>(1.0);
    }

    if (_nset) {
      return static_cast<similarity_type_t>(_nset) /
             static_cast<similarity_type_t>(_nset + 1);
    } else {
      return static_cast<similarity_type_t>(rhs._nset) /
             static_cast<similarity_type_t>(rhs._nset + 1);
    }
  }

  int c = ::bits_in_common((const unsigned int *)_bits, (const unsigned int *)rhs._bits,
                           _whole_words);

  if (0 == c) {
    return static_cast<similarity_type_t>(1.0);
  }

  // Use Bradshaw's notation

  float a = static_cast<float>(_nset - c);
  float b = static_cast<float>(rhs._nset - c);

  // cerr << "a = " << a << " b = " << b << " c = " << c << '\n';

  // similarity_type_t tv1 = float(c) / (tv.a() * a + tv.b() * b + float(c));
  // similarity_type_t tv2 = float(c) / (tv.b() * a + tv.a() * b + float(c));
  // similarity_type_t tnm = float(c) / static_cast<float>(_nset + rhs._nset - c);
  const similarity_type_t tv1 = Fraction(c, tv.a() * a + tv.b() * b + float(c));
  const similarity_type_t tv2 = Fraction(c, tv.b() * a + tv.a() * b + float(c));
  const similarity_type_t tnm = Fraction(c, static_cast<float>(_nset + rhs._nset - c));

  // These are similarities right now, so return 1.0 - the largest similarity

  // cerr << "Tv1 " << tv1 << ", tv2 " << tv2 << " tn " << tnm << '\n';
  if (tv1 >= tv2 && tv1 >= tnm) {
    return static_cast<similarity_type_t>(1.0) - tv1;
  }

  if (tv2 >= tv1 && tv2 >= tnm) {
    return static_cast<similarity_type_t>(1.0) - tv2;
  }

  return static_cast<similarity_type_t>(1.0) - tnm;
}

/*
  Sometimes it is more convenient to specify the tversky coefficients
  up front, and then call the function with no arguments
*/

static tversky_coeff_t global_tversky_alpha = 1.0;
static tversky_coeff_t global_tversky_beta = 1.0;

int
set_global_tversky_alpha(tversky_coeff_t na) {
  global_tversky_alpha = na;

  return 1;
}

int
set_global_tversky_beta(tversky_coeff_t nb) {
  global_tversky_beta = nb;

  return 1;
}

int
set_global_tversky_coefficients(tversky_coeff_t na, tversky_coeff_t nb) {
  global_tversky_alpha = na;
  global_tversky_beta = nb;

  return 1;
}

similarity_type_t
IWDYFP::tversky(IWDYFP &rhs) {
  assert(_whole_bytes == rhs._whole_bytes);

  if (_nset < 0) {
    (void)nset();
  }

  if (rhs._nset < 0) {
    (void)rhs.nset();
  }

  if (0 == _nset || 0 == rhs._nset) {
    if (0 == _nset && 0 == rhs._nset) {
      return static_cast<similarity_type_t>(1.0);
    }

    return static_cast<similarity_type_t>(0.0);
  }

  int c = ::bits_in_common((const unsigned int *)_bits, (const unsigned int *)rhs._bits,
                           _whole_words);

  // Use Bradshaw's notation

  int a = _nset - c;
  int b = rhs._nset - c;

  return (float(c) /
          (global_tversky_alpha * float(a) + global_tversky_beta * float(b) + float(c)));
}

similarity_type_t
IWDYFP::fraction_matched(IWDYFP &rhs) {
  (void)nbits();

  int nc;
  if (rhs.nbits() > _nbits) {
    nc = _bits_in_common(rhs);
  } else {
    nc = rhs._bits_in_common(*this);
  }

  return Fraction(nc, _nbits);
}

similarity_type_t
tanimoto(IWDYFP &fp1, IWDYFP &fp2) {
  return fp1.tanimoto(fp2);
}

similarity_type_t
fraction_matched(IWDYFP &fp1, IWDYFP &fp2) {
  return fp1.fraction_matched(fp2);
}

similarity_type_t
tversky(IWDYFP &fp1, IWDYFP &fp2) {
  return fp1.tversky(fp2);
}

/*
  The only complication is the need to maintain the _nset value
*/

void
IWDYFP::operator+=(const IWDYFP &rhs) {
  IW_Bits_Base::operator+=(rhs);

  if (IWDYFP_NSET_NOT_COMPUTED == _nset) {
    ;
  } else if (IWDYFP_NSET_NOT_COMPUTED == rhs._nset) {
    _nset = IWDYFP_NSET_NOT_COMPUTED;
  } else {
    _nset += rhs._nset;
  }

  return;
}

int
IWDYFP::fold(int nfold) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;

  _whole_words = _whole_words / 2;

  return IW_Bits_Base::fold(nfold);
}

/*
  We overload this method so we can compute _nset at the same time
*/

int
IWDYFP::construct_from_array_of_ints(const int *ii, int nb) {
  assert(nb > 0);

  if (_nbits) {  // zero out any pre-existing information
    clear();
  }

  allocate_space_for_bits(nb);
  _determine_whole_words();

  _nset = 0;

  for (int i = 0; i < _whole_bytes; i++) {
    for (int j = 0; j < IW_BITS_PER_BYTE; j++) {
      if (*ii) {
        _bits[i] |= one_bit_8[j];
        _nset++;
      }
      ii++;
    }
  }

  if (_extra_bits) {
    for (int i = 0; i < _extra_bits; i++) {
      if (*ii) {
        _bits[_whole_bytes - 1] |= one_bit_8[i];
        _nset++;
      }
      ii++;
    }
  }

  return 1;
}

int
IWDYFP::construct_from_descriptor_tdt_record(const IWString &buffer) {
  assert(buffer.ends_with('>'));

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

IWDYFP::construct_from_descriptor_record(const const_IWSubstring &buffer) {
  if (_nbits) {
    clear();
  }

  _nset = 0;

  int nb = buffer.nwords();

  allocate_space_for_bits(nb);
  _determine_whole_words();

  int j = 0;
  const_IWSubstring token;
  for (int i = 0; i < _nbits; i++) {
    if (!buffer.nextword(token, j)) {
      cerr << "IW_DY_Fingerprint::construct_from_descriptor_record: not enough tokens\n";
      cerr << buffer << '\n';
      return 0;
    }

    if ('0' != token) {
      set(i);
      _nset++;
    }
  }

  return 1;
}

int
IWDYFP::construct_from_hex(const const_IWSubstring &zhex) {
  _nset = IWDYFP_NSET_NOT_COMPUTED;
  if (!IW_Bits_Base::construct_from_hex(zhex)) {
    return 0;
  }

  _determine_whole_words();

  nset();

  return 1;
}

similarity_type_t
IWDYFP::fvb_modified_tanimoto(IWDYFP &rhs) {
  int ns1 = nset();
  int ns2 = rhs.nset();

  int n11;  // bits set in both

  if (0 == ns1 || 0 == ns2) {
    n11 = 0;
  } else {
    n11 = _bits_in_common(rhs);
  }

  int n00 = _nbits - ns1 - ns2 + n11;  // bits set in neither

  return fligner_verducci_blower(_nbits, _nset, rhs._nset, n00, n11);
}

similarity_type_t
IWDYFP::russel_rao(IWDYFP &rhs) {
  if (0 == nset() || 0 == rhs.nset()) {
    return static_cast<similarity_type_t>(0.0);
  }

  int n11 = _bits_in_common(rhs);

  return Fraction(n11, _nbits);
}

similarity_type_t
IWDYFP::forbes_similarity(IWDYFP &rhs) {
  int ns1 = nset();
  int ns2 = rhs.nset();

  if (0 == ns1 || 0 == ns2) {
    return static_cast<similarity_type_t>(0.0);
  }

  int a = _bits_in_common(rhs);
  if (0 == a) {
    return static_cast<similarity_type_t>(0.0);
  }

  // int b = ns1 - a;    // just set in just LHS
  // int c = ns2 - a;    // just set in just RHS

  return Fraction(a, ns1 * ns2);
}

similarity_type_t
IWDYFP::simple_matching(IWDYFP &rhs) {
  int ns1 = nset();
  int ns2 = rhs.nset();

  int a = _bits_in_common(rhs);

  int d = _nbits - ns1 - ns2 + a;  // bits set in neither

  return Fraction(a + d, _nbits);
}

similarity_type_t
IWDYFP::sorensendice(IWDYFP &rhs) {
  const int bic = _bits_in_common(rhs);

  return 1.0f - Fraction(bic, nset() + rhs.nset());
}

similarity_type_t
IWDYFP::overlap(IWDYFP &rhs) {
  const int bic = _bits_in_common(rhs);
  if (0 == bic) {
    return static_cast<similarity_type_t>(1.0);
  }

  const int ns1 = nset();
  const int ns2 = rhs.nset();

  if (ns1 < ns2) {
    return static_cast<similarity_type_t>(1.0) - Fraction(bic, ns1);
  } else {
    return static_cast<similarity_type_t>(1.0) - Fraction(bic, ns2);
  }
}

void *
IWDYFP::copy_to_contiguous_storage(void *p) const {
  // const unsigned char * initp = reinterpret_cast<const unsigned char *>(p);
  p = IW_Bits_Base::copy_to_contiguous_storage(p);

  size_t sbase = sizeof(IW_Bits_Base);

  memcpy(p, reinterpret_cast<const unsigned char *>(this) + sbase,
         sizeof(IWDYFP) - sbase);

  p = reinterpret_cast<unsigned char *>(p) + sizeof(IWDYFP) - sbase;

  return p;
}

void *
IWDYFP::copy_to_contiguous_storage_gpu(void *p) const {
  // const unsigned char * initp = reinterpret_cast<const unsigned char *>(p);
  p = IW_Bits_Base::copy_to_contiguous_storage_gpu(p);

  memcpy(p, &_nset, sizeof(int));

  p = reinterpret_cast<unsigned char *>(p) + sizeof(int);

  return p;
}

const void *
IWDYFP::build_from_contiguous_storage(const void *p, int allocate_arrays) {
  p = IW_Bits_Base::build_from_contiguous_storage(p, allocate_arrays);

  size_t sbase = sizeof(IW_Bits_Base);

  memcpy(reinterpret_cast<unsigned char *>(this) + sbase, p, sizeof(IWDYFP) - sbase);

  return reinterpret_cast<const unsigned char *>(p) + sizeof(IWDYFP) - sbase;
}

// The maximum similarity between two fingerprints is if the number of
// bits in common is the same as the number set in the least populated
// bitvector.
// Some simple math with the Tanimoto formula leads to this.
int
IWDYFP::SimilarityMightBeGreaterThan(const IWDYFP &rhs,
                                     similarity_type_t threshold) const {
  if (_nset == 0 || rhs._nset == 0) {
    return 0;
  }

  if (_nset <= rhs._nset) {
    return Fraction(_nset, rhs._nset) >= threshold;
  } else {
    return Fraction(rhs._nset, _nset) >= threshold;
  }
}
