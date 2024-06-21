#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <memory>
#include <nmmintrin.h>

#include "Foundational/iwbits/dy_fingerprint.h"
#include "Foundational/iwbits/iwbits.h"

#include "fixed_bit_vector.h"

namespace fixed_bit_vector {

using std::cerr;

/* Each single bit. Generate via
a = one(UInt64)
for o in 0:63
  b = IOBuffer()
  show(b, a << o)
  println(String(take!(b)))
end
*/

constexpr int kBitsPerWord = 64;
constexpr int kBytesPerWord = 8;

static uint64_t one_bit_64[] = {
  0x1,
  0x2,
  0x4,
  0x8,
  0x10,
  0x20,
  0x40,
  0x80,
  0x100,
  0x200,
  0x400,
  0x800,
  0x1000,
  0x2000,
  0x4000,
  0x8000,
  0x10000,
  0x20000,
  0x40000,
  0x80000,
  0x100000,
  0x200000,
  0x400000,
  0x800000,
  0x1000000,
  0x2000000,
  0x4000000,
  0x8000000,
  0x10000000,
  0x20000000,
  0x40000000,
  0x80000000,
  0x100000000,
  0x200000000,
  0x400000000,
  0x800000000,
  0x1000000000,
  0x2000000000,
  0x4000000000,
  0x8000000000,
  0x10000000000,
  0x20000000000,
  0x40000000000,
  0x80000000000,
  0x100000000000,
  0x200000000000,
  0x400000000000,
  0x800000000000,
  0x1000000000000,
  0x2000000000000,
  0x4000000000000,
  0x8000000000000,
  0x10000000000000,
  0x20000000000000,
  0x40000000000000,
  0x80000000000000,
  0x100000000000000,
  0x200000000000000,
  0x400000000000000,
  0x800000000000000,
  0x1000000000000000,
  0x2000000000000000,
  0x4000000000000000,
  0x8000000000000000,
};

#ifdef NOT_USED_MAYBE_USEFUL_SOMETIME
  static const uint64_t bx[] = {
    0xFFFFFFFF00000000,
    0xFFFF0000,
    0xFF00,
    0xF0,
    0xC,
    0x2,
    };

static const uint8_t BitReverseTable256[256] = 
{
#   define R2(n)     n,     n + 2*64,     n + 1*64,     n + 3*64
#   define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#   define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
    R6(0), R6(2), R6(1), R6(3)
};
#endif

// Return the number of 64 bit words needed for `nbits` bits.
int
words_for_bits(int nbits) {
  int result = nbits / kBitsPerWord;
  if (nbits % kBitsPerWord == 0)
    return result;
  return result + 1;
}

// Initialise our data structures to handle `nb` bits.
// If it is not a multiple of 64, it is rounded up.
void
FixedBitVector::_allocate_bits(int nb) {
  if (nb == 0) {
    resize(0);
    return;
  }

  _nwords = words_for_bits(nb);
  _bits = new uint64_t[_nwords];
  std::fill_n(_bits, _nwords, 0);
  _nbits = _nwords * kBitsPerWord;
}


void
FixedBitVector::_default_values() {
  _bits = nullptr;
  _nwords = 0;
  _nbits = 0;
}

FixedBitVector::FixedBitVector() {
  _default_values();
}

FixedBitVector::FixedBitVector(int nb) {
  if (nb == 0) {
    _default_values();
    return;
  }

  _nwords = words_for_bits(nb);
  _allocate_bits(nb);
}

FixedBitVector::FixedBitVector(const FixedBitVector& rhs) {
  _default_values();
  if (rhs._bits == nullptr) {
    return;
  }

  _allocate_bits(rhs._nwords * kBitsPerWord);
  _nwords = rhs._nwords;
  _nbits = rhs._nbits;
  std::copy_n(rhs._bits, _nwords, _bits);
}

FixedBitVector::FixedBitVector(FixedBitVector&& rhs) {
  _default_values();
  if (rhs._bits == nullptr) {
    return;
  }

  _nwords = rhs._nwords;
  _nbits = rhs._nbits;
  _bits = rhs._bits;

  delete [] rhs._bits;
  rhs._bits = nullptr;
  rhs._nwords = 0;
  rhs._nbits = 0;
}

int
FixedBitVector::DebugPrint(std::ostream& output) const {
  output << "FixedBitVector:_nwords " << _nwords;
  PrintOn(output);
  output << '\n';
  return output.good();
}

int
FixedBitVector::PrintOn(std::ostream& output, char t, char f, char sep) const {
  output << To01tf(t, f, sep);
  return output.good();
}

FixedBitVector
FixedBitVector::operator= (const FixedBitVector& rhs) {
  if (this == &rhs) {  // Self assignment is a no-op.
    return *this;
  }

  if (rhs._bits == nullptr) {
    _allocate_bits(0);
    return *this;
  }
  _allocate_bits(rhs._nwords * kBitsPerWord);
  _nwords = rhs._nwords;
  _nbits = rhs._nbits;
  std::copy_n(rhs._bits, _nwords, _bits);
  return *this;
}

FixedBitVector::~FixedBitVector() {
  if (_bits != nullptr)
    delete [] _bits;
}

void
FixedBitVector::resize(int nb) {
  if (_bits != nullptr)
    delete [] _bits;

  if (nb == 0) {
    _default_values();
    return;
  }

  _allocate_bits(nb);
}

bool
FixedBitVector::is_set(int b) const {
  return _bits[b / kBitsPerWord] & one_bit_64[b % kBitsPerWord];
}

void
FixedBitVector::set_bit(int b) {
  _bits[b / kBitsPerWord] |= one_bit_64[b % kBitsPerWord];
}

void
FixedBitVector::unset_bit(int b) {
  _bits[b / kBitsPerWord] &= ~one_bit_64[b % kBitsPerWord];
}

int
FixedBitVector::nset() const {
  int rc = 0;
  for (int i = 0; i < _nwords; ++i) {
    rc +=  _mm_popcnt_u64(_bits[i]);
  }
  return rc;
}

void
FixedBitVector::reset() {
  if (_bits == nullptr) {
    return;
  }
  std::fill_n(_bits, _nwords, 0);
}

// Could possibly be made more efficient with loop unrolling, see gfp_standard.cc
int
FixedBitVector::BitsInCommon(const FixedBitVector& rhs) const {
  int rc = 0;
  for (int i = 0; i < _nwords; ++i) {
    rc +=  _mm_popcnt_u64(_bits[i] & rhs._bits[i]);
  }
  return rc;
}

// Taken from http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious. Amazing!

int
most_significant_bit(uint64_t v)
{
  static const uint64_t b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000, 0xFFFFFFFF00000000};
  static const uint64_t S[] = {1, 2, 4, 8, 16, 32};

  uint64_t r = 0; // result of first_bit_set(v) will go here
  for (int i = 5; i >= 0; i--) // unroll for speed...
  {
    if (v & b[i])
    {
      v >>= S[i];
      r |= S[i];
    }
  }

  return r;
}

int
first_bit_set(uint64_t v) {
  static const uint64_t b[] = {
    0x00000000FFFFFFFF,
    0x0000FFFF,
    0x00FF,
    0x0F,
    0x3,
    0x1,
  };
  static const uint64_t S[] = {32, 16, 8, 4, 2, 1};

  if (v == 0) {
    return -1;
  }

  uint64_t result = 0;
  for (int i = 0; i < 6; ++i)
  {
    const uint64_t vnbi = v & b[i];
    if (vnbi == 0)
    {
      v >>= S[i];
      result |= S[i];
    }
    else if (vnbi == b[i]) {
      return result;
    }
  }

  return result;
}

// Note this could be optimized if we knew low order bits were preferentially set.
int
first_unset_bit(uint64_t v)
{

  static const uint64_t b[] = {
    0x00000000FFFFFFFF,
    0x0000FFFF,
    0x00FF,
    0x0F,
    0x3,
    0x1,
  };
  static const uint64_t S[] = {32, 16, 8, 4, 2, 1};

  if (v == 0) {
    return 0;
  }
  if (v == std::numeric_limits<uint64_t>::max()) {
    return -1;
  }

  uint64_t result = 0;
  for (int i = 0; i < 6; ++i)
  {
    const uint64_t vnbi = v & b[i];
    if (vnbi == 0) {
      return result;
    }
    if (vnbi == b[i])  // All bits set, shift
    {
      v >>= S[i];
      result |= S[i];
    }
  }

  return result;
}

int
FixedBitVector::FirstBitSet() const {
  for (int i = 0; i < _nwords; ++i) {
    if (_bits[i] == 0)
      continue;
    return i * kBitsPerWord + first_bit_set(_bits[i]);
  }
  return -1;
}

std::ostream&
operator<<(std::ostream& output, const FixedBitVector& b) {
  output << "FixedBitVector " << b.nbits();
  output << std::hex;
  for (int i = 0; i < b._nwords; ++i) {
    output << ' ' << b._bits[i];
  }
  output << std::dec;

  return output;
}

int
FixedBitVector::ConstructFromHex(const const_IWSubstring& hex) {
  if (hex.empty()) {
    resize(0);
    return 1;
  }

  // Two hex chars describe each byte, which is 8 bits.
  int nbits = hex.length() * 4;

  int nwords = nbits / kBitsPerWord;
  if (nwords * kBitsPerWord != nbits) {
    nwords++;
  }

  if (nwords != _nwords) {
    resize(nbits);
  } else {
    clear();  // Very important, we may not get a whole word of hex data.
  }

  return HexToBinary(hex, reinterpret_cast<unsigned char *>(_bits));
}

int
FixedBitVector::ConstructFromDaylightAsciiBitRep(const const_IWSubstring& ascii) {
  if (_bits) {
    resize(0);
  }

  unsigned int nbytes = (ascii.length() - 1) / 4;
  nbytes *= 3;
  int i = ascii[ascii.length() - 1] - '0';
  nbytes -= (3 - i);

  int number_bits = nbytes * IW_BITS_PER_BYTE;

  if (number_bits == 0) {
    return 1;
  }

  resize(number_bits);

  unsigned int bytes_read;
  if (! (du_ascii2bin(ascii.data(), ascii.length(),
                      reinterpret_cast<unsigned char *>(_bits), bytes_read)))
  {
    cerr << "FixedBitVector::ConstructFromDaylightAsciiBitRep: du_ascii2bin failed\n";
    cerr << ascii << '\n';
    return 0;
  }

  return 1;
}

bool
FixedBitVector::operator== (const FixedBitVector& rhs) const {
  if (_nwords != rhs._nwords) {
    return false;
  }

  return std::equal(_bits, _bits + _nwords, rhs._bits);
}

IWString
FixedBitVector::DaylightAsciiRepresentation() const {
  int nchars = 0;
  const int nbytes = _nwords * kBytesPerWord;
  char * b = du_bin2ascii(&nchars, nbytes, reinterpret_cast<char *>(_bits));
  IWString result;
  result.set_and_assume_ownership(b, nchars);
  return result;
}

IWString
FixedBitVector::DaylightAsciiRepresentationIncludingNsetInfo() const {
  IWString result = DaylightAsciiRepresentation();

  const int ns = nset();

  constexpr char sep = ';';

  result << sep << _nbits << sep << ns << sep << _nbits << sep << ns << ";1";

  return result;
}

int
FixedBitVector::Fold(int nfold) {
  assert (nfold > 0);

  int nw2 = _nwords / 2;

  for (int i = 0; i < nw2; i++)
  {
    _bits[i] = (_bits[i]) | (_bits[nw2 + i]);
  }

  _nbits = _nbits / 2;
  _nwords = _nwords / 2;

  if (1 == nfold)
    return 1;

  return Fold(nfold - 1);
}

FixedBitVector&
FixedBitVector::operator += (const FixedBitVector& rhs) {
  if (rhs.nbits() == 0) {
    return *this;
  }

  const int old_nwords = _nwords;
  int new_bits = _nbits + rhs._nbits;
  std::unique_ptr<uint64_t[]> save_previous_info(std::make_unique<uint64_t[]>(_nwords));
  std::copy_n(_bits, _nwords, save_previous_info.get());

  resize(new_bits);
  std::copy_n(save_previous_info.get(), old_nwords, _bits);
  std::copy_n(rhs._bits, rhs._nwords, _bits + old_nwords);

  return *this;
}

int
FixedBitVector::ConstructFromAscii01Representation(const const_IWSubstring& ascii,
                        char t) {
  resize(0);

  if (ascii.empty()) {
    return 1;
  }

  int new_words = ascii.length() / kBitsPerWord;
  if (ascii.length() / kBitsPerWord * kBitsPerWord != ascii.length()) {
    new_words++;
  }

  resize(new_words * kBitsPerWord);
  for (int i = 0; i < ascii.length(); ++i) {
    char c = ascii[i];
    if (c == t) {
      set_bit(i);
    }
  }

  return 1;
}

int
FixedBitVector::ConstructFromTdtRecordNset(const const_IWSubstring& tdt_record, int &nset) {
//assert (tdt_record.ends_with('>'));

  nset = 0;

  const int openangle = tdt_record.index('<');

  assert (openangle > 0);

  const int semicolon = tdt_record.index(';');
  int i = semicolon;
  const_IWSubstring token_nbits, token_nset;
  if (! tdt_record.nextword(token_nbits, i, ';') || 
      ! tdt_record.nextword(token_nset, i, ';'))  {
    cerr << "FixedBitVector::ConstructFromTdtRecordNset:cannot extract nbits/nset\n";
    return 0;
  }

  int file_nbits, file_nset;
  if (! token_nbits.numeric_value(file_nbits) || file_nbits <= 0 ||
      ! token_nset.numeric_value(file_nset) || file_nset < 0 || 
      file_nset > file_nbits) {
    cerr << "FixedBitVector::ConstructFromTdtRecordNset:invalid bit data " << token_nbits << ';' << token_nset << '\n';
    return 0;
  }

  resize(file_nbits);

  const_IWSubstring just_bits(tdt_record);
  just_bits.iwtruncate(semicolon);
  just_bits.remove_leading_chars(openangle + 1);

  if (! ConstructFromDaylightAsciiBitRep(just_bits)) {
    cerr << "FixedBitVector::ConstructFromDaylightAsciiBitRep: inner call failed\n";
    return 0;
  }

  nset = file_nset;

  return 1;
}

void
FixedBitVector::Flip() {
  for (int i = 0; i < _nwords; i++)
  {
    _bits[i] = ~ (_bits[i]);
  }
}

void
FixedBitVector::iwor(const FixedBitVector & rhs) {
  if (_nwords != rhs._nwords) {
    cerr << "FixedBitVector:iwor: size mismatch " << _nwords << " vs " << rhs._nwords << '\n';
    return;
  }

  for (int i = 0; i < _nwords; ++i) {
    _bits[i] |= rhs._bits[i];
  }
}

void
FixedBitVector::iwor(const FixedBitVector & rhs, int & changed) {
  if (_nwords != rhs._nwords) {
    cerr << "FixedBitVector:iwor: size mismatch " << _nwords << " vs " << rhs._nwords << '\n';
    return;
  }

  changed = 0;
  uint64_t aword;
  for (int i = 0; i < _nwords; ++i) {
    aword = _bits[i] | rhs._bits[i];
    if (aword != _bits[i]) {
      _bits[i] = aword;
      changed = 1;
    }
  }
}

void
FixedBitVector::iwand(const FixedBitVector & rhs) {
  if (_nwords != rhs._nwords) {
    cerr << "FixedBitVector:iwand: size mismatch " << _nwords << " vs " << rhs._nwords << '\n';
    return;
  }

  for (int i = 0; i < _nwords; ++i) {
    _bits[i] &= rhs._bits[i];
  }
}

void
FixedBitVector::iwand(const FixedBitVector & rhs, int& changed) {
  if (_nwords != rhs._nwords) {
    cerr << "FixedBitVector:iwand: size mismatch " << _nwords << " vs " << rhs._nwords << '\n';
    return;
  }

  changed = 0;
  uint64_t aword;
  for (int i = 0; i < _nwords; ++i) {
    aword = _bits[i] & rhs._bits[i];
    if (aword != _bits[i]) {
      _bits[i] = aword;
      changed = 1;
    }
  }
}

void
FixedBitVector::iwxor(const FixedBitVector & rhs) {
  if (_nwords != rhs._nwords) {
    cerr << "FixedBitVector:iwand: size mismatch " << _nwords << " vs " << rhs._nwords << '\n';
    return;
  }

  for (int i = 0; i < _nwords; ++i) {
    _bits[i] ^= rhs._bits[i];
  }
}

int
FixedBitVector::ConstructFromSparseRepresentation(const const_IWSubstring& ascii) {
  cerr << "FixedBitVector::ConstructFromSparseRepresentation:not implemented\n";
  return 0;
}

// NOte that this could be done cheaper by just looking at the length of the string
// and assuming each bit takes 1 bytes.
int
FixedBitVector::ConstructFromAscii01RepresentationWithSpaces(const const_IWSubstring& ascii, char sep) {
  int nb = ascii.nwords(sep);
  if (nb == 0) {
    resize(0);
    return 1;
  }

  _allocate_bits(nb);

  int i = 0;
  const_IWSubstring token;
  for (int bit_num = 0; ascii.nextword(token, i, sep); ++bit_num) {
    if (token == '1') {
      set_bit(bit_num);
    }
  }

  return 1;
}

IWString
FixedBitVector::To01tf(char t, char f, char sep) const {
  IWString result;
  if (_nwords == 0) {
    return result;
  }

  if (sep == '\0') {
    result.resize(_nbits);
  } else {
    result.resize(_nbits * 2);
  }

  if (is_set(0)) {
    result << t;
  } else {
    result << f;
  }

  if (sep == '\0') {
    for (int i = 1; i < _nbits; ++i) {
      if (is_set(i)) {
        result << t;
      } else {
        result << f;
      }
    }
  } else {
    for (int i = 1; i < _nbits; ++i) {
      result << sep;
      if (is_set(i)) {
        result << t;
      } else {
        result << f;
      }
    }
  }

  return result;
}

int
FixedBitVector::HexForm(IWString& destination) const {
  destination.resize_keep_storage(nbits() / 4);    // each hex character encodes 4 bits

  // InternalHexForm takes arguments in bytes
  iwbits::InternalHexForm(reinterpret_cast<const unsigned char *>(_bits), 
                          _nwords * kBytesPerWord, destination);

  return 1;
}

IWString
FixedBitVector::HexForm() const {
  IWString result;
  HexForm(result);
  return result;
}

int
FixedBitVector::NextOnBit(int& ndx) const {
  int word_index = ndx / kBitsPerWord;

  int istart = ndx % kBitsPerWord;

  for (word_index = ndx / kBitsPerWord; word_index < _nwords; ++word_index) {
    const uint64_t zword = _bits[word_index];
    for (int i = istart; i < kBitsPerWord; ++i, ++ndx) {
      if (one_bit_64[i] & zword) {
        ndx++;
        return ndx - 1;
      }
    }
    istart = 0;
  }

  return -1;
}

int
FixedBitVector::next_on_bit(int & ndx) const {
  return NextOnBit(ndx);
}

void *
FixedBitVector::CopyToContiguousStorage(void * destination) const {
  if (_nwords == 0) {
    return destination;
  }

  uint64_t * ptr = reinterpret_cast<uint64_t *>(destination);
  std::copy_n(_bits, _nwords, ptr);
  return ptr + _nwords;
}

void*
FixedBitVector::BuildFromContiguousStorage(const void * source, int nw) {
  if (_nwords) {
    resize(0);
  }

  if (nw == 0) {
    return (uint64_t*) source;
  }

  _allocate_bits(nw * kBitsPerWord);
  const uint64_t* ptr = reinterpret_cast<const uint64_t*>(source);
  std::copy_n(ptr, nw, _bits);
  _nwords = nw;
  _nbits = nw * kBitsPerWord;

  return (uint64_t*) (ptr + nw);
}

template <typename T>
int
FixedBitVector::ConstructFromArray(const T * bits, int nbits) {
  if (_nwords) {
    resize(0);
  }
  _allocate_bits(nbits);
  for (int i = 0; i < nbits; ++i) {
    if (bits[i]) {
      set_bit(i);
    }
  }
  return nbits;
}

template int FixedBitVector::ConstructFromArray(const int * bits, int nbits);

}  // namespace fixed_bit_vector
