#ifndef FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H
#define FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H

#include <cstdint>
#include <iostream>

#include "Foundational/iwstring/iwstring.h"

namespace fixed_bit_vector {
// A BitVector class that only operates on 64 bit words.
// Deliberately minimal functionality. Designed for speed.
class FixedBitVector {
  protected:
    uint64_t * _bits;
    int _nwords;
    int _nbits;

  // Private functions
  private:
    void _default_values();
    void _allocate_bits(int nb);

  friend
    std::ostream& operator<<(std::ostream& output, const FixedBitVector& b);

  public:
    FixedBitVector();
    FixedBitVector(int nb);
    FixedBitVector(const FixedBitVector& rhs);
    FixedBitVector(FixedBitVector&& rhs);
    ~FixedBitVector();

    FixedBitVector operator=(const FixedBitVector& rhs);

    int DebugPrint(std::ostream& output) const;
    int PrintOn(std::ostream& output, char t='1', char f='0', char sep='\0') const;

    // Allocate space to accommodate `nb` bits.
    // Note this will be rounded up to the nearest 64 bits.
    // Existing data is discarded.
    void resize(int nb);

    int nbits() const { return _nbits;}

    void set_bit(int b);
    // Set many bits. Should work for any iteratble.
    template <typename T>
    void set_bits(const T& bits) {
      for (auto bit : bits) {
        set_bit(bit);
      }
    }

    void unset_bit(int b);
    bool is_set(int b) const;

    // Turn off all bits.
    void reset();
    void clear() {
      reset();
    }

    // The number of bits set.
    int nset() const;

    // The number of 64 bit words.
    int nwords() const {
      return _nwords;
    }

    // Flip bits;
    void Flip();

    // Bits in common with another vector. Must be same size.
    int BitsInCommon(const FixedBitVector& rhs) const;

    // Returns -1 if no bits are set.
    int FirstBitSet() const;
    // Returns -1 if all bits are set.
    int FirstUnsetBit() const;

    // Return true if bit `i` is set.
    bool operator[](int i) const;

    // Concatenate `rhs` onto this.
    FixedBitVector& operator += (const FixedBitVector& rhs);

    // Shorten by folding itself `nfold` times.
    int Fold(int nfold);

    // Build from '123456789abcdef'.
    // If not a multiple of 64 bits, silently zero pads to the
    // next word.
    int ConstructFromHex(const const_IWSubstring& hex);
    // Generate Hex encoded forms.
    int HexForm(IWString & destination) const;
    IWString HexForm() const;

    int ConstructFromDaylightAsciiBitRep(const const_IWSubstring& ascii);

    // Daylight fingerprints have the number of set bits appended.
    int ConstructFromTdtRecordNset(const const_IWSubstring& tdt_record, int &nset);

    // Build from '00111000'
    // The true character can be over ridden.
    int ConstructFromAscii01Representation(const const_IWSubstring& ascii, char t='1');
    // Build from '0 0 1 1 1 0 0 0'
    // Argument `sep` can be used to use a different separator.
    int ConstructFromAscii01RepresentationWithSpaces(const const_IWSubstring& ascii, char sep=' ');

    // Return a string with bits set as `t`, unset as `f` and with
    // `sep` (of not null) between characters.
    IWString To01tf(char t, char f, char sep = '\0') const;

    // The sparse representation looks like:
    //  1,5-7,88;nbits
    int ConstructFromSparseRepresentation(const const_IWSubstring& ascii);

    // Return encoded in Daylight form.
    IWString DaylightAsciiRepresentation() const;
    IWString DaylightAsciiRepresentationIncludingNsetInfo() const;

    // Items in `bits` that are non zero, set the corresponding bit.
    template <typename T> int ConstructFromArray(const T * bits, int nbits);

    // Copy our contents to `destination`.
    // Return a pointer just past the end of what was written.
    void* CopyToContiguousStorage(void * destination) const;
    // Copy `nwords` 64 bit words from `source` to build.
    // Returns a pointer to just past the end of what was read.
    void* BuildFromContiguousStorage(const void * source, int nwords);

    // Dangerous to let this out, but needed.
    const uint64_t* bits() const { return _bits;}

    bool operator== (const FixedBitVector& rhs) const;

    // Or this with `rhs`.
    void iwor(const FixedBitVector & rhs);
    void iwand(const FixedBitVector & rhs);
    void iwxor(const FixedBitVector & rhs);

    // Variants of bit operators with an extra parameter indicating whether or not
    // the operation has changed `this`.
    void iwand(const FixedBitVector& rhs, int& changed);
    void iwor(const FixedBitVector& rhs, int& changed);

//  Sequentially fetch all the 1 bits. Pass 0 at first, Returns -1 when done
    int next_on_bit(int& ndx) const;
    int NextOnBit(int& ndx) const;
};

std::ostream&
operator<<(const FixedBitVector& b, std::ostream& output);

// Returns the first bit set in a uint64_t. Exposed just for testing.
// Note that this is not the most significant bit, but the least significant bit.
int first_bit_set(uint64_t);
// Returns the first unset bit in a uint64_t. Note that it returns the first
// unset least significant bit.
int first_unset_bit(uint64_t);

}  // namespace fixed_bit_vector
#endif  // FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H
