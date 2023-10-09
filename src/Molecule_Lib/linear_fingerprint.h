#ifndef MOLECULE_LIB_LINEAR_FINGERPRINT_H_
#define MOLECULE_LIB_LINEAR_FINGERPRINT_H_

#include <cstdint>

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "molecule.h"

namespace linear_fingerprint {

using atom_type_t = std::uint64_t;

namespace internal {

// There are two classes for linear fingerprinting. One public
// and one internal. Both need a common set of optional behaviours.
struct Options
{
  // Min and max path length - in bonds.
  int _min_length;
  int _max_length;

  // Are paths allowed to cross.
  bool _paths_can_cross;

  // When rings are detected, produce ring presence bits.
  bool _fingerprint_ring_presence;

  // See whether or not different paths produce the same bit.
  bool _check_for_collisions;

  // Generate labelled molecules showing times atoms used.
  bool _check_coverage;

  Options();

  Options& operator =(const Options&) = default;
};

class LinearFpStatus {
  private:
    // From the constructor.
    Options _options;

    // From the Molecule being processed.
    const int _matoms;
    const int _nedges;
    
    // From constructor.
    const Molecule& _m;

    // Extracted from the molecule during constructor.
    const Atom** _atom;

    // We just keep a pointer to what comes in.
    const atom_type_t* _atype;

    const Bond_list& _bond_list;

    // For each bond, indexed by bond number, the invariant used.
    uint64_t * _bond_constant;

    // For each bond, is it in the current path or not
    int * _bond_in_path;

    // For each atom, is it in the current path or not
    int * _atom_in_path;

    // The sequence of atom-bond-atom-bond-atom- that comprise the
    // current state of the path.
    int _path_length;
    std::uint64_t * _path;

    // An index for each item in the path. Even numbers will be atom
    // numbers, and odd numbers will be bond numbers

    int * _path_index;

    std::uint64_t _magic1;
    std::uint64_t _magic2;
    std::uint64_t _magic3;

    std::uint64_t _bond_magic1;
    std::uint64_t _bond_magic2;
    std::uint64_t _bond_magic3;
    std::uint64_t _bond_magic4;

    bool _need_to_examine_bits_formed;

    // Keep track of how many times each atom involved in fp.
    int * _coverage;

    Sparse_Fingerprint_Creator& _sfc;

    // Passed from caller, we do not own it.
    IWString_and_File_Descriptor * _stream_for_bit_meanings;

    // Private functions.

    void _AddBondToPath(const Bond & b, const atom_number_t next_atom);

    // A hash value for a bond.
    uint64_t _BondHash(const Bond & b) const;

    // Add 'atom_number' to the path.
    void _StartPath(const int atom_number);

    // If possible, recursively expand from the current state of _path.
    void _Expand();

    // Decrease path length by one bond.
    void _PopPath();

    // Depending in _min_length and _max_length, produce a bit.
    void _MaybeFormBit();

    // A bit has been produced, do we need to do anything special with it?
    void _ExamineBit(const std::uint64_t);

    // A ring path has been encountered, make a bit.
    void _FormRingBit();

    // Generate a bit from the currenate state of _path.
    void _FormFingerprintForward();
    void _FormFingerprintBackward();

    // Update the _coverage array based on current _path.
    void _DoCheckCoverage();

    // At end of run, set the atom map numbers from the _coverage array.
    void _AssignCoverageDataToMolecule() const;

    void _PrintPath(std::ostream& output) const;

    // Write a labelled molecule for the current _path.
    void _WriteBit(const uint64_t b);

    // Special case.
    void _HandleMoleculeWithNoBonds();

    // Called once for each molecule so we have ready access to atom numbers
    // during debugging.
    void _WriteLabelledSmiles() const;
    
  public:
    LinearFpStatus(const Options& options, const Molecule& m,
                   const int * include_atom,
                   const atom_type_t * atype,
                   Sparse_Fingerprint_Creator& sfp);
    ~LinearFpStatus();

    void SetStreamForBitMeanings(IWString_and_File_Descriptor* s) {
      _stream_for_bit_meanings = s;
      _need_to_examine_bits_formed = true;
    }

    int Fingerprint();
};

}  // namespace internal

class LinearFingerprintGenerator {
  private:
    internal::Options _options;

    // Will be passed as pointer to individual LinearFpStatus that
    // do the actual fingerprinting.
    IWString_and_File_Descriptor _stream_for_bit_meanings;

  public:
    LinearFingerprintGenerator();

    void set_min_length(int s) { _options._min_length = s;}
    void set_max_length(int s) { _options._max_length = s;}
    void set_paths_can_cross(int s) { _options._paths_can_cross = s;}
    void set_check_for_collisions(int s) { _options._check_for_collisions = s;}
    void set_check_coverage(bool s) {_options._check_coverage = s;}
    void set_fingerprint_ring_presence(bool s) { _options._fingerprint_ring_presence = s;}

    int  OpenStreamForBitMeanings(const char * s);

    int Fingerprint(Molecule& m, const int * include_atom,
        const atom_type_t * atype,
        Sparse_Fingerprint_Creator& sfp);

    int ReportCollisions(std::ostream& output) const;
};

}  // namespace linear_fingerprint

#endif  // MOLECULE_LIB_LINEAR_FINGERPRINT_H_
