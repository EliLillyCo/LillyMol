#ifndef MOLECULE_LIB_ATOM_PAIR_FINGERPRINT_H_
#define MOLECULE_LIB_ATOM_PAIR_FINGERPRINT_H_

#include <cstdint>
#include <tuple>
#include <unordered_map>

#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Molecule_Lib/molecule.h"

#define DISCERN_ATOM_PAIR_COLLISIONS

namespace atom_pair_fingerprint {

class AtomPairFingerprint {

  // In order to avoid passing long argument lists.
  class APairOfAtoms {
    public:
      Molecule& _m;
      const atom_type_t* _atype;
      atom_number_t _a1;
      atom_number_t _a2;
      int _d;
      uint64_t _bitnum;
      Sparse_Fingerprint_Creator& _sfc;
  
    public:
      APairOfAtoms(Molecule & m, const atom_type_t* atype, Sparse_Fingerprint_Creator& sfc) : _m(m), _atype(atype) , _sfc(sfc) {
      }

      void SetA1(const atom_number_t a) {
        _a1 = a;
      }

      int SetA2(const atom_number_t a2) {
        _a2 = a2;
        _d = _m.bonds_between(_a1, _a2);
        return _d;
      }

      const atom_number_t a1() const { return _a1;}
      const atom_number_t a2() const { return _a2;}
      Molecule& m() {
        return _m;
      }

      atom_type_t atype(const atom_number_t a) const {
        return _atype[a];
      }

      int d() const { return _d;}

      const Bond * BondBetweenAtoms() {
        return _m.bond_between_atoms(_a1, _a2);
      }

      std::tuple<atom_type_t, atom_type_t> OrderedAtomTypes() const {
        const atom_type_t t1 = _atype[_a1];
        const atom_type_t t2 = _atype[_a2];
        if (t1 > t2)
          return {t2, t1};
        return {t1, t2};
      }

      void HitBit(uint64_t b) {
        _bitnum = b;
        _sfc.hit_bit(b);
      }

      uint64_t bitnum() const {
        return _bitnum;
      }

      const IWString& smiles() {
        return _m.smiles();
      }
      const IWString& name() {
        return _m.name();
      }
  };

  public:
    using uint64_t = atom_type_t;
  private:
    // _min_separation will default to 1, which is 1 bond. Set to zero
    // to get bits corresponding to atom types themselves.
    int _min_separation;
    int _max_separation;

    // It is optional whether we handle bonded pairs specially.
    bool _fingerprint_bonded_atoms_with_btype;

    // 'Magic' numbers that govern just which bits get set.
    uint64_t _bond_magic1;
    uint64_t _bond_magic2;
    uint64_t _bond_magic3;
    uint64_t _bond_magic4;

    uint64_t _magic1;
    uint64_t _magic2;
    uint64_t _magic3;
    uint64_t _magic4;

    // If we are looking for bit collisions, writing the bits, or
    // looking for bit explanations, there will be a need to examine
    // every bit formed. Rather than check all the reasons, a single
    // bool covers all cases.
    bool _need_to_examine_bits_formed;

    // During tuning of magic parameters, it is useful to be able to
    // look for collisions.
    bool _check_for_collisions;

    uint64_t _collisions_found;

    // If check_for_collisions, a map from bit number to the components
    // that formed the bit, atype1 - dist - atype2
    std::unordered_map<uint64_t, std::tuple<uint64_t, int, uint64_t>> _bit_explanation;

    // It can be useful to just gather all atom pairs beyond _max_separation into a
    // single distance bucket.

    bool _include_out_of_range_separations;

    // We can accumulate how often each atom in the input is part of apair.
    bool _examine_coverage;

    int * _coverage;

    // A destination if we are writing bit meanings.
    IWString_and_File_Descriptor _stream_for_bit_meanings;

//  private functions.

    // Returns a hash for a bond - used when _fingerprint_bonded_atoms_with_btype.
    // No idea why std:: qualifier is needed!
    std::uint64_t _BondConstant(const Bond & b) const;

    // If we are checking for bit collisions, do that here.
    // Atom types t1 and t2 are separated by d bonds, and producted bit bitnum.
    void _DoCheckForCollisions(const APairOfAtoms& a_pair);

    void _UpdateCoverage(const APairOfAtoms& a_pair);
    // At the end of capturing coverage, transfer info to Molecule in a_pair.
    void _CoverageToMolecule(APairOfAtoms& a_pair);

    // In the case where all atoms are being processed, and there is only one fragment in the
    // molecule, call a version that avoids a lot of checks.
    int _FingerprintFast(Molecule& m, const atom_type_t* atom_type, Sparse_Fingerprint_Creator& sfc);

    // The more general form where all nprocess atoms in atoms_to_process will be
    // fingerprinted.
    int _Fingerprint(Molecule& m,
        const atom_number_t * atoms_to_process,
        const int nprocess,
        const atom_type_t* atom_type,
        Sparse_Fingerprint_Creator& sfc);

    // Bonded atoms a1 and a2 are to be fingerprinted, based on the bond between them.
    void _FingerprintBond(APairOfAtoms& a_pair);

    // Compute fingerprint for atoms a1 and a2, separated by d bonds.
    void _FingerprintPair(APairOfAtoms& a_pair);

    // All calculations come through this method.
    void _FormFingerprint(APairOfAtoms& a_pair);

    // If we need to examine the bits formed, b, that will be done here.
    // Atom number a1 and a2 are separated by d bonds.
    void _ExamineBit(APairOfAtoms& a_pair);

    // If writing bit meanings is enabled.
    void _DoWriteBitMeanings(APairOfAtoms& a_pair);

  public:
    AtomPairFingerprint();
    ~AtomPairFingerprint();

    void set_fingerprint_bonded_atoms_with_btype(bool s) {
      _fingerprint_bonded_atoms_with_btype = s;
    }

    void set_check_for_collisions(bool s) {
      _check_for_collisions = s;
      if (_check_for_collisions)
        _need_to_examine_bits_formed = true;
    }

    void set_include_out_of_range_separations(bool s) {
      _include_out_of_range_separations = s;
    }

    void set_examine_coverage(bool s) {
      _examine_coverage = s;
    }

    // note that we do not check these here, because it could
    // make things hard to use.
    void set_min_separation(const int s) {
      _min_separation = s;
    }

    void set_max_separation(const int s) {
      _max_separation = s;
    }

    int OpenStreamForBitMeanings(const char * fname);

    // The only reason this is non-const is because it might be
    // looking for bit collisions, writing bits or looking for
    // specific bits.
    int Fingerprint(Molecule& m, const int * include_atom,
        const atom_type_t * atom_type,
        Sparse_Fingerprint_Creator& sfc);

    uint64_t Collisions() const { return _collisions_found;}

    int ReportCollisions(std::ostream& output) const;
};

}  // namespace atom_pair_fingerprint

#endif  // MOLECULE_LIB_ATOM_PAIR_FINGERPRINT_H_
