#ifndef MOLECULE_TOOLS_DICER_API_H_
#define MOLECULE_TOOLS_DICER_API_H_

#include <unordered_map>
#include <iostream>
#include <string>
#include <tuple>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwbits/fixed_bit_vector.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/dicer_lib.h"

namespace dicer_api {

struct PerMoleculeData {
  int _matoms;

  // The size of a FixedBitVector with at least _matoms bits.
  int _vector_size;

  // Indexed by _matoms*_matoms
  int* _can_break;

  // for each atom in the starting molecule.
  dicer_lib::USPVPTR* _usvptr;

  uint32_t* _atype;

  // Used during create_subset.
  int* _xref;

  // keep track of the atoms that comprise the fragments already found.
  resizable_array_p<fixed_bit_vector::FixedBitVector> _found;

  public:
    PerMoleculeData(Molecule& m);
    ~PerMoleculeData();

    int vector_size() const {
      return _vector_size;
    }

    int* can_break() {
      return _can_break;
    }

    dicer_lib::USPVPTR* uspvptr() {
      return _usvptr;
    }

    uint32_t* atype() {
      return _atype;
    }
    int* xref() {
      return _xref;
    }


    // Determine the initial atom numbers and return the corresponding
    // entry in _bond_number
    int BondNumber(const Atom& at1, const Atom& at2) const;

    int AssignAtomTypes(Molecule& m, Atom_Typing_Specification& atom_typing);

    int InitialiseAtomPointers(Molecule& m);

    // int CanBreak(atom_number_t a1, atom_number_t a2) const;
    int CanBreak(const Atom& at1, const Atom& at2) const;
    int CanBreak(const Bond& b) const;

    // Return true if this bitvector is in _seen;
    // If not, transfer from `fp` to _seen and return false.
    int AlreadyFound(std::unique_ptr<fixed_bit_vector::FixedBitVector>& fp);
};

class Dicer {
  private:
    // Queries that define the breakable bonds.
    resizable_array_p<Substructure_Query> _query;

    int _max_bonds_to_break;

    int _min_fragment_size;
    int _max_fragment_size;

    // By default, we do not break amide bonds.
    int _break_amide_bonds;

    // Do we allow C-C bonds to break.
    int _break_cc_bonds;

    // Break any bond to a ring.
    int _break_ring_chain_bonds;

    isotope_t _label_join_points;

    Atom_Typing_Specification _atom_typing;

    // We can speed tings up by not identifying duplicates.
    int _determine_fragment_counts;

    // Inclusion and exclusion criteria on the fragments
    resizable_array_p<Substructure_Query> _fragments_must_contain;
    resizable_array_p<Substructure_Query> _fragments_cannot_contain;

    // Work like recap.
    // Recap just identifies the bonds to break and then breaks all of them.
    int _work_like_recap;

    // We can optionally keep track of all fragments found across calls to Dice()
    int _accumulate_global_fragment_count;
    std::unordered_map<std::string, uint32_t> _global;

  // private functions.

    int IdentifyBondsToBreak(Molecule& m, int* can_break);
    int IdentifyBondsToBreakViaQueries(Molecule& m, int* can_break);
    int IdentifyBondsToBreakHardCodedRules(Molecule& m, int* can_break);
    int Dice(Molecule& m,
            PerMoleculeData& pmd,
            int bonds_broken,
            std::unordered_map<std::string, uint32_t>& fragments);
    int MaybeStoreFragment(Molecule& m, std::unordered_map<std::string, uint32_t>& fragments);
    std::tuple<isotope_t, isotope_t> BreakTheBond(Molecule& m, const atom_number_t a1, const atom_number_t a2) const;
    void Restore(Molecule& m, atom_number_t a1, atom_number_t a2, const std::tuple<isotope_t, isotope_t>& iso) const;
    void UpdateGlobalFragmentCount(const std::unordered_map<std::string, uint32_t>& frags);
    int MakeTwoFragments(Molecule& m,
                        const int* fragment_membership,
                        PerMoleculeData& pmd,
                        int bonds_broken,
                        std::unordered_map<std::string, uint32_t>& fragments);
    int Recap(Molecule& m, PerMoleculeData& pmd,
            std::unordered_map<std::string, uint32_t>& fragments);

  public:
    Dicer();

    void set_max_bonds_to_break(int s) {
      _max_bonds_to_break = s;
    }

    // Should check that min and max are consistent with each other
    void set_min_fragment_size(int s) {
      _min_fragment_size = s;
    }
    void set_max_fragment_size(int s) {
      _max_fragment_size = s;
    }

    void set_break_cc_bonds(bool s) {
      _break_cc_bonds = s;
    }
    void set_break_amide_bonds(int s) {
      _break_amide_bonds = s;
    }

    void set_determine_fragment_counts(int s) {
      _determine_fragment_counts = s;
    }

    void set_label_join_points(int s) {
      _label_join_points = s;
    }

    void set_accumulate_global_fragment_count(bool s) {
      _accumulate_global_fragment_count = s;
    }

    void set_work_like_recap(int s) {
      _work_like_recap = s;
    }

    // Requirements on the fragments generated. Not implemented yet.
    int AddFragmentRequirementSmarts(const std::string& smarts);
    int AddFragmentDisqualifierSmarts(const std::string& smarts);

    // Not implemented yet.
    int set_atom_type(const std::string& s);

    // Note that adding external queries turns off all the built-in rules.
    int AddBreakBondSmarts(const std::string& smarts);
    int AddBreakBondQuery(const std::string& fname);

    // This is applied to each query in _query. Note make sure you call this only
    // after queries have been read.
    int set_perceive_symmetry_equivalent_matches(int s);

    int Dice(Molecule& m, std::unordered_map<std::string, uint32_t>& fragments);

    const std::unordered_map<std::string, uint32_t> global_fragment_count() const {
      return _global;
    }
};

}  // dicer_api

#endif  // MOLECULE_TOOLS_DICER_API_H_
