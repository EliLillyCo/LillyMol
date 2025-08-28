#ifndef MOLECULE_TOOLS_SAFE_GENERATE_LIB_H_
#define MOLECULE_TOOLS_SAFE_GENERATE_LIB_H_

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/mformula.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace safe_generate {

class SafeFragment {
  private:
    // the smiles of the fragment
    // Individual dot separated tokens from something like
    // [1Cl]%10.[1C]%11.[1C]%121=NO[1C]%13%11O1.[1C]%131CCCCN1.C1CC[1C]%10[1C]%12C1
    IWString _smiles;

    // The indices within _smiles where the first digit after the % signs are.
    resizable_array<int> _first_digit;

    // A molecule build from _smiles by dropping the %nn characters
    Molecule _m;

    mformula::MFormula _mformula;

    // Number of atoms in the fragment.
    int _natoms;

    // Number of rings in the fragment.
    int _nrings;

    // A list of ring sizes.
    resizable_array<int> _ring_sizes;

    // Number of connections (ring openings) in the fragment.
    int _ncon;

    // If two or more connections, the number of bonds between the closest two.
    int _distance;

    // Within _smiles the ring openings present.
    std::vector<int> _ring;

    // In a molecule consisting of individual SafeFragment's each fragment
    // can be marked as ok to select or not - depending on substructure
    // matches in the parent molecule.
    int _ok_to_select;

    // If atom type matching is in effect, we keep a list of the
    // isotopes (atom types). Note that this is complicated by
    // the fact that in order to actually match a fragment, we would
    // need to adjust the ring opening and closings to match the
    // atom types. Not sure we will do that with the SAFE generation
    // tool, but may wait for the Molecule based generator.
    // This is just a list of the isotopes at the join points. By
    // construction it will be same order as `_first_digit`.
    resizable_array<atom_type_t> _atypes;

  // Private functions
    int ProcessSquareBracket(const IWString& smi, int &i, int& next_ring,
                     IWString& destination);

  public:
    SafeFragment();

    int DebugPrint(std::ostream& output) const;

    const IWString& name() const {
      return _m.name();
    }
    void set_name(const IWString& s) {
      _m.set_name(s);
    }

    int natoms() const {
      return _natoms;
    }

    int ncon() const {
      return _ncon;
    }

    int nrings() const {
      return _nrings;
    }

    int distance() const {
      return _distance;
    }

    int Build(const const_IWSubstring& buffer);
    int Build(const dicer_data::DicerFragment& proto);

    int ok_to_select() const {
      return _ok_to_select;
    }
    void set_ok_to_select(int s) {
      _ok_to_select = s;
    }

    const IWString& smiles () const {
      return _smiles;
    }

    Molecule& mol() {
      return _m;
    }

    // Formula difference
    uint32_t FormulaDifference(const SafeFragment& rhs) const {
      return _mformula.Diff(rhs._mformula);
    }

    // Place into `new_smiles` the smiles of `f2` with the ring numbers
    // from `this`.
    int SameNumbers(const SafeFragment& f2, IWString& new_smiles) const;
};

// Examines all bonds in `m` and returns false if we find bonds that should
// be rejected.
int BondsOk(Molecule& m);

}  // namespace safe_generate

#endif  // MOLECULE_TOOLS_SAFE_GENERATE_LIB_H_
