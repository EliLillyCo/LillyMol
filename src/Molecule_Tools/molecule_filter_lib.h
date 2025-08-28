#ifndef MOLECULE_TOOLS_MOLECULE_FILTER_LIB_H_
#define MOLECULE_TOOLS_MOLECULE_FILTER_LIB_H_

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rotbond_common.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/xlogp.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/molecule_filter.pb.h"
#else
#include "molecule_filter.pb.h"
#endif

namespace molecule_filter_lib {

class MoleculeFilter {
  private:
    quick_rotbond::QuickRotatableBonds _rotbond;

    alogp::ALogP _alogp;

    molecule_filter_data::Requirements _requirements;

    bool _active;

  // private functions
    void InitialiseOptionalFeatures();

  public:
    MoleculeFilter();

    // Read textproto configuration file
    int Build(IWString& fname);

    int active() const {
      return _active;
    }

    // Copy `proto` to _requirements and initialise.
    int Build(const molecule_filter_data::Requirements& proto);

    // Return true if `m` is consistent with the constratins set in `_requirements`.
    int Ok(Molecule& m);
};


// Various supporting functions used in the calculation.
int CountHeteroatoms(const Molecule& m);
int AromaticRingCount(Molecule& m);
bool LargestFragment(const const_IWSubstring& smiles,
                const_IWSubstring& largest_frag,
                int& natoms, int& nrings);
std::tuple<int, int> MaxRingSystemSize(Molecule& m, std::unique_ptr<int[]>& tmp);
void RuleOfFive(Molecule & m, int& acceptor, int& donor);
int HalogenCount(const Molecule& m);
int Sp3Carbon(Molecule & m);

} //   namespace molecule_filter_lib

#endif // MOLECULE_TOOLS_MOLECULE_FILTER_LIB_H_
