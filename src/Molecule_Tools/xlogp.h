#ifndef MOLECULE_TOOLS_XLOGP_H_
#define MOLECULE_TOOLS_XLOGP_H_

#include <array>
#include <optional>

#ifdef BUILD_BAZEL
#include "Molecule_Tools/xlogp.pb.h"
#else
#include "xlogp.pb.h"
#endif

#ifdef DOES_THIS_WORK_WELL_ENOUGH
I AM NOT SURE WHETHER THIS IS GOOD ENOUGH FOR USE.
There are still some cases that are predicted very poorly,
primarily involving charged Nitrogen atoms. I have 
corrections for them, but that has not solved that problem.

But if we compare Biobyte and Marvin, we see that those
two methods compare with each other about the same as
this method compares with them. AE95 is very high for
the Biobyte comparison, but overall RMS is ok. The
comparison with the mean of Biobyte and Marvin is not
too bad at all.

Compare to Marvin
 RMS: 0.939
 R2: 0.801
 AE95: 1.842

Compare to Biobyte
 RMS: 1.154
 R2: 0.727
 AE95: 2.233

Compare Biobyte with Marvin
 RMS: 0.981
 R2: 0.808
 AE95: 1.923

Compare to the mean of BioByte and Marvin
 RMS: 0.931
 R2: 0.803
 AE95: 1.779

So maybe this is not so bad...

Compare to some measured logD values - not logD and not logP so
we expect significant differences.

xlogp:
    RMS: 2.411
    R2: 0.163
    AE95: 4.39

Marvin:
    RMS: 2.43
    R2: 0.322
    AE95: 4.04

BioByte:
    RMS: 2.00
    R2: 0.438
    QE95: 3.71

RDKit:
    RMS: 2.66
    R2: 0.211
    AE95: 4.408

So in terms of concordance with 'experimental' values, BioByte
is the winner. But this xlogp appears to be competitive with Marvin.

#endif // DOES_THIS_WORK_WELL_ENOUGH

namespace xlogp {

struct PerMoleculeData;

class XLogPCalc {
  private:
    friend struct PerMoleculeData;

    static constexpr int kNumberFragmentTypes = 88;
    std::array<double, kNumberFragmentTypes> _type_to_score;

    bool _apply_corrections = true;
    bool _apply_nitroxide = true;
    bool _display_unclassified_atom_messages = true;
    bool _display_assignments = false;

    // `m` must have no explicit Hydrogens, or be a copy that may be altered.
    std::optional<double> LogPHSuppressed(Molecule& m, int* status) const;

    int ProcessNewFragmentParameter(const XLogP::XlogpParameter& proto);
    int ReadNewFragmentParameters(const XLogP::XlogpParameters& proto);

  public:
    XLogPCalc();

    // Read updated parameter values from a textproto file.
    // Only values specified in the proto are overwritten.
    int ReadNewFragmentParameters(IWString& fname);

    void SetIssueUnclassifiedAtomMessages(bool s) {
      _display_unclassified_atom_messages = s;
    }
    void SetDisplayAtomAssignments(bool s) {
      _display_assignments = s;
    }

    // The calculation mutates `m`: explicit hydrogen atoms are removed,
    // aromaticity is recomputed with the Wang-Fu-Lai model, and then
    // recomputed after restoring the previous global aromaticity model.
    // Neither of these modifies `m`. The calculation needs the Hydrogen
    // suppressed graph, so a molecule carrying explicit Hydrogens is copied and
    // stripped first.
    //
    // `status` is filled per atom. When a copy is made its indices refer to the
    // Hydrogen suppressed molecule, not to `m` as the caller has it, since
    // removing Hydrogens renumbers the atoms. A caller wanting them aligned
    // should suppress the Hydrogens itself first - see xlogp_main.
    std::optional<double> LogP(Molecule& m, int* status) const;
    std::optional<double> LogP(Molecule& m) const;

    // For testing and diagnostics.
    void ForTestingSetApplyCorrections(bool s) {
      _apply_corrections = s;
    }
    void SetApplyNitroxide(bool s) {
      _apply_nitroxide = s;
    }
};

}  // namespace xlogp

#endif // MOLECULE_TOOLS_XLOGP_H_
