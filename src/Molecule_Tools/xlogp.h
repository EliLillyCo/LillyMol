#ifndef MOLECULE_TOOLS_XLOGP_H_
#define MOLECULE_TOOLS_XLOGP_H_

#include <optional>

#include "Molecule_Tools/xlogp.pb.h"

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

// Read updated parameter values from a textproto file.
// Note that these over-write only the values specified in the proto.
int ReadNewFragmentParameters(IWString& fname);

// For debugging it can be handy to see in detail which types get
// assigned to each atom.
void SetDisplayAtomAssignments(int s);

// It can be useful to suppress unclassified atom warning messages.
void SetIssueUnclassifiedAtomMessages(int s);

// If all atoms are classified, return xlogp.
std::optional<double> XLogP(Molecule& m);
// The user can have the atom types returned.
std::optional<double> XLogP(Molecule& m, int* status);

// For testing, it can be convenient to turn off the corrections.
// Do not use.

void ForTestingSetApplyCorrections(int s);
}  // namespace xlogp

#endif // MOLECULE_TOOLS_XLOGP_H_
