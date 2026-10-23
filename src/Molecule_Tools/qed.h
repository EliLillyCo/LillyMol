#ifndef MOLECULE_TOOLS_QED_H_
#define MOLECULE_TOOLS_QED_H_

#include <optional>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/nvrtspsa.h"

namespace qed {

struct QEDProperties {
  public:
    float amw;
    float alogp;
    uint32_t hba;
    uint32_t hbd;
    float psa;
    uint32_t rotb;
    uint32_t arom;
    uint32_t alerts;

  public:
    void Reset();
};

struct ADSparameter {
  float a;
  float b;
  float c;
  float d;
  float e;
  float f;
  float dmax;
};

class Qed {
  private:
    resizable_array_p<Substructure_Query> _queries;

    // We can either use the externally specified acceptor queries
    // or our own computation.
    resizable_array_p<Substructure_Query> _acceptor_queries;

    // We can use the donor queries from rdkit.
    resizable_array_p<Substructure_Query> _donor_queries;

    // We can use the RDKit strict rotbond definition.
    resizable_array_p<Substructure_Query> _rdkit_rotb_queries;

    alogp::ALogP _alogp;

    nvrtspsa::NovartisPolarSurfaceArea _nvrtspsa;

    quick_rotbond::QuickRotatableBonds _rotbond;

  // Private functions
    int Alerts(Molecule& m);
    // Used for externally supplied donor and/or acceptor queries
    int ExternalQueryCount(Molecule& m,
                resizable_array_p<Substructure_Query>& queries);
    int ReadQueryFile(const IWString& fname,
                resizable_array_p<Substructure_Query>& destination,
                const char* description);
    int FinishInitialisation();

  public:
    Qed();

    int Initialise(Command_Line& cl, char flag);
    int InitialiseFromDirectory(const IWString& qed_dir);
    int InitialiseFromEnvironment();

    // WARNING - THIS MODIFIES `m`. It suppresses Hydrogens and reduces to the
    // largest organic fragment. Intended for use via qed() below, which
    // arranges for a copy where one is needed.
    bool CalculateProperties(Molecule& m, QEDProperties& result);

    // The main entry point. Compute QED.
    //
    // Does NOT modify `m`. The calculation itself is destructive - it
    // suppresses Hydrogens, reduces to the largest organic fragment, and
    // computes alogp, which strips isotopes. A molecule that any of those would
    // alter is copied first. A single organic fragment with no removable
    // Hydrogen and no isotope is scored in place.
    std::optional<float> qed(Molecule& m);

    // If the user has computed the properties, convert that to a QED score.
    float ComputeQed(const QEDProperties& properties) const;

  private:
    // The destructive form. qed() arranges that `m` is either a copy or a
    // molecule that this will not alter.
    std::optional<float> QedDestructive(Molecule& m);
};

std::optional<float> qed(Molecule& m);

}  // namespace qed

#endif  // MOLECULE_TOOLS_QED_H_
