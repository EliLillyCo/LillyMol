#ifndef MOLECULE_TOOLS_ALOGP_H_
#define MOLECULE_TOOLS_ALOGP_H_

#include <optional>

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/alogp.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/alogp.pb.h"
#else
#include "alogp.pb.h"
#endif

namespace alogp {

// Computation of alogp
// Scott A Wildman, Gordon Crippen
// J Chem Inf Comput Sci 1999, 39, 868-873

// The atom types defined in Table 1

enum AlogPAtype {
  kC1 = 0,
  kC2 = 1,
  kC3 = 2,
  kC4 = 3,
  kC5 = 4,
  kC6 = 5,
  kC7 = 6,
  kC8 = 7,
  kC9 = 8,
  kC10 = 9,
  kC11 = 10,
  kC12 = 11,
  kC13 = 12,
  kC14 = 13,
  kC15 = 14,
  kC16 = 15,
  kC17 = 16,
  kC18 = 17,
  kC19 = 18,
  kC20 = 19,
  kC21 = 20,
  kC22 = 21,
  kC23 = 22,
  kC24 = 23,
  kC25 = 24,
  kC26 = 25,
  kC27 = 26,
  kCS = 27,

  kH1 = 28,
  kH2 = 29,
  kH3 = 30,
  kH4 = 31,
  kHS = 32,

  kN1 = 33,
  kN2 = 34,
  kN3 = 35,
  kN4 = 36,
  kN5 = 37,
  kN6 = 38,
  kN7 = 39,
  kN8 = 40,
  kN9 = 41,
  kN10 = 42,
  kN11 = 43,
  kN12 = 44,
  kN13 = 45,
  kN14 = 46,
  kNS = 47,

  kO1 = 48,
  kO2 = 49,
  kO3 = 50,
  kO4 = 51,
  kO5 = 52,
  kO6 = 53,
  kO7 = 54,
  kO8 = 55,
  kO9 = 56,
  kO10 = 57,
  kO11 = 58,
  kO12 = 59,
  kOS = 60,

  kF = 61,
  kCl = 62,
  kBr = 63,
  kI = 64,

  kHal = 67,
  kP = 66,

  kS1 = 67,
  kS2 = 68,
  kS3 = 69,

  kMe1 = 70,
  kMe2 = 71,

  kZwit = 72,

  kBias = 73,

  kLast = 74   // this is not used as a parameter, just to size the array below.
};

// When doing repeated scoring we save various attributes of the molecule to avoid
// recomputing them when repeatedly scoring a molecule.
struct ForFastScoring {
  int* atype;
  int* hcount;
  int* htype;
  int zwit;

  ForFastScoring();
  ~ForFastScoring();

  int resize(int s);
};

struct PerMoleculeData;

struct AlogPParams {
  float value[kLast + 1];
};

class ALogP {
  private:
    uint64_t _molecules_processed = 0;

    AlogPParams _params;

    // We can apply isotopic labels to the molecules being processed.
    int _label_with_atom_type = 0;

    // Observe that the contribution for an acid oxygen is very strongly
    // negative which seems unrealistic. If this is set, the alcohol atom
    // constant will be used instead, which might lead to more realistic
    // varlues.
    int _use_alcohol_for_acid = 0;

    // When RDKit processes a charge amine, it seems to count the Hydrogen
    // atoms for the uncharged form, 2 vs 3 for [N+H3].
    int _rdkit_charged_nitrogen = 0;

    // RDKit does not classify the Hydrogens on phosphoric acids as acid.
    // By default, we classify them as acidic, although strictly speaking,
    // the paper does not include P in the acidic hydrogen query. But they
    // do include Sulphur.
    // But nobody much cares about phosphorus...
    int _rdkit_phoshoric_acid_hydrogen = 0;

    // Note that many of the worst predictions involve Zwitterionic forms.
    int _add_zwitterion_correction = 0;

    // By default, we fail if we encounter an otherwise unclassified atom.
    int _fail_if_unclassified_atom = 1;

    // Whether or not to display error messages.
    int _display_error_messages = 1;

  // Private functions
    void DefaultParameters();

    int Carbon(PerMoleculeData& pmd, atom_number_t zatom);
    int AromaticCarbon(PerMoleculeData& pmd, atom_number_t zatom);
    int SaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom);
    int UnSaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom);
    int SaturatedPrimaryCarbon(PerMoleculeData& pmd, atom_number_t zatom);
    int SaturatedSecondaryCarbom(PerMoleculeData& pmd, atom_number_t zatom);

    int Nitrogen(PerMoleculeData& pmd, atom_number_t zatom);
    int AromaticNitrogen(PerMoleculeData& pmd, atom_number_t zatom);
    int SinglyConnectedSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom);
    int SaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom);
    int UnSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom);

    int Oxygen(PerMoleculeData& pmd, atom_number_t zatom);
    int AromaticOxygen(PerMoleculeData& pmd, atom_number_t zatom);
    int SaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom);
    int UnSaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom);

    int Fluorine(PerMoleculeData& pmd, atom_number_t zatom);
    int Chlorine(PerMoleculeData& pmd, atom_number_t zatom);
    int Bromine(PerMoleculeData& pmd, atom_number_t zatom);
    int Iodine(PerMoleculeData& pmd, atom_number_t zatom);

    int Phosphorus(PerMoleculeData& pmd, atom_number_t zatom);

    int Sulphur(PerMoleculeData& pmd, atom_number_t zatom);

    int IsHydrogenAcid(PerMoleculeData& pmd, atom_number_t zatom);
    int AddHydrogenContributions(PerMoleculeData& pmd, float& result);

    int AddZwitterionCorrection(PerMoleculeData& pmd, float& result);

    std::optional<float> LogPInner(Molecule& m, PerMoleculeData& pmd);

    std::optional<double> SingleAtomSpecialCase(Molecule& m);

    int ParametersFromProto(const alogp::AlogpParameters& proto);

    int GetHydrogenContributions(Molecule& m, PerMoleculeData& pmd, ForFastScoring& for_fast_scoring);

  public:
    ALogP();
    
    // Read an alogp::AlogpParameters textproto with the contributions.
    int ReadConfiguration(IWString& fname);

    // If a proto is available copies parameters from that.
    int ConfigFromProto(const alogp::AlogpConfiguration& proto);

    void set_display_error_messages(int s) {
      _display_error_messages = s;
    }

    void set_label_with_atom_type(int s) {
      _label_with_atom_type = s;
    }

    void set_use_alcohol_for_acid(int s) {
      _use_alcohol_for_acid = s;
    }

    void set_rdkit_charged_nitrogen(int s) {
      _rdkit_charged_nitrogen = s;
    }

    void set_rdkit_phoshoric_acid_hydrogen(int s) {
      _rdkit_phoshoric_acid_hydrogen = s;
    }

    void set_apply_zwitterion_correction(int s) {
      _add_zwitterion_correction = s;
    }

    void set_fail_if_unclassified_atom(int s) {
      _fail_if_unclassified_atom = s;
    }

    // Overwrite the contents of _params.value
    // Fails unless `n` matches the expected dimensionality.
    // Intended for use with alogp_optimise.
    template <typename T>
    int SetWeights(uint32_t n, const T* values);

    // Note that molecules must have formal charges assigned.
    // This class does not check that. Things will silently yield
    // bad values if charges have not been applied.
    // Note that explicit Hydrogen atoms are removed from `m`.
    std::optional<float> LogP(Molecule& m);

    // During alogp_optimise the same molecules are repeatedly
    // scored with different weights, so the atom assignment is
    // known.

    float LogP(Molecule& m, const int* atype, const int* hcount, const int* htype);

    // When repeatedly scoring the same molecule, but with different values for
    // the parameters, alogp_optimise, fill a struct that can then be used for
    // quickly re-evaluating a molecule, but presumably with different parameters.
    int FillForFastScoring(Molecule& m, ForFastScoring& for_fast_scoring);

    // Given a ForFastScoring structure, return the logp estimate.
    float LogP(Molecule& m, const ForFastScoring& for_fast_scoring);
};

}  //namespace alogp

#endif // MOLECULE_TOOLS_ALOGP_H_
