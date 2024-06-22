#ifndef MOLECULE_TOOLS_ALOGP_H_
#define MOLECULE_TOOLS_ALOGP_H_

#include <optional>

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/alogp.h"

namespace alogp {

// Computation of alogp
// Scott A Wildman, Gordon Crippen
// J Chem Inf Comput Sci 1999, 39, 868-873

// The atom types defined in Table 1

enum AlogPAtype {
  kNone = 0,
  kC1 = 1,
  kC2 = 2,
  kC3 = 3,
  kC4 = 4,
  kC5 = 5,
  kC6 = 6,
  kC7 = 7,
  kC8 = 8,
  kC9 = 9,
  kC10 = 10,
  kC11 = 11,
  kC12 = 12,
  kC13 = 13,
  kC14 = 14,
  kC15 = 15,
  kC16 = 16,
  kC17 = 17,
  kC18 = 18,
  kC19 = 19,
  kC20 = 20,
  kC21 = 21,
  kC22 = 22,
  kC23 = 23,
  kC24 = 24,
  kC25 = 25,
  kC26 = 26,
  kC27 = 27,
  kCS = 28,

  kH1 = 29,
  kH2 = 30,
  kH3 = 31,
  kH4 = 32,
  kHS = 33,

  kN1 = 34,
  kN2 = 35,
  kN3 = 36,
  kN4 = 37,
  kN5 = 38,
  kN6 = 39,
  kN7 = 40,
  kN8 = 41,
  kN9 = 42,
  kN10 = 43,
  kN11 = 44,
  kN12 = 45,
  kN13 = 46,
  kN14 = 47,
  kNS = 48,

  kO1 = 49,
  kO2 = 50,
  kO3 = 51,
  kO4 = 52,
  kO5 = 53,
  kO6 = 54,
  kO7 = 55,
  kO8 = 56,
  kO9 = 57,
  kO10 = 58,
  kO11 = 59,
  kO12 = 60,
  kOS = 61,

  kF = 62,
  kCl = 63,
  kBr = 64,
  kI = 65,

  kHal = 66,
  kP = 67,

  kS1 = 68,
  kS2 = 69,
  kS3 = 70,

  kMe1 = 71,
  kMe2 = 72
};

struct PerMoleculeData;

class ALogP {
  private:
    int _molecules_processed = 0;

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

    // Whether or not to display error messages.
    int _display_error_messages = 1;

  // Private functions
    int Carbon(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int AromaticCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int SaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int UnSaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int SaturatedPrimaryCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int SaturatedSecondaryCarbom(PerMoleculeData& pmd, atom_number_t zatom, float& result);

    int Nitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int AromaticNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int SinglyConnectedSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int SaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int UnSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result);

    int Oxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int AromaticOxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int SaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int UnSaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result);

    int Fluorine(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int Chlorine(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int Bromine(PerMoleculeData& pmd, atom_number_t zatom, float& result);
    int Iodine(PerMoleculeData& pmd, atom_number_t zatom, float& result);

    int Phosphorus(PerMoleculeData& pmd, atom_number_t zatom, float& result);

    int Sulphur(PerMoleculeData& pmd, atom_number_t zatom, float& result);

    int IsHydrogenAcid(PerMoleculeData& pmd, atom_number_t zatom);
    int AddHydrogenContributions(PerMoleculeData& pmd, float& result);

    int AddZwitterionCorrection(PerMoleculeData& pmd, float& result);

    std::optional<double> SingleAtomSpecialCase(Molecule& m);

  public:
    ALogP();
    
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

    // Note that molecules must have formal charges assigned.
    // This class does not check that. Things will silently yield
    // bad values if charges have not been applied.
    // Note that explicit Hydrogen atoms are removed from `m`.
    std::optional<float> LogP(Molecule& m);
};

}  //namespace alogp

#endif // MOLECULE_TOOLS_ALOGP_H_
