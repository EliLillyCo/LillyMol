#ifndef MOLECULE_TOOLS_NVRTSPSA_H_
#define MOLECULE_TOOLS_NVRTSPSA_H_

#include <optional>

#include "Molecule_Lib/molecule.h"

class Command_Line;

namespace nvrtspsa {

// Class for performing Novartis Polar Surface Area calculations.
class NovartisPolarSurfaceArea {
  private:
    int _display_psa_unclassified_atom_messages = 1;
    int _return_zero_for_unclassified_atoms = 0;
    // The paper is uncertain on how certain S atoms are handled. By default,
    // LillyMol keeps the historical non-zero contribution for [SD2] atoms.
    int _non_zero_contribution_for_SD2 = 1;
    // RDKit compatibility options.
    int _zero_for_all_sulphur_atoms = 0;
    int _zero_for_all_phosphorus_atoms = 0;
    int _convert_to_charge_separated = 0;

    // Remove explicit Hydrogens from a copy of `m` and evaluate that.
    std::optional<double> HydrogenSuppressedCopy(Molecule& m) const;

    // `m` must have no explicit Hydrogens. Builds the arrays and evaluates.
    std::optional<double> PolarSurfaceAreaHSuppressed(Molecule& m) const;

    // `m` must have no explicit Hydrogens and the arrays must describe it.
    std::optional<double> PolarSurfaceAreaWithArrays(Molecule& m,
                                           const atomic_number_t* z,
                                           const Atom** atom,
                                           const int* is_aromatic) const;

    double PolarSurfaceAreaInner(Molecule& m, const atomic_number_t* z,
                                 const Atom** atom, const int* is_aromatic) const;

  public:
    NovartisPolarSurfaceArea() = default;

    void set_display_psa_unclassified_atom_messages(int s) {
      _display_psa_unclassified_atom_messages = s;
    }
    void set_return_zero_for_unclassified_atoms(int s) {
      _return_zero_for_unclassified_atoms = s;
    }
    void set_non_zero_contribution_for_SD2(int s) {
      _non_zero_contribution_for_SD2 = s;
    }
    void set_zero_for_all_sulphur_atoms(int s) {
      _zero_for_all_sulphur_atoms = s;
    }
    void set_zero_for_all_phosphorus_atoms(int s) {
      _zero_for_all_phosphorus_atoms = s;
    }
    void set_convert_to_charge_separated(int s) {
      _convert_to_charge_separated = s;
    }

    int display_psa_unclassified_atom_messages() const {
      return _display_psa_unclassified_atom_messages;
    }
    int return_zero_for_unclassified_atoms() const {
      return _return_zero_for_unclassified_atoms;
    }
    int non_zero_contribution_for_SD2() const {
      return _non_zero_contribution_for_SD2;
    }
    int zero_for_all_sulphur_atoms() const {
      return _zero_for_all_sulphur_atoms;
    }
    int zero_for_all_phosphorus_atoms() const {
      return _zero_for_all_phosphorus_atoms;
    }
    int convert_to_charge_separated() const {
      return _convert_to_charge_separated;
    }

    void SetRDKitCompatibility();

    // A function that can process a command line option specifying TPSA
    // calculation options.
    int InitialiseOptions(const Command_Line& cl, char flag);

    // Does NOT modify `m`. The classification assumes the Hydrogen suppressed
    // graph - it reads implicit_hydrogens() and ncon() - so a molecule carrying
    // explicit Hydrogens is copied and stripped first. Without that a primary
    // amine written as N([H])[H] scores 3.24 rather than 26.02, the value for a
    // fully substituted nitrogen. See the commentary in nvrtspsa.cc.
    std::optional<double> PolarSurfaceArea(Molecule& m) const;

    // As above, for a caller that already has these arrays. They must describe
    // `m` as passed. If `m` carries explicit Hydrogens they cannot be used,
    // since removing the Hydrogens renumbers the atoms, and equivalent arrays
    // are built internally.
    std::optional<double> PolarSurfaceArea(Molecule& m, const atomic_number_t* z,
                                           const Atom** atom,
                                           const int* is_aromatic) const;
};

}  // namespace nvrtspsa

#endif  // MOLECULE_TOOLS_NVRTSPSA_H_
