#ifndef MOLECULE_TOOLS_NVRTSPSA_H_
#define MOLECULE_TOOLS_NVRTSPSA_H_

#include <optional>

#include "Molecule_Lib/molecule.h"

extern double novartis_polar_surface_area (Molecule & m);

extern double novartis_polar_surface_area (Molecule & m, 
                             const atomic_number_t * z,
                             const Atom ** atom,
                             const int * is_aromatic);
extern void   set_display_psa_unclassified_atom_mesages(int s);
extern void   set_return_zero_for_unclassified_atoms (int s);
// default is 1
extern void   set_non_zero_constribution_for_SD2(int s);

namespace nvrtspsa {

// Class for performing Novartis Polar Surface Area calculations.
// Note, currently this code just calls the standalone functions above.
// Calling any of the set_ methods also updates the file scope static
// settings. Not ideal, the code needs to be restructured.
class NovartisPolarSurfaceArea {
  private:
    int _display_psa_unclassified_atom_mesages;
    int _return_zero_for_unclassified_atoms;
    // default is 1
    int _non_zero_constribution_for_SD2;

  public:
    NovartisPolarSurfaceArea();

    void set_display_psa_unclassified_atom_mesages(int s) {
      _display_psa_unclassified_atom_mesages = s;
      ::set_display_psa_unclassified_atom_mesages(s);
    }
    void set_return_zero_for_unclassified_atoms(int s) {
      _return_zero_for_unclassified_atoms = s;
      ::set_return_zero_for_unclassified_atoms(s);
    }
    void set_non_zero_constribution_for_SD2(int s) {
      _non_zero_constribution_for_SD2 = s;
      ::set_non_zero_constribution_for_SD2(s);
    }

    // Long term, this should return nullopt if there is an unclassified
    // atom. But currently there is no clean way to do that, the entire
    // code base in the corresponding .cc file would need restructuring.
    std::optional<double> PolarSurfaceArea(Molecule& m);
};

}  // namespace nvrtspsa

#endif  // MOLECULE_TOOLS_NVRTSPSA_H_
