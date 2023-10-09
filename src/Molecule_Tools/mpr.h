#ifndef MOLECULE_LIB_MPR_H_
#define MOLECULE_LIB_MPR_H_

#define NPROPERTIES 8

#include "Molecule_Lib/molecule.h"

// Aug 2023. A more flexible way of selecting what gets computed, also
// extending the range of what can be computed. Should get better
// integration with iwdescr...

namespace mpr {
enum Feature {
  kNatoms = 0,
  kLargestRingSize = 1,
  kNrings = 2,
  kRingAtoms = 3,
  kAromaticAtoms = 4,
  kFusedRingAtoms = 5,
  kHeteroatomCount = 6,
  kUnsaturation = 7,
  kMaxDist = 8,
  kScaffoldAtoms = 9,
  kLinkerAtoms = 10,
  kLargestRingSys = 11,
  kRotBonds = 12
};

int GenerateMolecularProperties(Molecule& m,
                const resizable_array<Feature>& features,
                int* result);

class MolecularPropertiesGenerator {
  private:
    resizable_array<Feature> _features;

  // private functions.

    // If no features specified, do these.
    void AddDefaultFeatures();

  public:
    int Initialise(Command_Line& cl, char flag);

    int number_features() const {
      return _features.number_elements();
    }

    // Results can be generated as any int form.
    template <typename T>
    int GenerateMolecularProperties(Molecule& m, T* result);
};

}  // namespace mpr

class Molecular_Properties_Generator
{
  private:
    int _unsaturation_includes_aromatic;

// private functions

    template <typename T> int _molecular_properties_generation (Molecule & m, const int * ring_membership, T * properties) const;

  public:
    Molecular_Properties_Generator ();
    
    void set_unsaturation_includes_aromatic (int s) { _unsaturation_includes_aromatic = s;}

//  int operator () (Molecule &, unsigned char[]) const;
    template <typename T> int operator () (Molecule &, T[]) const;
};

#endif  // MOLECULE_LIB_MPR_H_
