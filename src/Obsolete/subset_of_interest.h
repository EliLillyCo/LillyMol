#ifndef MOLECULE_TOOLS_SUBSET_OF_INTEREST_H
#define MOLECULE_TOOLS_SUBSET_OF_INTEREST_H

class Molecule;

namespace subset_of_interest {

// Originally developed to support Qupkake.
// We have a molecule with an atom of interest marked.
// Identify bonds to break and reduce to the fragment
// containing the atom of interest.

class SubsetOfInterest {
  private:
    // if set, we apply this isotope to each break point.
    isotope_t _isotope_at_attachment_points;

    // By default, we do not break bonds to aromatic rings.
    // If set, this value will be added to the isotopic label
    // to the atom at the other end of the bond being broken.
    // Note that this is done after _isotope_at_attachment_points
    int _break_bonds_to_aromatic_rings;

  public:
    SubsetOfInterest();

    void set_isotope_at_attachment_points(isotope_t s) {
      _isotope_at_attachment_points = s;
    }

    void set_break_bonds_to_aromatic_rings(int s) {
      _break_bonds_to_aromatic_rings = s;
    }

    int ReduceToFragmentOfInterest(Molecule& m, atom_number_t zatom);
};

};  // namespace subset_of_interest

#endif  // MOLECULE_TOOLS_SUBSET_OF_INTEREST_H
