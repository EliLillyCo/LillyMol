#ifndef MOLECULE_LIB_ROTBOND_COMMON_H_
#define MOLECULE_LIB_ROTBOND_COMMON_H_

class Command_Line;

extern int part_of_otherwise_non_rotabable_entity (Molecule & m,
                                        atom_number_t a1,
                                        atom_number_t a2);

extern int triple_bond_at_either_end (const Molecule & m,
                           const Bond * b);

extern int is_non_rotatable_amide (Molecule & m,
                        atom_number_t a1,
                        atom_number_t a2);

extern int is_non_rotatable_sulphonamide (Molecule & m,
                               atom_number_t zatom1,
                               atom_number_t zatom2);
namespace quick_rotbond {
// A class that can be used for computing rotatable bonds.
// The cost/completeness of the calculation can be governed
// by settings.

class QuickRotatableBonds {
  public:

  enum class RotBond {
    kUndefined,
    kQuick,
    kExpensive
  };

  private:

    RotBond _calculation;

  // private functions.

    // The quickest calculation just counts the number of non ring
    // single bonds that are non terminal.

    int Quickest(Molecule& m);
    // the more expensive calculation handles CF3 and amides.
    int Expensive(Molecule& m, int* bond_rotatable);

    // If requested, we can label the rotatable bonds.
    isotope_t _isotope;

    // Even though we can label atoms as being part of a rotatable
    // bond, this leaves ambiguities, since if the same isotope is
    // used for every rotatable bond, we can have bonds that have
    // an isotope at each end, but the bond is not rotatable.
    // For that reason, we can optionally assign a unique isotope
    // to each rotatable bond.
    int _unique_isotope_each_bond = 0;

  public:
    QuickRotatableBonds();

    int Initialise(Command_Line& cl, char flag);

    int active () const {
      return _calculation != RotBond::kUndefined;
    }

    void set_isotope(isotope_t iso) {
      _isotope = iso;
    }
    void set_unique_isotope_each_bond(int s) {
      _unique_isotope_each_bond = s;
    }

    void set_calculation_type(RotBond rbtype) {
      _calculation = rbtype;
    }

    // Using whatever method is specified in _calculation, return the
    // number of rotatable bonds.
    // If bond_rotatable is not null each bond, set
    // the corresponding entry in `bond_rotatable` to 1 for rotatable bonds.
    int Process(Molecule& m, int* bond_rotatable = nullptr);

    // Use the current computation to determine the number of rotatable bonds
    // between atoms `a1` and `a2`. A greedy path is traced between the two
    // atoms and the first traversal from `a1` to `a2` is used.
    // Note that this could produce unstable results if there is a case of
    // multiple paths with different rotatable bonds, but not sure that
    // can happen.
    int RotatableBondsBetween(Molecule& m,
                atom_number_t a1,
                atom_number_t a2);

    // return an array of m.natoms()*m.natoms() with the rotatable bonds
    // between atoms.
    std::unique_ptr<int[]> RotatableBondsBetween(Molecule& m);
};

}  // namespace quick_rotbond

#endif  // MOLECULE_LIB_ROTBOND_COMMON_H_
