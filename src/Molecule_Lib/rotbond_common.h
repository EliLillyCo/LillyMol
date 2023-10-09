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
    int Expensive(Molecule& m);

    // If requested, we can label the rotatable bonds.
    isotope_t _isotope;

  public:
    QuickRotatableBonds();

    int Initialise(Command_Line& cl, char flag);

    int active () const {
      return _calculation != RotBond::kUndefined;
    }

    void set_calculation_type(RotBond rbtype) {
      _calculation = rbtype;
    }

    // Using whatever method is specified in _calculation, return the
    // number of rotatable bonds.
    int Process(Molecule& m);
};

}  // namespace quick_rotbond

#endif  // MOLECULE_LIB_ROTBOND_COMMON_H_
