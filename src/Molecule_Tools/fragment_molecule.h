#ifndef MOLECULE_TOOLS_FRAGMENT_MOLECULE_H_
#define MOLECULE_TOOLS_FRAGMENT_MOLECULE_H_

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

// A common operation is to identify a set of bonds in a molecule that are to be broken.
// Generally this is based on the rules in dicer, and long term, this should be incorproated
// into dicer.
namespace fragment_molecule {

class MoleculeFragmenter {
  private:
    int _add_queries_to_default_rules;

    // Do we break bonds between fully saturated carbon atoms.
    int _break_carbon_carbon_bonds;

    resizable_array_p<Substructure_Query> _to_break;
    resizable_array_p<Substructure_Query> _never_break;

  // private functions.
    int IdentifyBreakableBonds(Molecule& m, int* status);
    int IdentifyBreakableBond(Molecule& m, 
                        atom_number_t a1, atom_number_t a2,
                        int* status);

  public:
    MoleculeFragmenter();

    // Command line options -q and -Q are used.
    int Initialise(Command_Line& cl);

    void set_add_queries_to_default_rules(int s) {
      _add_queries_to_default_rules = s;
    }
    void set_break_carbon_carbon_bonds(int s) {
      _break_carbon_carbon_bonds = s;
    }

    // Return an array of bond numbers that can be broken. It will be up
    // to the caller to retrieve those bonds from the Molecule. Note that
    // we return a sorted list of bond numbers, highest to lowest, so
    // these can be processed in the order of the array.
    int IdentifyBreakableBonds(Molecule& m, resizable_array<int>& bonds_to_break);
};

}  // namespace fragment_molecule

#endif  // MOLECULE_TOOLS_FRAGMENT_MOLECULE_H_
