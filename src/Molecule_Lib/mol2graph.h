#ifndef MOLECULE_LIB_MOL2GERAPH_H_
#define MOLECULE_LIB_MOL2GERAPH_H_

#include "Foundational/cmdline/cmdline.h"
/*
  There are lots of choices when we produce a molecular graph
  In the code there are two means of doing this for historical reasons.
  This struct
*/

class Mol2Graph
{
  private:
    int _exclude_triple_bonds_from_graph_reduction;
    int _revert_all_directional_bonds_to_non_directional;
    int _preserve_cc_double_bonds_saturated;
    int _preserve_cc_double_bonds_no_heteroatoms;
    int _remove_chiral_centres;
    
    int _append_molecular_formula;

    int _aromatic_distinguishing_formula;

    bool _active;

  public:
    Mol2Graph();

    int construct (Command_Line & cl, const char flag, const int verbose);

    int debug_print (std::ostream &) const;

    bool active() const { return _active;}
    void set_active(bool s) { _active = s;}

    int exclude_triple_bonds_from_graph_reduction () const { return _exclude_triple_bonds_from_graph_reduction;}
    int revert_all_directional_bonds_to_non_directional() const { return _revert_all_directional_bonds_to_non_directional;}
    int preserve_cc_double_bonds_saturated () const { return _preserve_cc_double_bonds_saturated;}
    int preserve_cc_double_bonds_no_heteroatoms () const { return _preserve_cc_double_bonds_no_heteroatoms;}
    int append_molecular_formula () const { return _append_molecular_formula;}
    int aromatic_distinguishing_formula () const {return _aromatic_distinguishing_formula;}
    int remove_chiral_centres () const { return _remove_chiral_centres;}

    int some_kind_of_double_bond_preservation_active () const { return _preserve_cc_double_bonds_saturated + _preserve_cc_double_bonds_no_heteroatoms;}

    void set_exclude_triple_bonds_from_graph_reduction(int s) { _exclude_triple_bonds_from_graph_reduction = s;}
    void set_revert_all_directional_bonds_to_non_directional(int s) { _revert_all_directional_bonds_to_non_directional = s;}
    void set_preserve_cc_double_bonds_no_heteroatoms (int s) { _preserve_cc_double_bonds_no_heteroatoms = s;}
    void set_preserve_cc_double_bonds_saturated (int s) { _preserve_cc_double_bonds_saturated = s;}
    void set_append_molecular_formula (int s) { _append_molecular_formula = s;}
    void set_aromatic_distinguishing_formula(int s) { _aromatic_distinguishing_formula = s;}
    void set_remove_chiral_centres (int s) { _remove_chiral_centres = s;}

    // Activate what are likely the most useful options. May change over time.
    void TurnOnMostUsefulOptions();
};

// Do not use, prefer Mol2Graph setting.
extern void set_exclude_triple_bonds_from_graph_reduction (int);

#endif  // MOLECULE_LIB_MOL2GERAPH_H_
