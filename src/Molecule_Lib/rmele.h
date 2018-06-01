#ifndef IW_RMELE_H
#define IW_RMELE_H

/*
  When removing elements, we need a means of keeping track of the
  removal.
*/

class Molecule;
class Command_Line;

#include "ematch.h"

class Element_to_Remove : public Element_Matcher
{
  private:
    int _molecules_examined;
    int _molecules_changed;
    int _atoms_removed;

    int _maxcon_to_remove;
    int _add_bond_after_two_connected_removals;

    int _number_to_remove_per_molecule;

//  Private functions

    void _default_values ();
    int  _process (Molecule &);
    int  _process (Molecule &, const int *, int);

  public:
    Element_to_Remove (atomic_number_t);
    Element_to_Remove (const Element *);
    Element_to_Remove (const char *);
    Element_to_Remove (const IWString &);
    ~Element_to_Remove ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int report (std::ostream &) const;

    int maxcon_to_remove () const { return _maxcon_to_remove;}
    void set_maxcon_to_remove (int i) { _maxcon_to_remove = i;}

    int add_bond_after_two_connected_removals () const
         { return _add_bond_after_two_connected_removals;}
    void set_add_bond_after_two_connected_removals (int i)
         { _add_bond_after_two_connected_removals = i;}

    int number_to_remove_per_molecule () const { return _number_to_remove_per_molecule;}
    void set_number_to_remove_per_molecule (int i)
      { _number_to_remove_per_molecule = i;}

    int reset_counters ();

    int process (Molecule &);
    int process (Molecule &, const int *, int);
};

class Elements_to_Remove : public resizable_array_p<Element_to_Remove>
{
  private:
    int _remove_all_non_natural_elements;
    int _remove_all_isotopes;

  public:
    Elements_to_Remove();

    int active () const { return _number_elements;}

    int construct_from_command_line (Command_Line &, int = 0, char = 'X');

    int report (std::ostream &) const;

    int reset_counters ();

    void set_remove_all_non_natural_elements (int i) { _remove_all_non_natural_elements = i;}

    int process (Molecule &);
    int process (Molecule &, const int *, int);
};

#endif
