#ifndef IW_ELEMENT_TRANS_H
#define IW_ELEMENT_TRANS_H

#include <iostream>

#include "iwaray.h"

class Element;
class Molecule;
class IWString;
class Molecule_to_Match;

#include "ematch.h"

class Element_Transformation
{
  private:
    int   _transform_every_atom_type;
    Element_Matcher _from;
    const Element * _to;
    int _isotope;

    int _molecules_processed;
    int _molecules_changed;
    int _atoms_changed;

//  private functions

    void _default_values ();

  public:
    Element_Transformation ();

    int ok () const;
    int debug_print (std::ostream &) const;

//  int build (const char *);
    int build (const IWString &);

    int process (Molecule &);

    int process (Molecule_to_Match &);
};

class Element_Transformations : public resizable_array_p<Element_Transformation>
{
  private:
  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int active () const { return _number_elements;}

    int construct_from_command_line (Command_Line &, int = 0, char = 't');

    int process (Molecule *);

    int process (Molecule &);

    int process (Molecule_to_Match &);
};

class Command_Line;

extern int display_standard_etrans_options (std::ostream &, char = 't');

extern int process_element_transformations (Command_Line &,
                                            Element_Transformations &,
                                            int = 0,
                                            char = 't');
#endif
