#ifndef MOLECULE_LIB_ALLOWED_ELEMENTS_H_
#define MOLECULE_LIB_ALLOWED_ELEMENTS_H_

#include "element.h"

class Molecule;
class Command_Line;

class Allowed_Elements
{
  private:
    atomic_number_t _allowed_element[HIGHEST_ATOMIC_NUMBER + 1];

//  private functions

    void _default_values();

  public:
    Allowed_Elements();

    int build_from_command_line (Command_Line & c, char flag, int verbose);

    int contains_non_allowed_atoms (const Molecule & m) const;

    void set_allow (atomic_number_t, int);

    // Any element that is a metal will be set to exclude.
    void exclude_metals();

    void reset_to_defaults();

    // Set all elements to disallowed.
    void Clear();
};

#endif  // MOLECULE_LIB_ALLOWED_ELEMENTS_H_
