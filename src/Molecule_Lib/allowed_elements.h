#ifndef ALLOWED_ELEMENTS_H
#define ALLOWED_ELEMENTS_H

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
//  Allowed_Elements();

    int build_from_command_line (Command_Line & c, char flag, int verbose);

    int contains_non_allowed_atoms (const Molecule & m) const;

    void set_allow (atomic_number_t, int);

    void reset_to_defaults();
};

#endif
