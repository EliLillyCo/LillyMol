#ifndef IW_NUMBER_ASSIGNER_H
#define IW_NUMBER_ASSIGNER_H

#include <iostream>
#include <fstream>

class Molecule;
class Command_Line;

class Number_Assigner
{
  private:
    IWString _prefix_string;
    int _next_number_to_assign;

//  Some people don't want parentheses

    int _include_parentheses;

//  We can make fixed width names

    int _number_digits;

//  Optionally we can just replace the name with our unique form

    int _replace_name;

//  or we might want to replace just the first token of an existing name

    int _replace_first_token;

//  Aug 2014. We might want to merge with the name R(2)_existing_name

    IWString _separator_to_existing;

//  Often it is useful to create a cross reference file

    std::ofstream _cross_reference_file;

    int _only_apply_if_existing_name_is_empty;

  public:
    Number_Assigner ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int initialise (Command_Line &, char = 'n', int = 0);
    int initialise (int);

    int active () const { return _next_number_to_assign >= 0;}

    void deactivate() { _next_number_to_assign = -5;}

    int process (Molecule &);

    int process (IWString &);
};

int display_standard_number_assigner_options (std::ostream &, char = 'n');

#endif
