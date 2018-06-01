#ifndef E_MATCH_H
#define E_MATCH_H

#include "iwaray.h"
#include "iwcrex.h"

#include "element.h"

/*
  Atom matching capability. 
*/

class Element_Matcher
{
  private:
    const Element * _e;

    int _isotope;

    int _match_organic_only;
    int _match_non_organic_only;
    int _match_non_periodic_only;

    IW_Regular_Expression _symbol_rx;

//  private functions

    void _default_values ();

  public:
    Element_Matcher ();
    Element_Matcher (atomic_number_t);
    Element_Matcher (const Element *);
    Element_Matcher (const char *);
    Element_Matcher (const IWString &);

    int ok () const;
    int debug_print (std::ostream &) const;

    void set_element (const Element *);
    int  construct_from_string (const char *, int);
    int  construct_from_string (const char *);
    int  construct_from_string (const const_IWSubstring &);
    int  construct_from_string (const IWString &);

    int operator_less_less (std::ostream & os) const;

    const Element * element () const { return _e;}

    int isotope () const { return _isotope;}

    int matches (const Element *, int = 0);    // = 0 parameter is isotope
};

extern std::ostream &
operator<< (std::ostream &, const Element_Matcher &);

class Command_Line;

class Set_of_Element_Matches : public resizable_array_p<Element_Matcher>
{
  private:
  public:

    int construct_from_command_line (Command_Line &, int, char);

    int matches (const Element *, int = 0);    // = 0 parameter is isotope
};

extern void display_element_matcher_syntax (std::ostream & os);

#endif
