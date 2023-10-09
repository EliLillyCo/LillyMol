#ifndef MOLECULE_LIB_ATOM_ALIAS_H_
#define MOLECULE_LIB_ATOM_ALIAS_H_

/*
  From an ISIS reaction file we may have an atom alias 
*/

#include "Foundational/iwstring/iwstring.h"
#include "iwmtypes.h"

class Atom_Alias
{
  private:
    atom_number_t _atom;
    IWString      _alias;

//  private functions

    void _copy (const Atom_Alias &);

  public:
    Atom_Alias ();
    Atom_Alias (const Atom_Alias &);

    Atom_Alias & operator= (const Atom_Alias &);

    template <typename T> int build (const const_IWSubstring &, T &);

    atom_number_t atom_number () const { return _atom;}
    void set_atom_number (atom_number_t s) { _atom = s;}
    const IWString & alias () const { return _alias;}
};

#endif  // MOLECULE_LIB_ATOM_ALIAS_H_
