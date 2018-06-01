#ifndef ATOM_ALIAS_H
#define ATOM_ALIAS_H

/*
  From an ISIS reaction file we may have an atom alias 
*/

#include "iwmtypes.h"
#include "iwstring.h"
#include "iwstring_data_source.h"

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

#endif
