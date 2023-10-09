#ifndef SMILES_HASH_H
#define SMILES_HASH_H

#include "Foundational/iwstring/iwstring.h"

class Smiles_Hash_Linked_List
{
  private:

    IWString _smiles;

    Smiles_Hash_Linked_List * _next;

//  private functions

    int _add      (const const_IWSubstring &, Smiles_Hash_Linked_List *);

  public:
    Smiles_Hash_Linked_List (const const_IWSubstring &);
    ~Smiles_Hash_Linked_List ();

    IWString & smiles () { return _smiles;}

    Smiles_Hash_Linked_List * next () const { return _next;}
    void set_next (Smiles_Hash_Linked_List * n) { _next = n;}

    int number_smiles () const;

    int contains (const const_IWSubstring &) const;
    int add      (const const_IWSubstring &, Smiles_Hash_Linked_List **);
};

class Smiles_Hash
{
  private:
    int _primary_hash_size;

    Smiles_Hash_Linked_List ** _primary_hash;

  public:
    Smiles_Hash ();
    ~Smiles_Hash ();

    int ok () const;
    int debug_print (std::ostream & os) const;

    int active () const { return _primary_hash_size;}

    int smiles_stored () const;

    int report_hash_statistics (std::ostream &) const;

    int set_primary_hash_size (int ns);

    int contains (const const_IWSubstring &) const;

    int add (const const_IWSubstring &);
};

#endif
