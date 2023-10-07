#ifndef TT_SUPPORT_H
#define TT_SUPPORT_H

#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/molecule.h"

class TopoTorsion : public IWString
{
  private:
    int _count;

  public:
    TopoTorsion ();
    TopoTorsion (const char *);

    int  count () const { return _count;};
    unsigned char count_as_unsigned_char () const;
    void increment_count () { _count++;}
};

extern int sort_and_count_tts (resizable_array_p<TopoTorsion> & tt);

extern int initialise_single_character_atom_symbols (const Molecule & m, IWString & asymbol);

extern void
form_canonical_order (Molecule & m,
                      const atomic_number_t * z,
                      const int * ncon,
                      atom_number_t & a0,
                      atom_number_t & a1,
                      atom_number_t & a2,
                      atom_number_t & a3);

extern void
form_canonical_order (const int * invariant,
                      atom_number_t & a0,
                      atom_number_t & a1,
                      atom_number_t & a2,
                      atom_number_t & a3);

void
append_numbers_eight_chars (IWString & tt,
             const Molecule & m,
             const IWString & asymbol,
             atom_number_t a,
             int pi,
             int ncon);
void
append_numbers_twelve_chars (IWString & tt,
             const Molecule & m,
             atom_number_t a,
             int pi,
             int ncon);

void
append_numeric_torsion_names (IWString & tt,
             const Molecule & m,
             atom_number_t a,
             int pi,
             int ncon);

#endif
