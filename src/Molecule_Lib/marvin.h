#ifndef IW_MARVIN_H
#define IW_MARVIN_H

#include "iwstring.h"
#include "iw_stl_hash_map.h"

class Marvin_Structure_Information
{
  private:
    resizable_array_p<IWString> _atom_colour;
    resizable_array_p<IWString> _bond_colour;

    IW_Hash_Map<atom_number_t, unsigned int> _atom_number_colour_index;
    IW_Hash_Map<int, unsigned int> _bond_number_colour_index;

  public:

    void reset();
    void reset_atoms_and_bonds();

    int debug_print (std::ostream &) const;

    int add_atom_colour (const const_IWSubstring &);
    int add_bond_colour (const const_IWSubstring &);

    void set_atom_number_colour_index (atom_number_t a, unsigned int c)
          { _atom_number_colour_index[a] = c;}
    void set_bond_number_colour_index (atom_number_t b, unsigned int c)
          { _bond_number_colour_index[b] = c;}

//  These functions are used by the writer

    int write_atom_and_bond_colours (std::ostream &) const;

    int atom_colour_specifications_present() const { return _atom_number_colour_index.size();}
    int bond_colour_specifications_present() const { return _bond_number_colour_index.size();}

    int colour_index_for_atom (atom_number_t a) const;
    int colour_index_for_bond (int) const;
};

void set_marvin_structure_information_for_writing (const Marvin_Structure_Information *);

#endif
