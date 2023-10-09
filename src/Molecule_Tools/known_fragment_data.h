#ifndef Known_Fragment_Data_H
#define Known_Fragment_Data_H

/*
  We have a list of known parent molecules and known counterions.
*/

#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/molecule.h"

/*
  Data about an individual fragment.
*/

class Known_Fragment_Data
{
  private:
    int _only_check_molecular_formula;    // by default, we check the unique smiles

    int _remove_everything_if_all_fragments_match;   // if every fragment is a SALT

    extending_resizable_array<int> _natoms;

    typedef IWString_STL_Hash_Set _formula_usmi;

    _formula_usmi _known_salt_mf;
    _formula_usmi _known_salt_usmi;
    _formula_usmi _known_parent_mf;
    _formula_usmi _known_parent_usmi;

//  private functions

    int _common_read  (const const_IWSubstring & fname, _formula_usmi & mf, _formula_usmi & usmi);
    int _common_read  (data_source_and_type<Molecule> & input, _formula_usmi & mf, _formula_usmi & usmi);
    int _add_to_hash  (Molecule & m, _formula_usmi & mf, _formula_usmi & usmi);

    int _common_print_set (const _formula_usmi & h, std::ostream & os) const;

    int _scan_known_fragments (Molecule & m,
                                            const IWString & molecular_formula,
                                            resizable_array_p<Molecule> & fragments,
                                            _formula_usmi & mf,
                                            _formula_usmi & usmi,
                                            resizable_array<int> & fragments_identified) const;

    int _remove_soap (Molecule & m) const;

    int _remove_non_organic (Molecule & m) const;

    int _delete_set_of_fragments (Molecule & m, const resizable_array<int> & fragments_to_be_removed) const;

  public:
    Known_Fragment_Data ();
    ~Known_Fragment_Data ();

    int debug_print (std::ostream & os) const;

    int active () const { return (_known_salt_mf.size () > 0) || (_known_parent_mf.size () > 0);}
    void deactivate();

    void set_only_check_molecular_formula (int s) { _only_check_molecular_formula = s;}
    void set_remove_everything_if_all_fragments_match (int s) { _remove_everything_if_all_fragments_match = s;}

    int read_known_salts (const const_IWSubstring &);
    int read_known_parents (const const_IWSubstring &);

    int process (Molecule & m);
};

#endif
