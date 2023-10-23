#ifndef STORAGE_CONDITIONS_H
#define STORAGE_CONDITIONS_H

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/molecule.h"


// If we are storing all variants, there are some common settings shared
// between buildsmidb and in_database.
class AllVariants {
  private:
    int _active;

    // We always store the parent molecule.

    // These will be non empty if we are storing these variants.
    IWString _fragment_modified;
    IWString _chiral_modified;

  public:
    AllVariants();

    int Initialise(Command_Line& cl, char flag);

    int active() const {
      return _active;
    }
    void set_active(int s) {
      _active = s;
    }

    void PrependFragmentModified(Molecule& m) const;
    void PrependChiralModified(Molecule& m) const;
};

/*
  For buildsmidb and in_database, we need an easy way of comparing the
  storage conditions specified on the command line and what might
  already be in the database
*/

class Storage_Conditions
{
  private:
    int _reduce_to_largest_fragment;
    int _remove_chirality;
    int _remove_cis_trans_bonds;
    int _convert_isotopes;

    int _tautomer;
    int _use_aromatic_distinguishing_mf_in_tautomer;
    int _exclude_triple_bonds_from_graph_reduction;

    int _exclude_cc_double_bonds_saturated_from_graph_reduction;
    int _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction;

    // We store the exact parent, the largest fragment and largest fragment
    // without chirality in the same database.
    int _all_variants;

    int _initialised;

    int _good;

    int _ok_to_ignore_mismatched_store_conditions;
    
//  private functions

    int _store_conditions (Db & database) const;
    void _build_store_string (IWString & s) const;
    int  _build_from_string (const const_IWSubstring & s);

    int _ensure_compatible_with_currently_stored(Dbt & fromdb, int verbose);
    int _ensure_consistent_with_current_conditions (Dbt & zdata, int verbose);
    int _initialise_from_these_conditions (const const_IWSubstring & fromdb);
    int _check_compatible_with_currently_stored (const const_IWSubstring & fromdb, const IWString &) const;

  public:
    Storage_Conditions();

    int initialise (const IWString &, Db &);

    int initialised () const { return _initialised;}
    void deactivate () { _initialised = 0;}

    int debug_print(std::ostream &) const;

    int is_consistent () const;

    int initialise (Command_Line & cl, const char flag, Mol2Graph & m2g, const int verbose);

    int store_or_check_currently_stored (Db & db, int verbose);
    int ensure_consistent_with_current_conditions (Db & db, int verbose);

    void set_ok_to_ignore_mismatched_store_conditions (int s) { _ok_to_ignore_mismatched_store_conditions = s;}

    int  reduce_to_largest_fragment    () const { return _reduce_to_largest_fragment;}
    void set_reduce_to_largest_fragment (int s) { _reduce_to_largest_fragment = s; _initialised = 1;}

    int  remove_chirality () const { return _remove_chirality;}
    void set_remove_chirality (int s) { _remove_chirality = s; _initialised = 1;}

    int  remove_cis_trans_bonds () const { return _remove_cis_trans_bonds;}
    void set_remove_cis_trans_bonds (int s) { _remove_cis_trans_bonds = s; _initialised = 1;}

    int  convert_isotopes () const { return _convert_isotopes;}
    void set_convert_isotopes (int s) { _convert_isotopes = s; _initialised = 1;}

    int  tautomer () const { return _tautomer;}
    void set_tautomer (int s) { _tautomer = s; _initialised = 1;}

    int  use_aromatic_distinguishing_mf_in_tautomer () const { return _use_aromatic_distinguishing_mf_in_tautomer;}
    void set_use_aromatic_distinguishing_mf_in_tautomer (int s) { _use_aromatic_distinguishing_mf_in_tautomer = s; _initialised = 1;}

    int all_variants() const {
      return _all_variants;
    }
    void set_all_variants(int s) {
      _all_variants = s;
    }

    int  exclude_triple_bonds_from_graph_reduction () const { return _exclude_triple_bonds_from_graph_reduction;}
    void set_exclude_triple_bonds_from_graph_reduction (int s) { _exclude_triple_bonds_from_graph_reduction = s;}

    void set_exclude_cc_double_bonds_saturated_from_graph_reduction (int s) { _exclude_cc_double_bonds_saturated_from_graph_reduction = s;}
    int  exclude_cc_double_bonds_saturated_from_graph_reduction () const { return _exclude_cc_double_bonds_saturated_from_graph_reduction;}

    void set_exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction (int s) { _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction = s;}
    int  exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction () const { return _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction;}

    int form_key (Molecule & m, const Mol2Graph & mol2graph, int include_chirality, IWString & key) const;

    int transfer_to_mol2graph(Mol2Graph & mol2graph) const;
    void transfer_to_all_variants(AllVariants & all_variants) const {
      if (_all_variants) {
        all_variants.set_active(1);
      }
    }
};

#endif
