#include <stdlib.h>
#include <iostream>

#include "db_cxx.h"

#include "Foundational/cmdline/cmdline.h"

#include "storage_conditions.h"

static IWString store_info_key("_STORE_INFO");

using std::cerr;
using std::endl;

Storage_Conditions::Storage_Conditions()
{
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _remove_cis_trans_bonds = 0;
  _convert_isotopes = 0;

  _tautomer = 0;
  _use_aromatic_distinguishing_mf_in_tautomer = 0;

  _exclude_triple_bonds_from_graph_reduction = 0;

  _exclude_cc_double_bonds_saturated_from_graph_reduction = 0;
  _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction = 0;

  _all_variants = 0;

  _initialised = 0;

  _good = 1;

  _ok_to_ignore_mismatched_store_conditions = 0;

  return;
}

int
Storage_Conditions::debug_print(std::ostream & output) const
{
  output << "Storage_Conditions::debug_print:_initialised? " << _initialised << endl;
  output << _reduce_to_largest_fragment << " _reduce_to_largest_fragment\n";
  output << _remove_chirality << " _remove_chirality\n";
  output << _remove_cis_trans_bonds << " _remove_cis_trans_bonds\n";
  output << _convert_isotopes << " _convert_isotopes\n";
  output << _all_variants << " _all_variants\n";
  output << _tautomer << " _tautomer\n";
  output << _use_aromatic_distinguishing_mf_in_tautomer << " _use_aromatic_distinguishing_mf_in_tautomer\n";
  output << _exclude_triple_bonds_from_graph_reduction << " _exclude_triple_bonds_from_graph_reduction\n";
  output << _exclude_cc_double_bonds_saturated_from_graph_reduction << " _exclude_cc_double_bonds_saturated_from_graph_reduction\n";
  output << _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction << " _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction\n";

  return 1;
}

int
Storage_Conditions::initialise (const IWString & k,
                                Db & database)
{
  Dbt dkey ((char *)(k.rawchars()), k.length());

  Dbt zdata;

  if (0 != database.get(NULL, &dkey, &zdata, 0))
  {
    _initialised = 0;
    return 0;
  }

  const_IWSubstring fromdb((const char *)zdata.get_data(), zdata.get_size());

  if (! _build_from_string(fromdb))
    return 0;

  _initialised = 1;

  return 1;
}

void
display_storage_conditions_options (const char flag,
                                    std::ostream & output)
{
  output << "Storage conditions for buildsmidb and in_database\n";
  output << " -" << flag << " l         reduce to largest fragment\n";
  output << " -" << flag << " c         exclude chirality\n";
  output << " -" << flag << " z         exclude cis trans bonds\n";
  output << " -" << flag << " a         graph representation (tautomers)\n";
  output << " -" << flag << " q         exclude triple bonds graph reduction\n";
  output << " -" << flag << " d         exclude CC double bonds near saturated atoms from graph reduction\n";
  output << " -" << flag << " h         exclude CC double bonds no heteroatoms from graph reduction\n";
  output << " -" << flag << " w         use aromatic distinguishing molecular formula in graph reduction\n";
  output << " -" << flag << " i         discard isotopes\n";
  output << " -" << flag << " allg      use all graph specific options (a q d h w i)\n";

  return;
}

int
Storage_Conditions::initialise (Command_Line & cl,
                                const char flag,
                                Mol2Graph & mol2graph,
                                const int verbose)
{
  const_IWSubstring s;
  for (int i = 0; cl.value(flag, s, i); ++i)
  {
    if ('s' == s)
    {
      set_exclude_cc_double_bonds_saturated_from_graph_reduction(1);
      mol2graph.set_preserve_cc_double_bonds_saturated(1);

      if (verbose)
        cerr << "Will not convert C=C bonds adjacent to fully saturated Carbons in tautomer computations\n";
    }
    else if ('d' == s)
    {
      set_exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction(1);
      mol2graph.set_preserve_cc_double_bonds_no_heteroatoms(1);

      if (verbose)
        cerr << "Will not convert C=C bonds adjacent to all Carbon atoms in tautomer computations\n";
    }
    else if ('q' == s)
    {
      set_exclude_triple_bonds_from_graph_reduction(1);

      mol2graph.set_exclude_triple_bonds_from_graph_reduction(1);

      if (verbose)
        cerr << "Will not convert triple bonds in tautomer computations\n";
    }
    else if ('w' == s)
    {
      set_use_aromatic_distinguishing_mf_in_tautomer(1);

      if (verbose)
        cerr << "Will use formula_distinguishing_aromatic() for tautomer mf\n";
    }
    else if ("allg" == s)
    {
      set_tautomer(1);
      set_exclude_cc_double_bonds_saturated_from_graph_reduction(1);
      mol2graph.set_preserve_cc_double_bonds_saturated(1);
      set_exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction(1);
      mol2graph.set_preserve_cc_double_bonds_no_heteroatoms(1);
      set_exclude_triple_bonds_from_graph_reduction(1);
      mol2graph.set_exclude_triple_bonds_from_graph_reduction(1);
      set_use_aromatic_distinguishing_mf_in_tautomer(1);
      set_convert_isotopes(1);
      set_remove_cis_trans_bonds(1);
      set_remove_chirality(1);
//    mol2graph.debug_print(cerr);

      if (verbose)
        cerr << "All graph specific options turned on\n";
    }
    else if ('a' == s)
    {
      set_tautomer(1);
      if (verbose)
        cerr << "Will store molecular graphs\n";
    }
    else if ('i' == s)
    {
      set_convert_isotopes(1);
      if (verbose)
        cerr << "Will convert isotopes to non isotopic forms\n";
    }
    else if ('z' == s)
    {
      set_remove_cis_trans_bonds(1);

      if (verbose)
        cerr << "Cis trans bonding information suppressed\n";
    }
    else if ('l' == s)
    {
      set_reduce_to_largest_fragment(1);
      if (verbose)
        cerr << "Will strip to largest fragment before doing lookup\n";
    }
    else if ('c' == s)
    {
      set_remove_chirality(1);

      if (verbose)
        cerr << "Chirality information excluded from smiles\n";
    }
    else if ("help" == s)
    {
      display_storage_conditions_options(flag, cerr);
      exit(1);
    }
    else
    {
      cerr << "Storage_Conditions::initialise:unrecognised directive '" << s << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Storage_Conditions::store_or_check_currently_stored (Db & database,
                                                     int verbose)
{
  Dbt dkey((void *)(store_info_key.rawchars()), store_info_key.length());
  Dbt zdata;

  int rc = database.get(NULL, &dkey, &zdata, 0);

  if (0 == rc)
    return _ensure_compatible_with_currently_stored(zdata, verbose);

  if (verbose)
    cerr << "Stored info absent, storing...\n";

// Conditions not stored

  return _store_conditions(database);
}

void
Storage_Conditions::_build_store_string (IWString & s) const
{
  s.resize_keep_storage(0);

  if (_reduce_to_largest_fragment)
    s << 'L';
  if (_remove_chirality)
    s << 'C';
  if (_remove_cis_trans_bonds)
    s << 'Z';
  if (_convert_isotopes)
    s << 'I';
  if (_all_variants) {
    s << 'V';
  }
  if (_tautomer)
    s << 'A';
  if (_use_aromatic_distinguishing_mf_in_tautomer)
    s << 'W';
  if (_exclude_triple_bonds_from_graph_reduction)
    s << 'T';
  if (_exclude_cc_double_bonds_saturated_from_graph_reduction)
    s << 'S';
  if (_exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction)
    s << 'H';

  return;
}

int
Storage_Conditions::_store_conditions (Db & database) const
{
  IWString to_store;

  _build_store_string(to_store);

  Dbt dkey((char *)(store_info_key.rawchars()), store_info_key.length());
  Dbt zdata((char *)(to_store.rawchars()), to_store.length());

//cerr << "Storing '" << to_store << "'\n";

  int rc = database.put(NULL, &dkey, &zdata, 0);

  if (0 == rc)
    return 1;

  cerr << "Storage_Conditions::_store_conditions:cannot store key\n";
  database.err(rc, "");

  return 0;
}

int
Storage_Conditions::_ensure_compatible_with_currently_stored (Dbt & zdata,
                                                              int verbose)
{
  const_IWSubstring fromdb((const char *) zdata.get_data(), zdata.get_size());

  if (verbose > 1)
    cerr << "Currently stored conditions '" << fromdb << "'\n";

  IWString my_conditions;

  _build_store_string(my_conditions);

  if (my_conditions == fromdb)   // great, finished
    return 1;

  if (_initialised)    // some things came in from command line
    return _check_compatible_with_currently_stored(fromdb, my_conditions);
  else
    return _initialise_from_these_conditions(fromdb);
}

/*
  We had been initialised from the command line, are we compatible with what is
  in the database? Just compare strings
*/

int
Storage_Conditions::_check_compatible_with_currently_stored (const const_IWSubstring & fromdb,
                                                             const IWString & my_conditions) const
{
  if (fromdb == my_conditions)
    return 1;

  cerr << "Incompatible storage conditions, stored '" << fromdb << "', proposed '" << my_conditions << "'";

  if (_ok_to_ignore_mismatched_store_conditions)
    cerr << ", ignored\n";
  else
    cerr << endl;

  return _ok_to_ignore_mismatched_store_conditions;
}

int
Storage_Conditions::_build_from_string (const const_IWSubstring & s)
{
  int n = s.length();

  for (int i = 0; i < n; i++)
  {
    char c = s[i];

    if ('L' == c)
      _reduce_to_largest_fragment = 1;
    else if ('C' == c)
      _remove_chirality = 1;
    else if ('Z' == c)
      _remove_cis_trans_bonds = 1;
    else if ('I' == c)
      _convert_isotopes = 1;
    else if (c == 'V')
      _all_variants = 1;
    else if ('A' == c)
      _tautomer = 1;
    else if ('W' == c)
      _use_aromatic_distinguishing_mf_in_tautomer = 1;
    else if ('T' == c)
      _exclude_triple_bonds_from_graph_reduction = 1;
    else if ('S' == c)
      _exclude_cc_double_bonds_saturated_from_graph_reduction = 1;
    else if ('H' == c)
      _exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction = 1;
    else
    {
      cerr << "Storage_Conditions::_ensure_compatible_with_currently_stored:unrecognised attribute '" << s << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  Nothing on the command line, initialise ourselves
  with what came from the database
*/

int
Storage_Conditions::_initialise_from_these_conditions (const const_IWSubstring & fromdb)
{
  if (! _build_from_string(fromdb))
    return 0;

  return is_consistent();
}

int
Storage_Conditions::is_consistent () const
{
  if (_remove_chirality && _tautomer)
    return 0;

  if (_use_aromatic_distinguishing_mf_in_tautomer && 0 == _tautomer)
    return 0;

  if (_exclude_triple_bonds_from_graph_reduction && 0 == _tautomer)
    return 0;

  return 1;
}

/*
  We are reading a database and need to know whether anything that has
  come from the command line is compatible with what we have.
  If we have not been initialised, we initialise ourselves from
  what is stored
*/

int
Storage_Conditions::ensure_consistent_with_current_conditions (Db & db, int verbose)
{
  Dbt dkey((void *)(store_info_key.rawchars()), store_info_key.length());
  Dbt zdata;

  if (0 == db.get(NULL, &dkey, &zdata, 0))   // something stored, check it
    return _ensure_consistent_with_current_conditions(zdata, verbose);

// There is nothing stored in the database, we are by definition correct.
 
  return 1;
}

/*
  Storage conditions are stored in the database.
*/

int
Storage_Conditions::_ensure_consistent_with_current_conditions (Dbt & zdata,
                                                        int verbose)
{
  const_IWSubstring fromdb ((const char *)(zdata.get_data()), zdata.get_size());

  if (_initialised)
  {
    IWString my_conditions;
    _build_store_string(my_conditions);

    if (fromdb == my_conditions)
      return 1;

    cerr << "Incompatible storage conditions, specified '" << my_conditions << "', stored '" << fromdb << "'\n";
    return 0;
  }

// We are not initialised. Build ourselves from what was stored

  return _initialise_from_these_conditions(fromdb);
}

/*
  Slight problem, the include_chirality argument over-rides what is stored. Probably OK
*/

int
Storage_Conditions::form_key (Molecule & m,
                              const Mol2Graph & mol2graph,
                              int include_chirality,
                              IWString & key) const
{
//mol2graph.debug_print(cerr);

//cerr << "On entry " << m.unique_smiles() << endl;

  if ((_reduce_to_largest_fragment && m.number_fragments() > 1) ||
      _tautomer ||
      (_convert_isotopes && m.number_isotopic_atoms()) ||
      (0 == include_chirality && m.chiral_centres()) ||
      (_remove_cis_trans_bonds && m.cis_trans_bonds_present()) )
  {
    Molecule tmp(m);      // make a copy

    if (_reduce_to_largest_fragment)
      tmp.reduce_to_largest_fragment_carefully();   // note non-standard use - OK

    if (_tautomer || _convert_isotopes)    // the graph doesn't have isotopes in it
      tmp.transform_to_non_isotopic_form();

    if (! include_chirality)
      tmp.remove_all_chiral_centres();

    if (_remove_cis_trans_bonds)
      tmp.revert_all_directional_bonds_to_non_directional();

//  cerr << "Smiles for lookup " << tmp.unique_smiles() << " tautomer? " << _tautomer << endl;

    if (_tautomer)
    {
      IWString f;
      if (_use_aromatic_distinguishing_mf_in_tautomer)
        tmp.formula_distinguishing_aromatic(f);
      else
        tmp.molecular_formula(f);

      tmp.change_to_graph_form(mol2graph);

      key = tmp.unique_smiles();
      key << ':' << f;
    }
    else
      key = tmp.unique_smiles();

//  cerr << "After adding graph stuff " << key << endl;
  }
  else
  {
    key = m.unique_smiles();
//  cerr << "Default usmi " << key << endl;
    m.invalidate_smiles();        // we don't want (necessarily) to write a unique smiles
  }

  assert (key.length());

  return 1;
}

int
Storage_Conditions::transfer_to_mol2graph(Mol2Graph & mol2graph) const
{
  if (! _tautomer)
    return 1;

  mol2graph.set_remove_chiral_centres(_remove_chirality);
  mol2graph.set_exclude_triple_bonds_from_graph_reduction(_exclude_triple_bonds_from_graph_reduction);
  mol2graph.set_preserve_cc_double_bonds_saturated(_exclude_cc_double_bonds_saturated_from_graph_reduction);
  mol2graph.set_preserve_cc_double_bonds_no_heteroatoms(_exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction);
  mol2graph.set_aromatic_distinguishing_formula(_use_aromatic_distinguishing_mf_in_tautomer);
  mol2graph.set_revert_all_directional_bonds_to_non_directional(_remove_cis_trans_bonds);
  mol2graph.set_append_molecular_formula(1);

  return 1;
}

AllVariants::AllVariants() {
  _active = 0;
  _fragment_modified = "F/";
  _chiral_modified = "@";
}

static void
DisplayDashVOptions(std::ostream& output) {
  output << "The -V option controls building of 3 structural variants\n";
  output << "The starting molecule\n";
  output << "The starting molecule stripped to the largest fragment\n";
  output << "The starting molecule stripped to the largest fragment and chirality dropped\n";
  output << "The names associated with these variants will have a prefix indicating what\n";
  output << " -V def            Take defaults for everything\n";
  output << " -V lfrag=<s>      Prepend <s> to identifiers that have been stripped to the largest fragment\n";
  output << " -V chiral=<s>     Prepend <s> to identifiers that have lost chirality\n";
  output << " -V none           Do not prepend anything to identifiers\n";

  ::exit(0);
}


int
AllVariants::Initialise(Command_Line& cl, char flag) {
  int verbose = cl.option_present('v');

  const_IWSubstring v;
  for (int i = 0; cl.value('V', v, i); ++i) {
    if (v == "def") {
      
    } else if (v.starts_with("chiral=")) {
      v.remove_leading_chars(7);
      _chiral_modified = v;
    } else if (v.starts_with("lfrag=")) {
      v.remove_leading_chars(6);
      _fragment_modified = v;
    } else if (v == "none") {
      _chiral_modified.resize(0);
      _fragment_modified.resize(0);
    } else if (v == "help") {
      DisplayDashVOptions(cerr);
    } else {
      cerr << "Unrecognised -V qualifier '" << v << "'\n";
      DisplayDashVOptions(cerr);
    }
  }

  _active = 1;

  if (verbose) {
    cerr << "AllVariants activated\n";
  }

  return 1;
}

static void
PrependToName(const IWString& prefix, Molecule& m) {
  IWString new_name;
  new_name.resize(prefix.size() + m.name().size());
  new_name << prefix << m.name();
  m.set_name(new_name);
}

void
AllVariants::PrependFragmentModified(Molecule& m) const {
  PrependToName(_fragment_modified, m);
}

void
AllVariants::PrependChiralModified(Molecule& m) const {
  PrependToName(_chiral_modified, m);
}
