/*
  We often want to know whether of not a given molecule is in a database
*/

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "db_cxx.h"
#include "storage_conditions.h"

using std::cerr;

const char* prog_name = nullptr;

static Chemical_Standardisation chemical_standardisation;

static Charge_Assigner charge_assigner;

static Element_Transformations element_transformations;

static int verbose = 0;

static Molecule_Output_Object stream_for_in_database;

static Molecule_Output_Object stream_for_not_in_database;

// this is not necessary, ignoring errors can be done with -i ICTE
static bool SkipSmilesErrors = false;

//  Sometimes we just want molecules written with a marker as to
//  whether or not they are in the database

static Molecule_Output_Object stream_for_marked;

// Just check things, so no output is OK

static int ok_no_output = 0;

static IWString marker("IN DATABASE");

static int append_database_identifier = 0;

// static int strip_to_largest_fragment = 0;

// static int do_graph_lookup = 0;

// static int use_aromatic_distinguishing_mf_in_tautomer = 0;

static int lookup_chiral_and_non_chiral_smiles = 0;

static int lookup_multi_fragment_molecules_as_well_as_largest_fragment = 0;

// static int ignore_isotopes = 0;

// static int ignore_cis_trans_bonds = 0;

static int molecules_read = 0;

static int molecules_in_database = 0;

static int molecules_not_in_database = 0;

// The number of databases we open

static int ndb = 0;

//  Pointers to the opened databases

static Db** database;  //(NULL, DB_CXX_NO_EXCEPTIONS);

static Db formula_db(NULL, DB_CXX_NO_EXCEPTIONS);

static int formula_database_active = 0;

// The name of each database

static IWString* dbname = nullptr;

// How many items found in each database

static int* found_in_database = nullptr;

// It can be interesting to see how the molecules found in the
// database vary with natoms.

static extending_resizable_array<int> atoms_in_input;
static extending_resizable_array<int> atoms_in_matches;

//  By default, we stop once we have found a match

static int examine_all_databases = 0;

static char examine_all_databases_separator = ' ';

static int echo_db_key = 0;  // useful for debugging

static IWString_and_File_Descriptor stream_for_tabular_output;

// Jul 2023. Look up all variants until a match found.
static AllVariants all_variants;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Checks for molecules in a Berkeley database by looking up the unique smiles\n";
  cerr << "Usage: " << prog_name << " -d <dbname> <options> <input_file>...\n";
  cerr << "  -d <dbname>    specify Berkeley database\n";
  cerr << "  -j             discern lookup transformations from stored info (recommended)\n";
//cerr << "  -a             base lookup on graph form\n";
//cerr << "  -w             use aromatic distinguishing molecular formula for tautomer\n";
//cerr << "  -q             exclude triple bonds from graph reduction\n";
//cerr << "  -c             exclude chirality\n";
  cerr << "  -H ...         specify storage conditions, enter '-H help' for info\n";
  cerr << "  -V ...         store all structural variants, enter '-V help' for info\n";
//cerr << "  -z             exclude cis-trans bonding information\n";
  cerr << "  -b             look for both chiral and non-chiral smiles\n";
  cerr << "  -l             strip to largest fragment\n";
//cerr << "  -I             transform to non-isotopic form before doing lookup\n";
  cerr << "  -F <file>      stream for molecules found in the database\n";
  cerr << "  -p             append database identifier to name\n";
  cerr << "  -e <sep>       search all databases, insert <sep> string between database finds\n";
  cerr << "  -U <file>      stream for molecules NOT found in the database\n";
  cerr << "  -M <file>      stream for molecules with database marker\n";
  cerr << "  -D <file>      write tabular output of search results\n";
  cerr << "  -m <string>    specify database marker (default '" << marker << "')\n";
  cerr << "  -W <fname>     write per atom count lookup rates to <fname>\n";
  cerr << "  -n             ok to have no output - informational lookup only\n";
  cerr << "  -Y             echo the DB key (useful for debugging)\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  // cerr << "  -K             skip smiles errors\n";
//(void) display_standard_charge_assigner_options(cerr, 'N');
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -T ...          standard element transformation options\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  The database might have FRAG: strings prepended on stored identifiers.
  If we are looking up a largest fragment, then we want to keep that
  info. Otherwise we want to remove it
*/

static int
do_append_database_identifier(IWString& tmp, const Dbt& from_database,
                              int looking_up_largest_fragment, const char separator) {
  if (0 == tmp.length()) {
    cerr << "No name for input molecule, appending DB name\n";
  } else {
    tmp += separator;
  }

  if (0 == from_database.get_size()) {
    cerr << "Yipes, zero length field in database\n";
  } else if (looking_up_largest_fragment) {
    tmp.strncat(reinterpret_cast<const char*>(from_database.get_data()),
                from_database.get_size());
  } else {
    //  cerr << "Contents '";
    //  cerr.write(reinterpret_cast<const char *>(from_database.get_data()),
    //  from_database.get_size()); cerr << '\n';

    IWString s(reinterpret_cast<const char*>(from_database.get_data()),
               from_database.get_size());
    s.gsub("FRAG:", "");
    tmp << s;
  }

  return 1;
}

static int
do_append_database_identifier(Molecule& m, const Dbt& from_database,
                              int looking_up_largest_fragment, const char separator) {
  IWString tmp(m.name());

  int rc = do_append_database_identifier(tmp, from_database, looking_up_largest_fragment,
                                         separator);

  m.set_name(tmp);

  return rc;
}

static int
check_database(const IWString& key, int& which_database, Dbt& from_database) {
  Dbt dkey;

  dkey.set_data((void*)(key.rawchars()));  // loss of const OK
  dkey.set_size(key.nchars());

  for (int i = 0; i < ndb; i++) {
    if (0 != database[i]->get(NULL, &dkey, &from_database, 0)) {
      continue;
    }

    which_database = i;
    found_in_database[i]++;
    return 1;
  }

  if (verbose > 2) {
    cerr << "No match for '" << key << "'\n";
  }

  return 0;
}

static int
do_check_all_databases(const IWString& skey, 
                       const int looking_up_largest_fragment, IWString& db_identifiers) {
  Dbt dkey;

  dkey.set_data((void*)(skey.rawchars()));  // loss of const OK
  dkey.set_size(skey.nchars());

  int rc = 0;

  for (int i = 0; i < ndb; ++i) {
    Dbt from_database;

    if (0 != database[i]->get(NULL, &dkey, &from_database, 0)) {
      continue;
    }

    found_in_database[i]++;
    rc++;

    if (append_database_identifier) {
      do_append_database_identifier(db_identifiers, from_database,
                                    looking_up_largest_fragment,
                                    examine_all_databases_separator);
    }
  }

  if (rc) {
    if (append_database_identifier) {
    }
  }

  return rc;
}

/*
  Given the Molecule, form the key that we will use for looking up the molecule
  We don't perturb our input, so if we need to lookup a different form we must make a copy
*/

static int
form_key(Molecule& m, const Storage_Conditions& sc, const Mol2Graph& mol2graph,
         int include_chirality, IWString& key) {
  sc.form_key(m, mol2graph, include_chirality, key);

  assert(key.length());

  if (echo_db_key) {
    cerr << "Key " << key << " " << m.name() << "\n";
  }

  return 1;
}

static void
preprocess(Molecule& m, const Storage_Conditions& sc) {
  if (element_transformations.number_elements()) {
    (void)element_transformations.process(m);
  }

  (void)chemical_standardisation.process(m);

  if (charge_assigner.active()) {  // do this after chemical standardisation
    (void)charge_assigner.process(m);
  }

  m.revert_all_directional_bonds_to_non_directional();  // until I get cis-trans bonding
                                                        // fixed

  m.unset_unnecessary_implicit_hydrogens_known_values();

  return;
}

static int
handle_molecule_not_in_database(Molecule& m) {
  molecules_not_in_database++;

  if (verbose > 1) {
    cerr << m.name() << " not in database\n";
  }

  if (stream_for_not_in_database.active()) {
    stream_for_not_in_database.write(m);
  }

  if (stream_for_marked.active()) {
    stream_for_marked.write(m);
  }

  return 1;
}

static int
handle_molecule_found(Molecule& m, const int which_database = -1) {
  molecules_in_database++;

  if (verbose) {
    atoms_in_matches[m.natoms()]++;
    if (verbose > 1 && which_database >= 0) {
      cerr << m.name() << " in database '" << dbname[which_database] << "'\n";
    }
  }

  if (stream_for_in_database.active()) {
    stream_for_in_database.write(m);
  }

  if (stream_for_marked.active()) {
    IWString tmp = m.name();
    tmp << ' ' << marker;
    m.set_name(tmp);
    stream_for_marked.write(m);
  }

  return 1;
}

static int
in_database_examine_all_databases(Molecule& m,
                        const IWString& key,
                        int looking_up_largest_fragment) {
  IWString new_name;
  if (append_database_identifier) {
    new_name = m.name();
  }

  if (do_check_all_databases(key, looking_up_largest_fragment, new_name)) {
    if (append_database_identifier) {
      m.set_name(new_name);
    }

    return handle_molecule_found(m);
  } else {
    return 0;
  }
}

static int
in_database(Molecule& m, const Storage_Conditions& sc, int include_chirality,
            const Mol2Graph& mol2graph, int looking_up_largest_fragment = 0) {
  IWString key;

  form_key(m, sc, mol2graph, include_chirality, key);

  // cerr << "From " << m.name() << " key " << key << " include chirality? " <<
  // include_chirality << '\n';

  if (examine_all_databases) {
    return in_database_examine_all_databases(m, key, looking_up_largest_fragment);
  }

  // Stop once a database match is found

  int which_database;
  Dbt from_database;
  if (check_database(key, which_database, from_database)) {

    if (append_database_identifier) {
      do_append_database_identifier(m, from_database, looking_up_largest_fragment, ' ');
    }

    return handle_molecule_found(m, which_database);
  } else {
    return 0;
  }
}

// First lookup `m.unique_smiles`. If that fails, try stripping to
// the largest fragment, and then removing chirality in order to find
// a match.
// Returns 1 if a match is found, 0 otherwise.
int
LookupAllStructuralVariants(Molecule& m,
                const AllVariants& all_variants,
                const Storage_Conditions& sc,
                const Mol2Graph& mol2graph) {

  static constexpr int kIncludeChirality = 1;
  static constexpr int kLargestFragment = 0;

  IWString original_name(m.name());

  if (in_database(m, sc, kIncludeChirality, mol2graph, kLargestFragment)) {
    return 1;
  }

  if (m.number_fragments() > 1) {
    m.reduce_to_largest_fragment_carefully();
    all_variants.PrependFragmentModified(m);
    if (in_database(m, sc, kIncludeChirality, mol2graph, kLargestFragment)) {
      return 1;
    }
  }

  if (m.chiral_centres() > 0) {
    m.remove_all_chiral_centres();
    all_variants.PrependChiralModified(m);
    if (in_database(m, sc, kIncludeChirality, mol2graph, kLargestFragment)) {
      return 1;
    }
  }

  m.set_name(original_name);

  return 0;
}

static int
formula_in_database(Molecule& m, Db& formula_db) {
  IWString f;
  m.formula_distinguishing_aromatic(f);

  // cerr << "Checking formula '" << f << "'\n";

  Dbt zkey;
  zkey.set_data((void*)(f.rawchars()));
  zkey.set_size(f.length());

  Dbt zdata;

  if (0 != formula_db.get(NULL, &zkey, &zdata, 0)) {
    return 0;
  }

  delete reinterpret_cast<char*>(zdata.get_data());

  return 1;
}

static int
in_database(Molecule& m, const Storage_Conditions& sc, const Mol2Graph& mol2graph) {
  if (verbose) {
    atoms_in_input[m.natoms()]++;
  }

  if (!formula_database_active) {
    ;
  } else if (formula_in_database(m, formula_db)) {
    ;
  } else {
    return 1;
  }

  int found = 0;
  if (sc.tautomer()) {
    found = in_database(m, sc, 0, mol2graph);
  } else if (lookup_chiral_and_non_chiral_smiles) {
    found = in_database(m, sc, 1, mol2graph);
    if (!found && m.chiral_centres()) {
      m.invalidate_canonical_ordering_information();
      found = in_database(m, sc, 0, mol2graph);
    }
  } else if (all_variants.active()) {
    found = LookupAllStructuralVariants(m, all_variants, sc, mol2graph);
  } else {
    found = in_database(m, sc, !sc.remove_chirality(), mol2graph);
  }

  if (stream_for_tabular_output.is_open()) {
    stream_for_tabular_output << m.name() << ' ' << found << '\n';
    stream_for_tabular_output.write_if_buffer_holds_more_than(4096);
  }

  if (found) {
    return 1;
  }

  if (lookup_multi_fragment_molecules_as_well_as_largest_fragment &&
      m.number_fragments() > 1) {
    m.reduce_to_largest_fragment_carefully();
    // Last arg means leave FRAG: prefixes on retrieved identifiers.
    static constexpr int kLeaveFrag = 1;
    found = in_database( m, sc, 1, mol2graph, kLeaveFrag);
  }

  if (found) {
    return 1;
  }

  return handle_molecule_not_in_database(m);
}

static int
in_database(data_source_and_type<Molecule>& input, const Storage_Conditions& sc,
            const Mol2Graph& mol2graph) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m, sc);

    if (!in_database(*m, sc, mol2graph)) {
      return 0;
    }
  }

  return 1;
}

static int
in_database(const char* fname, FileType input_type, const Storage_Conditions& sc,
            const Mol2Graph& mol2graph) {
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }
  if (SkipSmilesErrors) {
    input.set_connection_table_errors_allowed(std::numeric_limits<int>::max());
  }

  return in_database(input, sc, mol2graph);
}

static int
handle_file_opening(Command_Line& cl, char flag, FileType output_type,
                    Molecule_Output_Object& zfile) {
  if (output_type) {
    zfile.add_output_type(output_type);
  } else if (!zfile.determine_output_types(cl)) {
    cerr << "Cannot determine output types\n";
    return 0;
  }

  IWString fname;
  cl.value(flag, fname);

  if (!zfile.new_stem(fname, 1)) {
    cerr << "Cannot open -" << flag << " file '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Molecules for -" << flag << " written to '" << fname << "'\n";
  }

  return 1;
}

template <typename T>
int
write_per_atom_lookup_rates(const extending_resizable_array<int>& atoms_in_input,
                            const extending_resizable_array<int>& atoms_in_matches,
                            T& os) {
  for (int i = 1; i < atoms_in_input.number_elements(); i++) {
    if (0 == atoms_in_input[i]) {
      continue;
    }

    int matches;
    if (i >= atoms_in_matches.number_elements()) {
      matches = 0;
    } else {
      matches = atoms_in_matches[i];
    }

    float ratio = static_cast<float>(matches) / static_cast<float>(atoms_in_input[i]);

    os << "Looked up " << atoms_in_input[i] << " molecules w/ " << i << " atoms, found "
       << matches << " ratio " << ratio << '\n';
  }

  return 1;
}

static int
write_per_atom_lookup_rates(const extending_resizable_array<int>& atoms_in_input,
                            const extending_resizable_array<int>& atoms_in_matches,
                            const char* fname) {
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "Cannot open per atom count file name '" << fname << "'\n";
    return 0;
  }

  return write_per_atom_lookup_rates(atoms_in_input, atoms_in_matches, output);
}

static int
in_database(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:o:A:g:N:h:T:d:alE:F:U:pM:m:wIncbzO:W:Ljqye:YH:D:KV:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      usage(6);
    }
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  // Just too many opportunities for error if we don't enforce a uniform
  // aromaticity definition

  set_global_aromaticity_type(Daylight);

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot determine charge assigner from command line\n";
      usage(77);
    }
  }

  if (cl.option_present('T')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T')) {
      usage(8);
    }
  }

  if (cl.option_present('b') && (cl.option_present('c') || cl.option_present('a'))) {
    cerr << "The -b option is inconsistent with the -c and -a options\n";
    usage(6);
  }

  Mol2Graph mol2graph;

  Storage_Conditions sc;

  if (cl.option_present('H')) {
    if (!sc.initialise(cl, 'H', mol2graph, verbose)) {
      cerr << "Cannot initialise storage conditions (-H)\n";
      return 1;
    }
  }

  if (cl.option_present('z')) {
    sc.set_remove_cis_trans_bonds(1);

    if (verbose) {
      cerr << "Cis trans bonding information suppressed\n";
    }
  }

  sc.set_remove_cis_trans_bonds(1);  // until I get cis-trans unique smiles working
  bool storage_conditions_on_command_line = false;

  if (cl.option_present('K')) {
    SkipSmilesErrors = true;
  }

  if (cl.option_present('a')) {
    sc.set_tautomer(1);
    if (verbose) {
      cerr << "Lookup will be based on a graph lookup\n";
    }

    if (!cl.option_present('g')) {
      cerr << "WARNING, graph lookup without chemical standardisation, very dangerous!\n";
    }

    if (cl.option_present('w')) {
      sc.set_use_aromatic_distinguishing_mf_in_tautomer(1);

      if (verbose) {
        cerr << "Will use formula_distinguishing_aromatic() for tautomer mf\n";
      }
    }

    if (cl.option_present('q')) {
      sc.set_exclude_triple_bonds_from_graph_reduction(1);

      mol2graph.set_exclude_triple_bonds_from_graph_reduction(1);

      set_exclude_triple_bonds_from_graph_reduction(1);

      if (verbose) {
        cerr << "Will not convert triple bonds in tautomer computations\n";
      }
    }

    if (cl.option_present('h')) {
      const_IWSubstring h = cl.string_value('h');

      if ('s' == h) {
        sc.set_exclude_cc_double_bonds_saturated_from_graph_reduction(1);
        mol2graph.set_preserve_cc_double_bonds_saturated(1);

        if (verbose) {
          cerr << "Will not convert C=C bonds adjacent to fully saturated Carbons in "
                  "tautomer computations\n";
        }
      } else if ('c' == h) {
        sc.set_exclude_cc_double_bonds_no_heteroatoms_from_graph_reduction(1);
        mol2graph.set_preserve_cc_double_bonds_no_heteroatoms(1);

        if (verbose) {
          cerr << "Will not convert C=C bonds adjacent to all Carbon atoms in tautomer "
                  "computations\n";
        }
      } else {
        cerr << "Unrecognised -h qualifier '" << h << "'\n";
        usage(1);
      }
    }
    storage_conditions_on_command_line = true;
  } else if (cl.option_present('c')) {
    sc.set_remove_chirality(1);

    if (verbose) {
      cerr << "Chirality information excluded from smiles\n";
    }
    storage_conditions_on_command_line = true;
  } else if (cl.option_present('b')) {
    lookup_chiral_and_non_chiral_smiles = 1;

    if (verbose) {
      cerr << "Will look up both chiral and non-chiral smiles\n";
    }
    storage_conditions_on_command_line = true;
  } else if (cl.option_present('V')) {
    if (! all_variants.Initialise(cl, 'V')) {
      cerr << "Cannot initialise all structure variants (-V)\n";
      return 1;
    }
    storage_conditions_on_command_line = false;
    sc.set_all_variants(1);
  } else {
    if (verbose) {
      cerr << "Chirality information included in smiles\n";
    }
  }

  if (cl.option_present('l')) {
    sc.set_reduce_to_largest_fragment(1);
    if (verbose) {
      cerr << "Will strip to largest fragment before doing lookup\n";
    }
    storage_conditions_on_command_line = true;
  }

  if (cl.option_present('L')) {
    lookup_multi_fragment_molecules_as_well_as_largest_fragment = 1;

    if (verbose) {
      cerr << "Will lookup multi fragment molecules as well as largest fragment\n";
    }

    sc.set_reduce_to_largest_fragment(0);
    storage_conditions_on_command_line = true;
  }

  if (cl.option_present('I')) {
    sc.set_convert_isotopes(1);
    if (verbose) {
      cerr << "Isotopes ignored for lookups\n";
    }
    storage_conditions_on_command_line = true;
  }

  ndb = cl.option_count('d');

  if (0 == ndb) {
    cerr << "Must specify one or more databases via -d option\n";
    usage(8);
  }

  int storage_conditions_must_match_database = 0;

  if (cl.option_present('j')) {
    storage_conditions_must_match_database = 1;

    if (verbose) {
      cerr << "Will check for storage condition incompatibilities\n";
    }
  }

  database = new Db*[ndb];
  std::unique_ptr<Db*[]> free_database(database);
  dbname = new IWString[ndb];
  std::unique_ptr<IWString[]> free_dbname(dbname);
  found_in_database = new_int(ndb);
  std::unique_ptr<int[]> free_found_in_database(found_in_database);

  if (cl.option_present('d')) {
    int i = 0;
    IWString d;
    while (cl.value('d', d, i)) {
      database[i] = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

      int rc = database[i]->open(NULL, d.null_terminated_chars(), NULL, DB_UNKNOWN,
                                 DB_RDONLY, 0);

      if (0 != rc) {
        cerr << "Cannot open database '" << d << "'\n";
        database[i]->err(rc, "");
        return i + 1;
      }

      if (!storage_conditions_must_match_database) {  // nothing to check
        ;
      } else if (!storage_conditions_on_command_line)  // get our storage conditions from
                                                       // the DB
      {
        sc.initialise("_STORE_INFO", *(database[i]));
        //      sc.debug_print(cerr);
        sc.transfer_to_mol2graph(mol2graph);
        sc.transfer_to_all_variants(all_variants);
      } else if (!sc.ensure_consistent_with_current_conditions(*(database[i]), verbose)) {
        cerr << "Database '" << d << " incompatible with lookup conditions\n";
        if (storage_conditions_must_match_database) {
          exit(3);
        }
      }

      dbname[i] = d;

      if (verbose) {
        cerr << "Database " << i << " is '" << d << "'\n";
      }

      i++;
    }
  }

  if (cl.option_present('O')) {
    const char* dbname = cl.option_value('O');

    int rc = formula_db.open(NULL, dbname, NULL, DB_UNKNOWN, DB_RDONLY, 0);
    if (0 != rc) {
      cerr << "Cannot open formula database '" << dbname << "'\n";
      formula_db.err(rc, "");
      return 4;
    }

    if (verbose) {
      cerr << "Will lookup formula in '" << dbname << "'\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (0 == input_type && 1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  }

  if (0 == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (!cl.option_present('F') && !cl.option_present('U') && !cl.option_present('M') &&
      !cl.option_present('n') && !cl.option_present('D')) {
    cerr << "Must specify -F, -U or -M for some output, or -n for no output\n";
    usage(18);
  }

  FileType output_type = FILE_TYPE_INVALID;
  if (!cl.option_present('o')) {
    output_type = FILE_TYPE_SMI;
    if (verbose) {
      cerr << "Output type defaults to SMI\n";
    }
  }

  if (cl.option_present('F')) {
    if (!handle_file_opening(cl, 'F', output_type, stream_for_in_database)) {
      return 13;
    }
  }

  if (cl.option_present('p')) {
    if (!cl.option_present('F') && !cl.option_present('M')) {
      cerr << "The -p option only makes sense with the -F or -M option\n";
      usage(13);
    }

    append_database_identifier = 1;
    if (verbose) {
      cerr << "Will append database identifiers\n";
    }
  }

  if (cl.option_present('U')) {
    if (!handle_file_opening(cl, 'U', output_type, stream_for_not_in_database)) {
      return 14;
    }
  }

  if (cl.option_present('m') && !cl.option_present('M')) {
    cerr << "The -m option (append string to database members) must be used with the -M "
            "option\n";
    usage(73);
  }

  if (cl.option_present('M')) {
    if (!handle_file_opening(cl, 'M', output_type, stream_for_marked)) {
      return 15;
    }
  }

  if (cl.option_present('m')) {
    cl.value('m', marker);
    if (verbose) {
      cerr << "Molecules in the database marked with '" << marker << "'\n";
    }
  }

  if (cl.option_present('D')) {
    const char* d = cl.option_value('D');

    if (!stream_for_tabular_output.open(d)) {
      cerr << "Cannot open stream for tabular output '" << d << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Tabular output written to '" << d << "'\n";
    }

    stream_for_tabular_output << "ID InDB\n";
  }

  if (cl.option_present('n')) {
    ok_no_output = 1;

    if (verbose) {
      cerr << "No structure output, just a lookup\n";
    }
  }

  if (cl.option_present('y')) {
    set_display_abnormal_valence_messages(0);
    set_display_strange_chemistry_messages(0);

    if (verbose) {
      cerr << "Will suppress invalid valence messages\n";
    }
  }

  if (cl.option_present('e')) {
    examine_all_databases = 1;

    IWString e = cl.string_value('e');
    if (! char_name_to_char(e)) {
      cerr << "Unrecognised or invalid -e qualifier '" << e << "'\n";
      return 1;
    }
    examine_all_databases_separator = e[0];
  }

  if (cl.option_present('Y')) {
    echo_db_key = 1;

    if (verbose) {
      cerr << "Will echo the database key as it is formed\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  assert(cl.ok());
  int rc = 0;
  for (const char* fname : cl) {
    if (!in_database(fname, input_type, sc, mol2graph)) {
      rc = 1;
      break;
    }
  }

  if (stream_for_in_database.active()) {
    stream_for_in_database.do_close();
  }

  if (verbose || ok_no_output) {
    cerr << "Examined " << molecules_read << " molecules\n";
    cerr << "Found " << molecules_in_database << " in database, "
         << molecules_not_in_database << " not in database\n";

    for (int i = 0; i < ndb; i++) {
      cerr << found_in_database[i] << " molecules found in database '" << dbname[i] << "'\n";
    }

    if (stream_for_not_in_database.active()) {
      cerr << "wrote " << stream_for_not_in_database.molecules_written()
           << " unique molecules to -U stream\n";
    }

    if (cl.option_present('W')) {
      const char* w = cl.option_value('W');
      write_per_atom_lookup_rates(atoms_in_input, atoms_in_matches, w);
    }
  }

  assert(cl.ok());

  if (nullptr != database) {
    for (int i = 0; i < ndb; i++) {
      database[i]->close(0);
      delete database[i];
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = in_database(argc, argv);

  return rc;
}
