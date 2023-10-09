/*
  Load smiles into a BerkeleyDb database
*/
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "db_cxx.h"
#include "storage_conditions.h"

#include "sys/stat.h"
#include "sys/types.h"

using std::cerr;

static const char *prog_name = nullptr;

static Chemical_Standardisation chemical_standardisation;

static Charge_Assigner charge_assigner;

static Element_Transformations element_transformations;

static int verbose = 0;

static int store_multi_fragment_molecules_as_well_as_largest_fragment = 0;

static int multi_fragment_molecules_stored = 0;

// static int store_chiral_smiles = 1;

static int store_chiral_and_non_chiral_forms = 0;

static int molecules_read = 0;

static int molecules_stored = 0;

static IWString name_separator(':');

static int store_duplicates = 1;

static int duplicates_not_stored = 0;

static Molecule_Output_Object stream_for_duplicates;
static Molecule_Output_Object stream_for_unique_molecules;

static int data_appended = 0;

static int identical_entries_not_stored = 0;

static Report_Progress report_progress;

static int nostructs_ignored = 0;

static Db database(NULL, DB_CXX_NO_EXCEPTIONS);

static Db formula_db(NULL, DB_CXX_NO_EXCEPTIONS);

static int formula_database_active = 0;

static int formulae_stored = 0;

/*
  If we are doing delta-type processing, we need to keep track
  of all the keys that have changed during a load
*/

static IW_STL_Hash_Set changed_keys;

static IWString_and_File_Descriptor stream_for_changed_keys;

/*
  If we are writing to a hash instead of a database
*/

static int store_in_hash = 0;

static IW_STL_Hash_Map_String db_hash;

// Apr 2022. Allow for only storing some of the input stream of molecules.
// Drop any that have too many atoms.
static int max_atoms = std::numeric_limits<int>::max();
static int discarded_for_max_atoms = 0;

// Jul 2023.
// Allow building a database that contains all structural variants
// other than the graph.
AllVariants all_variants;

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
  cerr << "Builds a BerkeleyDb database with the unique smiles as the key\n";
  cerr << "Usage: " << prog_name << " <options> <input_file>...\n";
  cerr << "  -d <dbname>    specify BerkeleyDb database\n";
  cerr << "  -S <fname>     write to a hash and then write hash\n";
  cerr << "  -O <dbname>    optional database for molecular formulae\n";
  cerr << "  -H ...         specify storage conditions, enter '-H help' for info\n";
//cerr << "  -a             store graph form\n";
//cerr << "  -w             use aromatic distinguishing molecular formula for tautomer\n";
//cerr << "  -q             exclude triple bonds from graph reduction\n";
//cerr << "  -h s           exclude carbon-carbon double bonds from graph reduction, fully saturated nbrs\n";
//cerr << "  -h c           exclude carbon-carbon double bonds from graph reduction, all carbon nbrs\n";
//cerr << "  -c             exclude chirality\n";
//cerr << "  -z             exclude cis-trans bonding information\n";
  cerr << "  -b             store chiral and non-chiral smiles\n";
  cerr << "  -l             strip to largest fragment\n";
//cerr << "  -I             remove isotopes before storing\n";
//cerr << "  -u             ignore no-structs\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -n <separator> separator for when storing duplicate entries\n";
  cerr << "  -p             don't store duplicates\n";
  cerr << "  -D <file>      write duplicates to <file>\n";
  cerr << "  -U <file>      write non-duplicates to <file>\n";
  cerr << "  -o <type>      specify output type(s) for -D and -U options\n";
  cerr << "  -r <number>    report progress every <number> molecules\n";
  cerr << "  -C <fname>     write changed keys to <fname>\n";
  cerr << "  -x <natoms>    discard molecules containing more than <natoms> atoms\n";
  cerr << "  -V ...         store all structural variants, enter '-V help' for info\n";
  cerr << "  -y             suppress messages about bad valences\n";
//(void) display_standard_charge_assigner_options (cerr, 'N');
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -T ...         standard element transformation options\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
process_changed_key(const IWString &key) {
  if (changed_keys.contains(key)) {
    return 1;
  }

  changed_keys.insert(key);

  return 1;
}

static int
write_changed_key(Db &database, const IWString &key,
                  IWString_and_File_Descriptor &stream_for_changed_keys) {
  Dbt dkey((char *)key.rawchars(), key.length());
  Dbt fromdb;

  if (0 != database.get(NULL, &dkey, &fromdb, 0)) {
    cerr << "Cannot retrieve key '" << key << "' from database\n";
    return 0;
  }

  const_IWSubstring tmp(reinterpret_cast<const char *>(fromdb.get_data()),
                        fromdb.get_size());

  stream_for_changed_keys << key << ' ' << tmp << '\n';

  stream_for_changed_keys.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
write_changed_keys(Db &database, const IW_STL_Hash_Set &changed_keys,
                   IWString_and_File_Descriptor &stream_for_changed_keys) {
  for (IW_STL_Hash_Set::const_iterator i = changed_keys.begin(); i != changed_keys.end();
       ++i) {
    if (!write_changed_key(database, (*i), stream_for_changed_keys)) {
      return 0;
    }
  }

  if (verbose) {
    cerr << "Wrote " << changed_keys.size() << " changed keys\n";
  }

  return 1;
}

static int
do_store(Db &database, Dbt &dbkey, Dbt &dbdata) {
  const int rc = database.put(NULL, &dbkey, &dbdata, 0);

  if (0 != rc) {
    const_IWSubstring key(reinterpret_cast<const char *>(dbkey.get_data()),
                          dbkey.get_size());
    const_IWSubstring zdata(reinterpret_cast<const char *>(dbdata.get_data()),
                            dbdata.get_size());

    cerr << "Error, cannot store '" << key << "' value '" << zdata << "'\n";
    database.err(rc, "");
    return 0;
  }

  molecules_stored++;

  return 1;
}

static int
do_store(Db &database, Dbt &dbkey, const IWString &dbdata) {
  Dbt tostore;

  tostore.set_data((char *)dbdata.rawchars());
  tostore.set_size(dbdata.length());

  return do_store(database, dbkey, tostore);
}

static void
preprocess(Molecule &m, const Storage_Conditions &sc) {
  if (element_transformations.number_elements()) {
    (void)element_transformations.process(&m);
  }

  (void)chemical_standardisation.process(m);

  if (charge_assigner.active()) {  // do this after chemical standardisation
    (void)charge_assigner.process(m);
  }

  if (sc.reduce_to_largest_fragment()) {
    m.reduce_to_largest_fragment_carefully();  // note non-standard use - OK
  }

  if (sc.convert_isotopes()) {
    (void)m.transform_to_non_isotopic_form();
  }

  if (sc.remove_cis_trans_bonds()) {
    m.revert_all_directional_bonds_to_non_directional();
  }

  m.unset_unnecessary_implicit_hydrogens_known_values();

  return;
}

static int
store_formula(Molecule &m) {
  IWString f;
  m.formula_distinguishing_aromatic(f);

  Dbt zkey;
  zkey.set_data(const_cast<char *>(f.rawchars()));
  zkey.set_size(f.length());

  Dbt zdata;

  if (0 ==
      formula_db.get(NULL, &zkey, &zdata, 0))  // formula already there, don't do anything
  {
    delete reinterpret_cast<char *>(zdata.get_data());
    return 1;
  }

  // We just store an arbitrary string in the formula db

  zdata.set_data(const_cast<char *>("1"));
  zdata.set_size(1);

  int rc = formula_db.put(NULL, &zkey, &zdata, 0);
  if (0 == rc) {
    formulae_stored++;
    return 1;
  }

  cerr << "Cannot store into formula db '" << f << "'\n";

  return 0;
}

static int
do_store_in_hash(Molecule &m, const IWString &key, const IWString &mname,
                 IW_STL_Hash_Map_String &db_hash) {
  IW_STL_Hash_Map_String::iterator f = db_hash.find(key);

  if (f == db_hash.end()) {
    if (verbose > 1) {
      cerr << "New entry in database '" << mname << "'\n";
    }

    if (stream_for_unique_molecules.active()) {
      stream_for_unique_molecules.write(m);
    }

    molecules_stored++;

    db_hash[key] = mname;

    return db_hash.size();
  }

  // An entry with this key already in the database.

  if (stream_for_duplicates.active()) {
    stream_for_duplicates.write(m);
  }

  if (!store_duplicates) {
    duplicates_not_stored++;

    if (verbose > 1) {
      cerr << "Already in database, not stored\n";
    }
    return 1;
  }

  // If this is the same as what's already in the database, don't store it again

#ifdef DEBUG_NAME_COMPARISON
  cerr << "already_in_database.size is " << already_in_database.get_size() << '\n';
  const_IWSubstring junk((const char *)(already_in_database.get_data()),
                         already_in_database.get_size());
  cerr << "Stored '" << junk << "', trying to store '" << mname << "'\n";
#endif

  if ((*f).second == mname) {
    if (verbose > 1) {
      cerr << "Identical to current contents, not stored\n";
    }

    identical_entries_not_stored++;

    return 1;
  }

  data_appended++;

  IWString tmp = (*f).second;

  if (verbose > 1) {
    cerr << "Already in database, appending '" << mname << "' to '" << tmp << "'\n";
  }

  tmp << name_separator << mname;

  if (stream_for_changed_keys.active()) {
    process_changed_key(key);
  }

  molecules_stored++;

  db_hash[key] = tmp;

  return db_hash.size();
}

static int
buildsmidb(Molecule &m, const Storage_Conditions &sc, int include_chirality,
           const Mol2Graph &mol2graph, char prepend = ' ') {
  IWString key;
  IWString molecular_formula;

  if (formula_database_active) {
    // even if we don't store our molecule because it is already present, we
    // would want the formula database updated if the entry is absent
    store_formula( m);
  }

  // If we are storing the graph, we need to be careful to not mess up any output

  sc.form_key(m, mol2graph, include_chirality, key);

  assert(key.length() > 0);

  IWString mname = m.name();
  if (' ' != prepend) {
    mname.insert_at_beginning(prepend);
  }

  if (store_in_hash) {
    return do_store_in_hash(m, key, mname, db_hash);
  }

  Dbt dbkey;

  dbkey.set_data(const_cast<char *>(key.rawchars()));  // loss of const OK
  dbkey.set_size(key.length());

  Dbt already_in_database;

  if (0 != database.get(NULL, &dbkey, &already_in_database, 0)) {
    if (verbose > 1) {
      cerr << "New entry in database '" << mname << "'\n";
    }

    if (stream_for_unique_molecules.active()) {
      stream_for_unique_molecules.write(m);
    }

    return do_store(database, dbkey, mname);
  }

  // An entry with this key already in the database.

  IWString tmp;
  // tmp.set_and_assume_ownership(reinterpret_cast<char
  // *>(already_in_database.get_data()), already_in_database.get_size());
  tmp.set(reinterpret_cast<char *>(already_in_database.get_data()),
          static_cast<int>(already_in_database.get_size()));

  if (stream_for_duplicates.active()) {
    stream_for_duplicates.write(m);
  }

  if (!store_duplicates) {
    duplicates_not_stored++;

    if (verbose > 1) {
      cerr << "Already in database, not stored\n";
    }
    return 1;
  }

  // If this is the same as what's already in the database, don't store it again

#ifdef DEBUG_NAME_COMPARISON
  cerr << "already_in_database.size is " << already_in_database.get_size() << '\n';
  const_IWSubstring junk((const char *)(already_in_database.get_data()),
                         already_in_database.get_size());
  cerr << "Stored '" << junk << "', trying to store '" << mname << "'\n";
#endif

  if (already_in_database.get_size() != static_cast<unsigned int>(mname.length())) {
    ;
  } else if (0 ==
             mname.strncmp(reinterpret_cast<const char *>(already_in_database.get_data()),
                           already_in_database.get_size())) {
    if (verbose > 1) {
      cerr << "Identical to current contents, not stored\n";
    }

    identical_entries_not_stored++;

    return 1;
  }

  data_appended++;

  if (verbose > 1) {
    cerr << "Already in database, appending '" << mname << "' to '" << tmp << "'\n";
  }

  tmp.resize(tmp.length() + name_separator.length() + mname.length() + 3);

  tmp << name_separator << mname;

  if (stream_for_changed_keys.active()) {
    process_changed_key(key);
  }

  return do_store(database, dbkey, tmp);
}

void
PrependToName(const IWString& prefix, Molecule& m) {
  IWString new_name;
  new_name.resize(prefix.size() + m.name().size());
  new_name << prefix << m.name();
  m.set_name(new_name);
}

// Store the unique smiles of `m` in `database`. Value is
// `m.name`.
static int
do_store(Molecule& m, Db& database) {
  const IWString& usmi = m.unique_smiles();
  Dbt dbkey((void*) usmi.data(), usmi.length());

  Dbt already_in_database;

  if (0 != database.get(NULL, &dbkey, &already_in_database, 0)) {
    if (verbose > 1) {
      cerr << "New entry in database '" << m.name() << "'\n";
    }

    if (stream_for_unique_molecules.active()) {
      stream_for_unique_molecules.write(m);
    }

    return do_store(database, dbkey, m.name());
  }

  IWString new_name;
  new_name.strncat(reinterpret_cast<const char*>(already_in_database.get_data()),
                   already_in_database.get_size());
  new_name << name_separator << m.name();

  return do_store(database, dbkey, new_name);
}

int
BuildAllStructuralVariants(Molecule& m,
                        const AllVariants& all_variants,
                        Db& database) {
  if (! do_store(m, database)) {
    return 0;
  }

  if (m.number_fragments() > 1) {
    m.reduce_to_largest_fragment_carefully();
    all_variants.PrependFragmentModified(m);
    if (! do_store(m, database)) {
      return 0;
    }
  }

  if (m.chiral_centres() > 0) {
    m.remove_all_chiral_centres();
    all_variants.PrependChiralModified(m);
    if (! do_store(m, database)) {
      return 0;
    }
  }

  return 1;
}

static int
buildsmidb(Molecule &m, const Storage_Conditions &sc, const Mol2Graph &mol2graph) {
  int rc;
  if (all_variants.active()) {
    rc = BuildAllStructuralVariants(m, all_variants, database);
  } else if (sc.tautomer()) {
    rc = buildsmidb(m, sc, 0, mol2graph);
  } else if (store_chiral_and_non_chiral_forms) {
    rc = buildsmidb(m, sc, 1, mol2graph);
    if (rc && m.chiral_centres()) {
      rc = buildsmidb(m, sc, 0, mol2graph, '@');
    }
  } else if (!sc.remove_chirality()) {
    rc = buildsmidb(m, sc, 1, mol2graph);
  } else {
    rc = buildsmidb(m, sc, 0, mol2graph);
  }

  if (0 == rc) {
    return 0;
  }

  if (store_multi_fragment_molecules_as_well_as_largest_fragment &&
      m.number_fragments() > 1) {
    m.reduce_to_largest_fragment_carefully();
    IWString tmp;
    tmp << "FRAG:" << m.name();
    m.set_name(tmp);
    multi_fragment_molecules_stored++;
    return buildsmidb(m, sc, mol2graph);
  }

  return rc;
}

static int
buildsmidb(data_source_and_type<Molecule> &input, const Storage_Conditions &sc,
           const Mol2Graph &mol2graph) {
  Molecule *m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (report_progress()) {
      cerr << "read " << molecules_read << " molecules, stored " << molecules_stored << '\n';
    }

    preprocess(*m, sc);

    if (0 == m->natoms()) {
      nostructs_ignored++;
      continue;
    }

    if (m->natoms() > max_atoms) {
      ++discarded_for_max_atoms;
      continue;
    }

    if (!buildsmidb(*m, sc, mol2graph)) {
      return 0;
    }
  }

  return 1;
}

static int
buildsmidb(const char *fname, FileType input_type, const Storage_Conditions &sc,
           const Mol2Graph &mol2graph) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  int rc = buildsmidb(input, sc, mol2graph);

  if (verbose) {
    cerr << molecules_read << " molecules read, " << molecules_stored
         << " molecules stored\n";
  }

  return rc;
}

#ifdef NOT_USED
static void
foobar(Db &database) {
  Dbt dkey((void *)"foo", 3);
  Dbt zdata;

  int jjj = database.get(NULL, &dkey, &zdata, 0);
  cerr << "Return code is " << jjj << '\n';
  if (0 != jjj) {
    database.err(jjj, "");
  }

  return;
}
#endif

static int
do_write_hash(const IW_STL_Hash_Map_String &db_hash,
              IWString_and_File_Descriptor &output) {
  for (IW_STL_Hash_Map_String::const_iterator i = db_hash.begin(); i != db_hash.end();
       ++i) {
    output << (*i).first << ' ' << (*i).second << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

static int
do_write_hash(const IW_STL_Hash_Map_String &db_hash, const char *fname) {
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "Cannot open database output hash '" << fname << "'\n";
    return 0;
  }

  return do_write_hash(db_hash, output);
}

static int
buildsmidb(int argc, char **argv) {
  Command_Line cl(argc, argv, "vi:A:g:N:h:T:d:alE:n:D:U:po:Icr:bzuO:tLwC:S:qyH:x:V:");

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

  Storage_Conditions sc;
  sc.set_remove_cis_trans_bonds(1);  // until I get cis-trans unique smiles working
  Mol2Graph mol2graph;

  if (cl.option_present('l')) {
    sc.set_reduce_to_largest_fragment(1);
    if (verbose) {
      cerr << "Will strip to largest fragment before doing lookup\n";
    }
  }

  if (cl.option_present('I')) {
    sc.set_convert_isotopes(1);
    if (verbose) {
      cerr << "Isotopes stripped\n";
    }
  }

  if (cl.option_present('x')) {
    if (!cl.value('x', max_atoms) || max_atoms < 1) {
      cerr << "The max atoms optin (-x) must be a whole +ve number\n";
      return 1;
    }
    if (verbose) {
      cerr << "Will not store molecules having more than " << max_atoms << " atoms\n";
    }
  }

  if (cl.option_present('H')) {
    if (!sc.initialise(cl, 'H', mol2graph, verbose)) {
      cerr << "Cannot initialise storage conditions (-H)\n";
      return 1;
    }
  }

  if (cl.option_present('S')) {
    store_in_hash = 1;
    if (verbose) {
      cerr << "Will store in a hash then write\n";
    }
  } else if (!cl.option_present('d')) {
    cerr << "Must specify database via -d option\n";
    usage(8);
  }

  if (cl.option_present('d')) {
    const char *dbname = cl.option_value('d');

    int flags;
    DBTYPE dbtype;
    int mode;

    if (dash_s(dbname)) {
      dbtype = DB_UNKNOWN;
      flags = 0;
      mode = 0;
    } else {
      dbtype = DB_BTREE;
      flags = DB_CREATE;
      mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;
    }

    int rc = database.open(NULL, dbname, NULL, dbtype, flags, mode);

    if (0 != rc) {
      cerr << "Cannot open database '" << dbname << "'\n";
      database.err(rc, "");
      return 2;
    }

    if (verbose) {
      cerr << "Smiles will be written to database '" << dbname << "'\n";
    }
  }

  // The formula database idea really does not work well.
  if (cl.option_present('O')) {
    const char *dbname = cl.option_value('O');

    int flags = DB_CREATE;

    int mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

    DBTYPE dbtype = DB_BTREE;

    int rc = formula_db.open(NULL, dbname, NULL, dbtype, flags, mode);

    if (0 != rc) {
      cerr << "Cannot open formula database '" << dbname << "'\n";
      formula_db.err(rc, "");
      return 2;
    }

    if (verbose) {
      cerr << "Molecular formulae will be written to database '" << dbname << "'\n";
    }

    formula_database_active = 1;
  }

  if (cl.option_present('p')) {
    store_duplicates = 0;
    if (verbose) {
      cerr << "Will not store any duplicate entry\n";
    }
  }

  if (cl.option_present('y')) {
    set_display_abnormal_valence_messages(0);
    set_display_strange_chemistry_messages(0);

    if (verbose) {
      cerr << "Will suppress invalid valence messages\n";
    }
  }

  if (cl.option_present('V')) {
    if (! all_variants.Initialise(cl, 'V')) {
      cerr << "Cannot initialise all structure variants (-V)\n";
      return 1;
    }

    sc.set_all_variants(1);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (!cl.option_present('i') && 1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (0 == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (cl.option_present('n')) {
    (void)cl.value('n', name_separator);
    if (verbose) {
      cerr << "Will store names separated by '" << name_separator << "'\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('D')) {
    const_IWSubstring fname;
    cl.value('D', fname);

    stream_for_duplicates.add_output_type(FILE_TYPE_SMI);

    if (!stream_for_duplicates.new_stem(fname, 1)) {
      cerr << "Cannot open duplicates file '" << fname << "'\n";
      return 87;
    }

    if (verbose) {
      cerr << "Duplicates written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('U')) {
    const_IWSubstring fname;
    cl.value('U', fname);

    stream_for_unique_molecules.add_output_type(FILE_TYPE_SMI);

    if (!stream_for_unique_molecules.new_stem(fname, 1)) {
      cerr << "Cannot open unique molecules file '" << fname << "'\n";
      return 87;
    }

    if (verbose) {
      cerr << "Unique molecules written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('C')) {
    const char *c = cl.option_value('C');

    if (!stream_for_changed_keys.open(c)) {
      cerr << "Cannot open stream for changed keys '" << c << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Will write changed keys to '" << c << "'\n";
    }
  }

  if (cl.option_present('b') && (cl.option_present('c') || cl.option_present('a'))) {
    cerr << "The -b option is inconsistent with the -c and -a options\n";
    usage(6);
  }

  if (cl.option_present('z')) {
    sc.set_remove_cis_trans_bonds(1);

    if (verbose) {
      cerr << "Cis trans bonding information suppressed\n";
    }
  }

  sc.set_remove_cis_trans_bonds(1);  // until I get those unique smiles working

  if (cl.option_present('a')) {
    sc.set_tautomer(1);
    if (verbose) {
      cerr << "Will store molecular graphs\n";
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
  } else if (cl.option_present('c')) {
    sc.set_remove_chirality(1);

    if (verbose) {
      cerr << "Chirality information excluded from smiles\n";
    }
  } else if (cl.option_present('b')) {
    sc.set_remove_chirality(0);
    store_chiral_and_non_chiral_forms = 1;

    if (verbose) {
      cerr << "Both chiral and non-chiral unique smiles stored\n";
    }
  } else if (cl.option_present('H')) {  // specified
    ;
  } else {
    sc.set_remove_chirality(0);

    if (verbose) {
      cerr << "All structural information stored\n";
    }
  }

  if (store_in_hash) {
    ;
  } else if (!sc.store_or_check_currently_stored(database, verbose)) {
    cerr
        << "Warning, storage conditions being used incompatible with stored conditions\n";
  }

  if (cl.option_present('L')) {
    store_multi_fragment_molecules_as_well_as_largest_fragment = 1;

    sc.set_reduce_to_largest_fragment(0);

    if (verbose) {
      cerr << "Will store multi fragment molecules as well as largest fragment\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(3);
    }
  }

  if (0 == verbose) {
    set_display_strange_chemistry_messages(0);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!buildsmidb(cl[i], input_type, sc, mol2graph)) {
      cerr << "Error processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules, stored " << molecules_stored
         << '\n';
    if (discarded_for_max_atoms > 0) {
      cerr << discarded_for_max_atoms << " molecules discarded for > " << max_atoms
           << " atoms\n";
    }
    if (name_separator.length()) {
      cerr << data_appended << " identifiers appended to existing database entries\n";
    }
    if (identical_entries_not_stored) {
      cerr << identical_entries_not_stored << " identical entries not stored\n";
    }
    if (!store_duplicates) {
      cerr << duplicates_not_stored << " duplicate molecules not stored\n";
    }

    if (formula_database_active) {
      cerr << formulae_stored << " molecular formulae stored\n";
    }

    if (store_multi_fragment_molecules_as_well_as_largest_fragment) {
      cerr << multi_fragment_molecules_stored << " multi fragment molecules stored\n";
    }
  }

  if (stream_for_changed_keys.active()) {
    write_changed_keys(database, changed_keys, stream_for_changed_keys);
  }

  if (store_in_hash) {
    if (verbose) {
      cerr << "hash contains " << db_hash.size() << " items\n";
    }

    do_write_hash(db_hash, cl.option_value('S'));
  } else {
    database.close(0);
  }

  if (formula_database_active) {
    formula_db.close(0);
  }

  return rc;
}

int
main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = buildsmidb(argc, argv);

  return rc;
}
