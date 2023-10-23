/*
  Looks up molecules in a BerkeleyDB structure database.
  What is special about this variant is that it works in conjunction with
  inventory database(s).
  It is overly complex due to requirements that have emerged over the years.
*/

#include <iostream>
#include <limits>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "db_cxx.h"
#include "storage_conditions.h"

using std::cerr;

const char* prog_name = NULL;

static Chemical_Standardisation chemical_standardisation;

static int verbose = 0;

static Molecule_Output_Object stream_for_in_database;

static Molecule_Output_Object stream_for_not_in_database;

// this is not necessary, the right way to do this is '-i ICTE'.
static bool SkipSmilesErrors = false;

// This used to be an option, but now always set. The inventory databases
// have keys that do not have leading zero's.
static int remove_leading_zeros = 1;

// Historically there have been negative amounts in the inventory db's.
// Currently this is not being used.
static int negative_amounts_ignored = 0;

static int remove_directional_bonding = 1;  // until we start storing it

// write every molecule with an annotation of its fate. The -M option.

static Molecule_Output_Object stream_for_everything;

// Usually this will be desirable.
static int append_database_identifier = 0;

static int write_amounts_as_mg = 0;

static int append_all_matched_identifiers = 0;

// By default, LSN's that are found in the db come out with 'LSN-' prepended.
// That was not a good choice, but cannot be changed becuase it would be
// a breaking change.
static int prepend_lly_dash = 1;

static int write_all_inventory_database_matches = 0;

#ifdef DASH_A_OPTION_NOW_OBSOLETE
// Graph lookup options now done via the -H option
static int do_graph_lookup = 0;
static int use_aromatic_distinguishing_mf_in_tautomer = 0;
#endif  // DASH_A_OPTION_NOW_OBSOLETE

static int need_chirality = 1;

// If there are more than one inventory database active, we need
// a name for each one. These things should be put into an
// InventoryDatabase class.
// Then there should be a resizable_array_p of them.

static Db** inventory = NULL;
static IWString* dbname = nullptr;
static int number_inventory_databases = 0;
static int* inventory_is_liquid = nullptr;

static int lookup_chiral_and_non_chiral_smiles = 0;

static float solid_cutoff = static_cast<float>(0.0);
static int liquid_cutoff = 0;

static int molecules_read = 0;

static int molecules_in_database = 0;

static int molecules_not_in_database = 0;

static int molecules_in_database_but_below_cutoff = 0;

static int molecules_passing_constraints = 0;

// the structure database.
static Db database(NULL, DB_CXX_NO_EXCEPTIONS);

static int malformed_database_records = 0;

static int echo_key = 0;

/*
  It can be interesting to see how the molecules found in the
  database vary with natoms.
*/

static extending_resizable_array<int> atoms_in_input;
static extending_resizable_array<int> atoms_in_matches;

static IWString fname_for_per_natoms_lookup_report;

static int return_on_first_match = 1;

/*
  Since we are getting inventory data from multiple databases, and they
  may have different units, we need an object to hold fetch data
*/

template <typename T>
class Source_Data {
 private:
  T _amt;
  IWString _lly;
  IWString _dbname;

  //  private functions

  int _append_inventory_information(IWString&, T, const char*) const;

 public:
  Source_Data();

  int debug_print(std::ostream&) const;

  int new_amt(T, const IWString&, const IWString&);

  int append_inventory_information(IWString&, T multiplier, const char*) const;

  int append_inventory_information(IWString&, const char*) const;

  int
  has_data() const {
    return _lly.length();
  }

  void
  invalidate() {
    _lly.resize(0);
  }

  int satisfies_minimum_amounts(T) const;
};

/*
  We keep track of our best value from both the liquid and solid archives
  But we may also have the situation where there is no inventory data at
  all. In that case, we keep track of a default lilly number that matches
  the given structure
*/

class Inventory_Data {
 private:
  Source_Data<float> _liquid;
  Source_Data<float> _solid;

  IWString _default;  // if no sample data found

  // Keep track of the largest amount from each database.

  iwmaxid<float, IWString>** _maxval;

 public:
  Inventory_Data();
  ~Inventory_Data();

  int debug_print(std::ostream&) const;

  int new_mg(float, const IWString&, const IWString&);
  int new_ul(float, const IWString&, const IWString&);

  int has_data() const;

  int
  has_solid_data() const {
    return _solid.has_data();
  }

  int
  has_liquid_data() const {
    return _liquid.has_data();
  }

  int append_inventory_information(IWString&) const;

  int satisfies_minimum_amounts(float, int);

  int set_default_match_if_needed(const IWString& s);

  void amount_from_inventory(int, float, const IWString&);

  int write_tabular_file_for_eilas(const IWString&, std::ostream&) const;
};

template <typename T>
Source_Data<T>::Source_Data() {
  _amt = static_cast<T>(0);

  return;
}

template <typename T>
int
Source_Data<T>::debug_print(std::ostream& output) const {
  output << "amt " << _amt << "\n";

  return 1;
}

Inventory_Data::Inventory_Data() {
  _maxval = new iwmaxid<float, IWString>*[number_inventory_databases];
  for (int i = 0; i < number_inventory_databases; i++) {
    _maxval[i] = new iwmaxid<float, IWString>(static_cast<float>(0.0), "*");
  }
}

Inventory_Data::~Inventory_Data() {
  if (nullptr != _maxval) {
    for (int i = 0; i < number_inventory_databases; i++) {
      delete _maxval[i];
    }
    delete[] _maxval;
  }

  return;
}

int
Inventory_Data::debug_print(std::ostream& output) const {
  output << "Inventory_Data::debug_print\n";
  output << "Liquid: ";
  _liquid.debug_print(output);
  output << "Solid: ";
  _solid.debug_print(output);

  return 1;
}

int
Inventory_Data::new_mg(float mg, const IWString& lly, const IWString& dname) {
  return _solid.new_amt(mg, lly, dname);
}

int
Inventory_Data::new_ul(float ul, const IWString& lly, const IWString& dname) {
  return _liquid.new_amt(ul, lly, dname);
}

int
Inventory_Data::has_data() const {
  if (_solid.has_data()) {
    return 1;
  }

  if (_liquid.has_data()) {
    return 1;
  }

  return 0;
}

int
Inventory_Data::append_inventory_information(IWString& s) const {
  if (_solid.has_data()) {
    if (write_amounts_as_mg) {
      return _solid.append_inventory_information(s, static_cast<float>(1000.0), "mg");
    } else {
      return _solid.append_inventory_information(s, "gm");
    }
  }

  if (_liquid.has_data()) {
    return _liquid.append_inventory_information(s, "ul");
  }

  return 1;
}

int
Inventory_Data::satisfies_minimum_amounts(float sld,  // args are the global cutoffs
                                          int lq) {
  //  there were no constraints, we automatically match
  if (static_cast<float>(0.0) == sld && 0 == lq) {
    return 1;
  }

  if (static_cast<float>(0.0) == sld) {  // don't check
    ;
  } else if (_solid.has_data() && _solid.satisfies_minimum_amounts(sld)) {
    return 1;
  }

  // don't check
  if (0 == lq) {
    ;
  } else if (_liquid.has_data() && _liquid.satisfies_minimum_amounts(lq)) {
    _solid.invalidate();
    return 1;
  }

  return 0;
}

int
Inventory_Data::set_default_match_if_needed(const IWString& s) {
  if (0 == _default.length()) {
    _default = s;
  }

  return 1;
}

void
Inventory_Data::amount_from_inventory(int i, float a, const IWString& id) {
  _maxval[i]->try_this(a, id);

  return;
}

int
Inventory_Data::write_tabular_file_for_eilas(const IWString& mname,
                                             std::ostream& output) const {
  assert(nullptr != dbname);

  output << mname;

  for (int i = 0; i < number_inventory_databases; i++) {
    float m = _maxval[i]->maxval();
    const IWString& id = _maxval[i]->which_is_max();

    output << ' ' << dbname[i] << " " << id << ' ';

    if (inventory_is_liquid[i]) {
      output << m << " ul";
    } else if (write_amounts_as_mg) {
      output << static_cast<int>(m * 1000.0) << " mg";
    } else {
      output << m << " gm";
    }
  }

  output << '\n';

  return output.good();
}

template <typename T>
int
Source_Data<T>::append_inventory_information(IWString& s, T multiplier,
                                             const char* units) const {
  return _append_inventory_information(s, static_cast<T>(_amt * multiplier), units);
}

template <typename T>
int
Source_Data<T>::_append_inventory_information(IWString& s, T v, const char* units) const {
  s << _lly << ' ' << v << ' ' << units;

  if (_dbname.length()) {
    s << ' ' << _dbname;
  }

  return 1;
}

template <typename T>
int
Source_Data<T>::append_inventory_information(IWString& s, const char* units) const {
  return _append_inventory_information(s, _amt, units);
}

template <typename T>
int
Source_Data<T>::new_amt(T nam, const IWString& lly, const IWString& dname) {
  if (nam < _amt) {
    return 0;
  }

  _amt = nam;

  _lly = lly;
  _dbname = dname;

  return 1;
}

template <typename T>
int
Source_Data<T>::satisfies_minimum_amounts(T a) const {
  return _amt >= a;
}

template class Source_Data<int>;
template class Source_Data<float>;

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

  cerr << DB_VERSION_STRING << '\n';

  cerr << R"(Checks for molecules in a Berkeley DB database of unique smiles and associated inventory db's
Usage: " << prog_name << " <options> <input_file>...
  -d <dbname>    specify Berkeley database
  -I <dbname>    specify inventory database (typically gsat_inv.db)
                 liquid databases start with 'LIQUID:'
  -N <name>      name of the corresponding -I database
  -c <mg>        specify minimum amount for solid sources
  -l <ul>        specify minimum amount for liquid sources
  -H ...         specify storage conditions, enter '-H help' for info
  -b             lookup both chiral and non-chiral forms of chiral molecules
  -y             discard chirality
  -F <file>      stream for molecules found in the database
  -p             append database identifier to name
  -f             append all database identifiers to name - colon separated list
  -j             write all database identifiers on separate lines
  -U <file>      stream for molecules NOT found in the database
  -M <file>      write all molecules together with database id
  -m             write database amounts as mg rather than grams
  -u             do not display abnormal valence errors
  -R <fname>     produce tabular file with details of find rate by atom
  -e             omit the LLY- prefix prepended to found LLY numbers
  -C             when looking chiral and non chiral forms, fetch all forms present (not just first found)
  -z             remove leading zero's from LSN's reported
  -i <type>      specify input file type
  -o <type>      specify output file type(s)
  -g ...         chemical standardisation
  -A ...         aromaticity
  -v             verbose output
)";
  // clang-format on

  exit(rc);
}

static IWString
MaybeWithLlyPrepended(const IWString& mname) {
  IWString result(mname);
  if (prepend_lly_dash) {
    result << " LLY-";
  } else {
    result << ' ';
  }

  return result;
}

// Set the name of `m` to reflect both prepend_lly_dash
// and any inventory data. If there is no inventory data
// just append `token` instead.
static int
SetMoleculeName(Molecule& m,
                const Inventory_Data& invdata,
                const IWString& token) {
  IWString tmp = MaybeWithLlyPrepended(m.name());

  if (invdata.has_data()) {
    invdata.append_inventory_information(tmp);
  } else {
    tmp << token;
  }

  m.set_name(tmp);

  return 1;
}

int
MaybeWriteStreamForInDatabase(Molecule& m) {
  if (stream_for_in_database.active()) {
    return stream_for_in_database.write(m);
  }

  return 0;
}

int
MaybeWriteStreamForEverything(Molecule& m) {
  if (stream_for_everything.active()) {
    return stream_for_everything.write(m);
  }

  return 0;
}

int
MaybeWriteStreamForNotInDatabase(Molecule& m) {
  if (stream_for_not_in_database.active()) {
    return stream_for_not_in_database.write(m);
  }

  return 0;
}

static int
check_database(const IWString& key, IWString& dvalue) {
  Dbt dkey(const_cast<char*>(key.rawchars()), key.length());
  Dbt fdb;

  if (verbose > 2) {
    cerr << key << " is key\n";
  }

  if (0 != database.get(NULL, &dkey, &fdb, 0)) {
    return 0;
  }

  dvalue.set(reinterpret_cast<const char*>(fdb.get_data()),
             static_cast<int>(fdb.get_size()));

  return 1;
}

// Return the numeric value of 'x' in something that looks 
// like 'x y z'.
// `lly` is the LSN of the inventory record and is only used
// when reporting errors.
std::optional<float>
InventoryRecordToAmount(const Dbt& fromdb, const const_IWSubstring& lly) {

  const_IWSubstring amt;
  amt.set(reinterpret_cast<const char*>(fromdb.get_data()),
          static_cast<int>(fromdb.get_size()));

  amt.truncate_at_first(' ');

  float result = 0.0f;
  if (!amt.numeric_value(result) || result < 0.0f) {
    cerr << "Warning, invalid amount '" << amt << "' in database for '" << lly << "'\n";
    return std::nullopt;
  }

  return result;
}

// Fetch inventory data for `key`. Note that inventory databases do not have
// leading zero's in their keys, so we take a local copy and change it.
// `mname` is what was returned from the structure database lookup, and so
// might have leading @ characters or various other things.
static int
fetch_inventory_data(const_IWSubstring key,  // note local copy
                     const const_IWSubstring& mname,
                     Inventory_Data& invdata) {
  key.remove_leading_chars('0');
  Dbt zkey(const_cast<char*>(key.rawchars()), key.length());

  int rc = 0;

  for (int i = 0; i < number_inventory_databases; i++) {
    Dbt fromdb;
    if (0 != inventory[i]->get(NULL, &zkey, &fromdb, 0)) {
      invdata.set_default_match_if_needed(mname);
      continue;
    }

    rc++;

    std::optional<float> amt = InventoryRecordToAmount(fromdb, mname);
    if (! amt) {
      continue;
    }

    if (inventory_is_liquid[i]) {
      invdata.new_ul(*amt, mname, dbname[i]);
    } else {
      invdata.new_mg(*amt, mname, dbname[i]);
    }
  }

  return rc;
}

// If global variable append_all_matched_identifiers is set,
// append `token` to `matched_lsns`.
void
MaybeAppend(const IWString& mname,
            const_IWSubstring& token,
            IWString& matched_lsns) {
  if (! append_all_matched_identifiers) {
    return;
  }

  if (matched_lsns.empty()) {
    matched_lsns = mname;
    matched_lsns << ' ';
    if (prepend_lly_dash) {
      matched_lsns << "LSN-";
    }
  } else {
    matched_lsns << ':';
  }
  
  if (remove_leading_zeros) {
    const_IWSubstring tmp(token);
    tmp.remove_leading_chars('0');
    matched_lsns << tmp;
  } else {
    matched_lsns << token;
  }
}

// #define DEBUG_DO_WRITE_ALL_DATABASE_MATCHES

// Examine all inventory databases for the lsn's concatenated in `dvalue`.
// Append data to m.name
static int
do_write_all_database_matches(Molecule& m, const IWString& dvalue) {
  int molecule_written = 0;

  const IWString mname(m.name());  // make a copy of the name

#ifdef DEBUG_DO_WRITE_ALL_DATABASE_MATCHES
  cerr << "Parsing structure data '" << dvalue << "'\n";
#endif

  const_IWSubstring token;
  for (int i = 0; dvalue.nextword(token, i, ':');) {
    // TOKEN is our database key, LLY is what we report to the user
    const const_IWSubstring lly(token);

    // cerr << "lsn '" << lly << "'\n";
    // If this is a match only due to the loss of chirality, skip.
    if (token.starts_with('@')) {
      if (need_chirality) {
        continue;
      }
      token.remove_leading_chars(1);
    }

    Inventory_Data invdata;
    fetch_inventory_data(token, lly, invdata);

    if (!invdata.satisfies_minimum_amounts(solid_cutoff, liquid_cutoff)) {
      continue;
    }

    if (stream_for_in_database.active() || stream_for_everything.active()) {
      SetMoleculeName(m, invdata, lly);
      MaybeWriteStreamForInDatabase(m);
      MaybeWriteStreamForEverything(m);
      molecule_written = 1;
      m.set_name(mname);  // reset name so we do not accumulate inventory info.
    }
  }

  if (molecule_written) {
    return 1;
  }

  if (stream_for_not_in_database.active()) {
    stream_for_not_in_database.write(m);
  }

  return 1;
}

#ifdef NO_LONGER_NEEDED___
// Return true if either solid_cutoff or liquid_cutoff is positive.
static bool
InventoryCutoffsSpecified() {
  if (solid_cutoff > 0.0f) {
    return true;
  }

  if (liquid_cutoff > 0) {
    return true;
  }

  return false;
}
#endif

static int
process_molecule_found_in_database(Molecule& m, const IWString& dvalue) {
  molecules_in_database++;

  if (verbose) {
    atoms_in_matches[m.natoms()]++;
  }

  if (write_all_inventory_database_matches) {
    return do_write_all_database_matches(m, dvalue);
  }

  // The molecule is in the database. Parse the database record - a collection of Lilly
  // numbers, possibly with '@' characters prepended.

  Inventory_Data invdata;

  // cerr << "Token from db '" << dvalue << "'\n";

  // If we are writing all database identifiers, we want to only write those that match
  // the lookup conditions.
  IWString matched_lsns;

  if (number_inventory_databases) {
    int i = 0;
    const_IWSubstring token;
    while (dvalue.nextword(token, i, ':')) {
      // TOKEN is our database key, LLY is what we report to the user
      const_IWSubstring lly(token);

      if (token.starts_with('@')) {
        if (need_chirality) {
          continue;
        }
        token.remove_leading_chars('@');
      }

      MaybeAppend(m.name(), lly, matched_lsns);

      if (verbose > 2) {
        cerr << " lly '" << lly << "' key '" << token << "'\n";
      }

      fetch_inventory_data(token, lly, invdata);
    }
  }

  if (append_all_matched_identifiers) {
    if (matched_lsns.size() > 0) {
      m.set_name(matched_lsns);
    }
  } else if (append_database_identifier) {
    SetMoleculeName(m, invdata, dvalue);
  }

  MaybeWriteStreamForEverything(m);

  // If the amount is less than the cutoff, it is as if it isn't in the database.

  if (0 == number_inventory_databases) {
    ;
  } else if (!invdata.satisfies_minimum_amounts(solid_cutoff, liquid_cutoff)) {
    molecules_in_database_but_below_cutoff++;

    MaybeWriteStreamForNotInDatabase(m);

    return 1;
  }

  molecules_passing_constraints++;

  // We have an amount >= cutoff, or no cutoff specified

  MaybeWriteStreamForInDatabase(m);

  return 1;
}

static int
in_database(Molecule& m, Mol2Graph& mol2graph, Storage_Conditions& sc) {
  if (verbose) {
    atoms_in_input[m.natoms()]++;
  }

  IWString key;

  // not remove chirality means include it
  sc.form_key(m, mol2graph, !sc.remove_chirality(), key);

  assert(key.length() > 0);

  if (echo_key) {
    cerr << "Key " << key << " from " << m.smiles() << ' ' << m.name() << '\n';
  }

  IWString dvalue;

  const int found = check_database(key, dvalue);

  if (found && return_on_first_match) {
    return process_molecule_found_in_database(m, dvalue);
  }

  if (lookup_chiral_and_non_chiral_smiles && m.chiral_centres()) {
    sc.form_key(m, mol2graph, 0, key);

    //  cerr << "Looking up '" << key << "'\n";

    IWString s;
    if (check_database(key, s)) {
      dvalue.append_with_spacer(s, ':');

      return process_molecule_found_in_database(m, dvalue);
    }
  }

  if (found) {
    return process_molecule_found_in_database(m, dvalue);
  }

  molecules_not_in_database++;

  if (verbose > 1) {
    cerr << "Not in database\n";
  }

  if (stream_for_not_in_database.active()) {
    stream_for_not_in_database.write(m);
  }

  if (stream_for_everything.active()) {
    stream_for_everything.write(m);
  }

  return 1;
}

static void
preprocess(Molecule& m) {
  (void)chemical_standardisation.process(m);

  if (remove_directional_bonding) {
    m.revert_all_directional_bonds_to_non_directional();
  }

  m.unset_unnecessary_implicit_hydrogens_known_values();

  return;
}

static int
in_database(data_source_and_type<Molecule>& input, Mol2Graph& mol2graph,
            Storage_Conditions& sc) {
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!in_database(*m, mol2graph, sc)) {  // fatal error
      return 0;
    }
  }

  molecules_read += input.molecules_read();

  return 1;
}

static int
in_database(const char* fname, FileType input_type, Mol2Graph& mol2graph,
            Storage_Conditions& sc) {
  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (SkipSmilesErrors) {
    input.set_connection_table_errors_allowed(std::numeric_limits<int>::max());
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return in_database(input, mol2graph, sc);
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
in_database(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:o:A:g:d:abyE:F:U:pc:l:M:I:N:zfmR:LwqjeuH:CKZ");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  // TEMPORARY HACK TO BE REMOVED
  if (cl.option_present('Z')) {
    set_unique_smiles_legacy_atom_ordering(1);
  }

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

  set_global_aromaticity_type(Daylight);

  if (cl.option_present('e')) {
    prepend_lly_dash = 0;

    if (verbose) {
      cerr << "No LLY- prepended to found Lilly numbers\n";
    }
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify database via -d option\n";
    usage(8);
  }

  number_inventory_databases = cl.option_count('I');

  if (cl.option_present('z')) {
    remove_leading_zeros = 1;

    if (verbose) {
      cerr << "Leading 0's will be stripped from database records\n";
    }
  }

  if (cl.option_present('d')) {
    IWString fname;
    cl.value('d', fname);

    //  cerr << "Opening '" << fname.null_terminated_chars() << "'\n";

    int rc = database.open(NULL, fname.null_terminated_chars(), NULL, DB_UNKNOWN,
                           DB_RDONLY, 0);

    if (0 != rc) {
      cerr << "Cannot open database '" << fname << "'\n";
      database.err(rc, "");
      return 2;
    }

    if (verbose) {
      cerr << "Structure lookup from '" << fname << "'\n";
    }
  }

  // Was never needed, do not use.
  if (cl.option_present('K')) {
    cerr << "The -K option is deprecated, use '-i ICTE' instead\n";
    SkipSmilesErrors = true;
  }

  if (cl.option_present('a') && cl.option_present('b')) {
    cerr << "Sorry, cannot use both -a and -b options\n";
    usage(8);
  }

  Storage_Conditions sc;

  Mol2Graph mol2graph;

  if (cl.option_present('H')) {
    if (!sc.initialise(cl, 'H', mol2graph, verbose)) {
      cerr << "Cannot initialise storage conditions (-H)\n";
      return 1;
    }
  }

  if (cl.option_present('L')) {
    sc.set_reduce_to_largest_fragment(1);

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('a')) {
    cerr << "Warning the -a option is obsolete, use the -H option for graph/tautomer lookups\n";
#ifdef DASH_A_OPTION_NOW_OBSOLETE
    do_graph_lookup = 1;
    if (verbose) {
      cerr << "Lookup will be based on a graph lookup\n";
    }

    if (cl.option_present('w')) {
      use_aromatic_distinguishing_mf_in_tautomer = 1;

      if (verbose) {
        cerr << "Will use aromatic distinguishing formula for tautomer matches\n";
      }
    }
    if (cl.option_present('q')) {
      set_exclude_triple_bonds_from_graph_reduction(1);

      if (verbose) {
        cerr << "Triple bonds not reduced during graph formation\n";
      }
    }
#endif  // DASH_A_OPTION_NOW_OBSOLETE
  } else if (cl.option_present('y')) {
    need_chirality = 0;
    if (verbose) {
      cerr << "Chirality will be excluded from lookups\n";
    }
  } else if (cl.option_present('b')) {
    need_chirality = 1;
    lookup_chiral_and_non_chiral_smiles = 1;

    if (verbose) {
      cerr << "Will look up both chiral and non-chiral smiles\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', solid_cutoff) || solid_cutoff < 0) {
      cerr << "The solid cutoff option (-c) must be followed by a valid whole number "
              "(mg)\n";
      usage(14);
    }

    if (verbose) {
      cerr << "Solid cutoff is " << solid_cutoff << " mg\n";
    }

    solid_cutoff = solid_cutoff / 1000.0F;  // convert to grams
  }

  if (cl.option_present('C')) {
    return_on_first_match = 0;
    need_chirality = 1;
    lookup_chiral_and_non_chiral_smiles = 1;
    if (verbose) {
      cerr << "Will fetch all instances found in the database\n";
    }
  }

  if (cl.option_present('l')) {
    if (!cl.value('l', liquid_cutoff) || liquid_cutoff < 0) {
      cerr << "The liquid cutoff option (-l) must be followed by a valid whole number "
              "(mg)\n";
      usage(14);
    }

    if (verbose) {
      cerr << "Liquid cutoff is " << liquid_cutoff << " mg\n";
    }
  }

  if (cl.option_present('u')) {
    set_display_abnormal_valence_messages(0);
    if (verbose) {
      cerr << "Will suppress messages about invalid valences\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (0 == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (!cl.option_present('F') && !cl.option_present('U') && !cl.option_present('M') &&
      !cl.option_present('T')) {
    cerr << "Must specify -F, -M, -U or -T for some output\n";
    usage(18);
  }

  if (cl.option_present('T') &&
      (cl.option_present('F') || cl.option_present('M') || cl.option_present('U'))) {
    cerr << "The -T option is incompatible with the -F, -M and -U options\n";
    usage(6);
  }

  FileType output_type = FILE_TYPE_INVALID;
  if (!cl.option_present('o')) {
    if (verbose) {
      cerr << "Output type defaults to SMI\n";
    }
    output_type = FILE_TYPE_SMI;
  }

  if (cl.option_present('F')) {
    if (!handle_file_opening(cl, 'F', output_type, stream_for_in_database)) {
      return 13;
    }
  }

  if (cl.option_present('p') && cl.option_present('f')) {
    cerr << "Sorry, the -p and -f options are mutually incompatible\n";
    usage(4);
  }

  if (cl.option_present('j')) {
    write_all_inventory_database_matches = 1;

    if (verbose) {
      cerr << "Will write all inventory matches\n";
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
  } else if (cl.option_present('f')) {
    append_all_matched_identifiers = 1;

    if (verbose) {
      cerr << "Will append all matching Lilly numbers to the output\n";
    }
  }

  if (cl.option_present('m')) {
    write_amounts_as_mg = 1;

    if (verbose) {
      cerr << "Will write database amounts as MG rather than Grams\n";
    }

    if (!append_database_identifier) {
      append_database_identifier = 1;
    }
  }

  if (cl.option_present('U')) {
    if (!handle_file_opening(cl, 'U', output_type, stream_for_not_in_database)) {
      return 14;
    }
  }

  if (cl.option_present('M')) {
    if (!handle_file_opening(cl, 'M', output_type, stream_for_everything)) {
      return 15;
    }
  }

  if (cl.option_present('R')) {
    fname_for_per_natoms_lookup_report = cl.option_value('R');

    if (verbose) {
      cerr << "Report on find rate as a function of atom count written to '"
           << fname_for_per_natoms_lookup_report << "'\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (number_inventory_databases) {
    inventory = new Db*[number_inventory_databases];
    dbname = new IWString[number_inventory_databases];
    inventory_is_liquid = new int[number_inventory_databases];

    for (int i = 0; i < number_inventory_databases; i++) {
      IWString fname = cl.string_value('I', i);

      if (fname.starts_with("LIQUID:")) {
        inventory_is_liquid[i] = 1;
        fname.remove_leading_chars(7);
      } else {
        inventory_is_liquid[i] = 0;
      }

      inventory[i] = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

      int rc = inventory[i]->open(NULL, (char*)fname.null_terminated_chars(), NULL,
                                  DB_UNKNOWN, DB_RDONLY, 0);

      if (0 != rc) {
        cerr << "Cannot open database '" << fname << "' '";
        inventory[i]->err(rc, "");
        return i + 1;
      }
    }
  }

  MDL_File_Supporting_Material* mdlfsm = global_default_MDL_File_Supporting_Material();
  mdlfsm->set_ignore_unrecognised_mdl_m_records(1);

  int n = cl.option_count('N');
  if (0 == n) {
    ;
  } else if (n != number_inventory_databases) {
    cerr << "If the -N option is used, there must be the same number of -N as -I "
            "options\n";
    usage(5);
  } else {
    for (int i = 0; i < n; i++) {
      cl.value('N', dbname[i], i);
    }
  }

  set_display_abnormal_valence_messages(0);  // error messages not informative here
  set_display_strange_chemistry_messages(0);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!in_database(cl[i], input_type, mol2graph, sc)) {
      rc = i + 1;
      break;
    }
  }

  for (int i = 0; i < number_inventory_databases; i++) {
    inventory[i]->close(0);
    delete inventory[i];
  }

  delete [] inventory;
  delete [] dbname;
  delete [] inventory_is_liquid;

  if (verbose) {
    cerr << "Examined " << molecules_read << " molecules\n";
    cerr << "Found " << molecules_in_database << " in database, "
         << molecules_not_in_database << " not in database\n";

    if (solid_cutoff > static_cast<float>(0.0) || liquid_cutoff > 0) {
      cerr << molecules_in_database_but_below_cutoff
           << " molecules in the database but with amounts less than cutoff(s)\n";
    }

    cerr << molecules_passing_constraints << " molecules passed the constraints\n";

    write_per_atom_lookup_rates(atoms_in_input, atoms_in_matches, cerr);
  }

  if (malformed_database_records) {
    cerr << malformed_database_records << " erroneous database records\n";
  }

  if (verbose && negative_amounts_ignored) {
    cerr << "Ignored " << negative_amounts_ignored
         << " inventory records with negative amounts\n";
  }

  if (! fname_for_per_natoms_lookup_report.empty()) {
    IWString_and_File_Descriptor output;

    if (!output.open(fname_for_per_natoms_lookup_report)) {
      cerr << "Cannot open file for by atom find rate '"
           << fname_for_per_natoms_lookup_report << "'n";
      return 3;
    }

    set_default_iwstring_float_concatenation_precision(3);

    write_per_atom_lookup_rates(atoms_in_input, atoms_in_matches, output);
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = in_database(argc, argv);

  return rc;
}
