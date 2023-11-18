// Look in a BerkeleyDB inventory database for molecules with enough sample.

#include <stdlib.h>

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwre2.h"

#include "db_cxx.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static float solid_cutoff = static_cast<float>(0.0);

static float liquid_cutoff = static_cast<float>(0.0);

static float solid_micromoles_needed = static_cast<float>(0.0);

static int records_read = 0;

static int records_with_enough_inventory = 0;

static int records_with_not_enough_inventory = 0;

static int records_with_no_inventory_data = 0;

static int chop_solid_amounts_to_milligrams = 1;

static IWString_and_File_Descriptor stream_for_not_enough_inventory;

static int append_amount_to_output = 0;

static int append_molecular_weight = 0;

//static int write_amount_as_initial_lookup = 0;

static int identifier_column = -1;

static int break_after_finding_a_match = 1;

static std::unique_ptr<re2::RE2> id_rx;

static int input_is_tdt = 0;

static IWString identifier_tag("PCN<");

static int auto_strip_leading_zeros = 1;

class Inventory_Database {
 protected:
  Db _database;

  IWString _name;  // "gsat", "flexstore" "wetlab" etc..

  int _samples_looked_up;

  int _samples_found_in_database;

  int _samples_found_in_database_with_enough_inventory;

  float _cutoff;

  IWString _units;  // gm, ul, ...

  //  protected functions

  int _common_lookup(Dbt&, int&, IWString&);
  int _common_lookup(Dbt& dkey, int& found, IWString& retrieved,
                     const_IWSubstring& amount_as_string, float& amt);

 public:
  Inventory_Database();
  virtual ~Inventory_Database();

  virtual int report(std::ostream&) const;

  int build(IWString&, const IWString&);

  const IWString&
  id() const {
    return _name;
  }

  const IWString&
  units() const {
    return _units;
  }

  virtual int lookup(Dbt&, int&, IWString&) = 0;
};

Inventory_Database::Inventory_Database() : _database(NULL, DB_CXX_NO_EXCEPTIONS) {
  _cutoff = static_cast<float>(0.0);

  _samples_looked_up = 0;

  _samples_found_in_database = 0;

  _samples_found_in_database_with_enough_inventory = 0;

  return;
}

Inventory_Database::~Inventory_Database() {
  _database.close(0);

  return;
}

int
Inventory_Database::build(IWString& dbname, const IWString& dname) {
  int rc = _database.open(NULL, dbname.null_terminated_chars(), NULL, DB_UNKNOWN,
                          DB_RDONLY, 0);

  if (0 != rc) {
    cerr << "Cannot open database '" << dbname << "' '";
    _database.err(rc, "");
    return 0;
  }

  _name = dname;

  return 1;
}

int
Inventory_Database::report(std::ostream& os) const {
  os << "Database '" << _name << "' looked up " << _samples_looked_up << " samples.\n";
  os << _samples_found_in_database << " samples found in database\n";
  os << _samples_found_in_database_with_enough_inventory << " samples with >= " << _cutoff
     << ' ' << _units << '\n';

  return os.good();
}

class Liquid_Inventory_Database : public Inventory_Database {
 private:
 public:
  Liquid_Inventory_Database();

  int lookup(Dbt&, int&, IWString&);
};

Liquid_Inventory_Database::Liquid_Inventory_Database() {
  _cutoff = liquid_cutoff;

  _units = "ul";

  return;
}

class Echo_Inventory_Database : public Inventory_Database {
 private:
 public:
  Echo_Inventory_Database();

  int lookup(Dbt&, int&, IWString&);
};

Echo_Inventory_Database::Echo_Inventory_Database() {
  _cutoff = 0.0f;

  _units = "ul";

  return;
}

class Solid_Inventory_Database : public Inventory_Database {
 private:
  int _need_to_convert_to_mg;

  float _micromoles_needed;

  //  private functions

  int _enough_sample(float amt, const IWString& retrieved) const;

 public:
  Solid_Inventory_Database();

  int report(std::ostream&) const;

  void
  set_need_to_convert_to_mg(const int s) {
    _need_to_convert_to_mg = s;
  }

  int lookup(Dbt&, int&, IWString&);
};

Solid_Inventory_Database::Solid_Inventory_Database() {
  _micromoles_needed = solid_micromoles_needed;

  _cutoff = solid_cutoff;

  _need_to_convert_to_mg = 1;
  if (_need_to_convert_to_mg) {
    _cutoff = _cutoff / 1000.0f;
  }

  _units = "mg";

  return;
}

int
Inventory_Database::_common_lookup(Dbt& dkey, int& found, IWString& retrieved) {
  _samples_looked_up++;

#ifdef DEBUG_LOOKUP
  cerr << "Looking for " << cerr.write(dkey.dptr, dkey.dsize) << '\n';
#endif

  Dbt zdata;

#define BUFFER_SIZE 512
  char buffer[BUFFER_SIZE];  // much larger than any inventory data stored

  zdata.set_data(&buffer);
  zdata.set_size(BUFFER_SIZE);
  zdata.set_ulen(BUFFER_SIZE);
  zdata.set_flags(DB_DBT_USERMEM);

  if (0 != _database.get(NULL, &dkey, &zdata, 0)) {
    return 0;
  }

  retrieved.set(reinterpret_cast<const char*>(zdata.get_data()),
                static_cast<int>(zdata.get_size()));

  found++;

  _samples_found_in_database++;

  return 1;
}

int
Inventory_Database::_common_lookup(Dbt& dkey, int& found, IWString& retrieved,
                                   const_IWSubstring& amount_as_string, float& amt) {
  if (!_common_lookup(dkey, found, retrieved)) {
    return 0;
  }

  amount_as_string = retrieved;

  amount_as_string.truncate_at_first(' ');

  if (!amount_as_string.numeric_value(amt) || amt < static_cast<float>(0.0)) {
    cerr << "Yipes, invalid amount '" << retrieved << "', ignored\n";
    return 0;
  }

  if (_cutoff > static_cast<float>(0.0) && amt < _cutoff) {
    return 0;
  }

  _samples_found_in_database_with_enough_inventory++;

  return 1;
}

int
Solid_Inventory_Database::report(std::ostream& os) const {
  os << "Database '" << _name << "' looked up " << _samples_looked_up << " samples.\n";
  os << _samples_found_in_database << " samples found in database\n";
  os << _samples_found_in_database_with_enough_inventory << " samples with >= ";
  if (_cutoff > static_cast<float>(0.0)) {
    os << _cutoff << ' ' << _units << '\n';
  } else {
    os << _micromoles_needed << " uM\n";
  }

  return os.good();
}

int
Liquid_Inventory_Database::lookup(Dbt& dkey, int& found, IWString& output_buffer) {
  IWString retrieved;
  const_IWSubstring amount_as_string;
  float amt;

  if (!_common_lookup(dkey, found, retrieved, amount_as_string, amt)) {
    return 0;
  }

  if (append_amount_to_output) {
    output_buffer << ' ' << _name << ' ' << amount_as_string << ' ' << _units;
  }

  return 1;
}

int
Solid_Inventory_Database::_enough_sample(float amt, const IWString& retrieved) const {
  if (1 == retrieved.nwords()) {
    cerr << "Solid_Inventory_Database::_enough_sample:no AMW '" << retrieved << "'\n";
    return 0;
  }

  const_IWSubstring tmp;
  retrieved.word(1, tmp);

  float amw;

  if (!tmp.numeric_value(amw)) {
    cerr << "Solid_Inventory_Database::_enough_sample:invalid AMW '" << retrieved
         << "'\n";
    return 0;
  }

  const float amount_per_micromole = amw / 1000.0;

  float micromoles_present = amt / amount_per_micromole;

  return micromoles_present >= _micromoles_needed;
}

int
Solid_Inventory_Database::lookup(Dbt& dkey, int& found, IWString& output_buffer) {
  IWString retrieved;
  const_IWSubstring amount_as_string;
  float amt;

  // This works because _cutoff will have been set to zero

  if (!_common_lookup(dkey, found, retrieved, amount_as_string, amt)) {
    return 0;
  }

#ifdef DEBUG_LOOKUP
  cerr << "Looking for " << cerr.write(dkey.dptr, dkey.dsize) << '\n';
#endif

  if (_need_to_convert_to_mg) {
    amt *= 1000.0f;
  }

  if (chop_solid_amounts_to_milligrams) {
    amt = static_cast<float>(static_cast<int>(amt));
  }

  if (_cutoff >= 0.0f) {  // already testing in _common_lookup
    ;
  } else if (!_enough_sample(amt, retrieved)) {
    return 0;
  }

  if (append_amount_to_output) {
    output_buffer << ' ' << _name << ' ' << amt << ' ' << _units;
  }

  if (append_molecular_weight) {
    retrieved.remove_leading_words(1);
    if (!output_buffer.ends_with(' ')) {
      output_buffer << ' ';
    }
    output_buffer << retrieved;
  }

  return 1;
}

/*
  Exactly the same as a Liquid inventory database, should
*/

int
Echo_Inventory_Database::lookup(Dbt& dkey, int& found, IWString& output_buffer) {
  IWString retrieved;
  const_IWSubstring amount_as_string;
  float amt;

  if (!_common_lookup(dkey, found, retrieved, amount_as_string, amt)) {
    return 0;
  }

  if (append_amount_to_output) {
    output_buffer << ' ' << _name << ' ' << amount_as_string << ' ' << _units;
  }

  return 1;
}

static int number_databases = 0;
static Inventory_Database** database = NULL;

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
  cerr << "Checks against inventory databases\n";
  cerr << " -d L:<dbname>  liquid database\n";
  cerr << " -d S:<dbname>  solid database\n";
  cerr << " -N name        name for -d database(s)\n";
  cerr << " -c <mg>        cutoff in mg\n";
  cerr << " -m <ul>        cutoff in ul\n";
  cerr << " -C <col>       identifier column - use -C 2 for a smiles file\n";
  cerr << " -B <fname>     write records with no or not enough inventory to <fname>\n";
  cerr << " -a             append the inventory amount and source to each record written\n";
  cerr << " -k             look in all databases, even after a match is found\n";
  cerr << " -Z             do NOT automatically strip leading zero's from identifiers\n";
  cerr << " -u <uM>        minimum number of micromoles in Solid inventory\n";
  cerr << " -y             do NOT chop fractional milligrams to whole numbers\n";
  cerr << " -g             input is a fingerprint file\n";
  cerr << " -v             verbose output\n";
  cerr << "Note to maintainers:amounts in the DB are expected to be in grams\n";
  // clang-format on

  exit(rc);
}

static int
enough_inventory_3(const const_IWSubstring& lly, int& found_match,
                   IWString& amount_data) {
  found_match = 0;

  if (id_rx && !iwre2::RE2PartialMatch(lly, *id_rx)) {
    cerr << "Invalid identifier '" << lly << "' does not match '" << id_rx->pattern()
         << "'\n";
    return 0;
  }

  Dbt zkey;
  zkey.set_data((void*)(lly.rawchars()));  // dangerous cast, but should be OK
  zkey.set_size(lly.length());

  int nfound = 0;  // how many of the databases have a key for LLY

  for (int i = 0; i < number_databases; i++) {
    if (database[i]->lookup(zkey, nfound, amount_data)) {
      found_match++;

      if (break_after_finding_a_match) {
        break;
      }
    }
  }

  if (found_match) {  // found in at least one database with enough inventory
    records_with_enough_inventory++;
  } else if (0 == nfound) {  // not found in any of the databases
    records_with_no_inventory_data++;
  } else {
    records_with_not_enough_inventory++;
  }

  return 1;
}

static int
enough_inventory_2(const const_IWSubstring& buffer, const const_IWSubstring& lly,
                   IWString_and_File_Descriptor& output) {
  int found_match;

  IWString amount_data;

  if (!lly.starts_with('0')) {
    if (!enough_inventory_3(lly, found_match, amount_data)) {  // ghastly error
      return 0;
    }
  } else {
    const_IWSubstring mylly(lly);
    if (auto_strip_leading_zeros) {
      mylly.remove_leading_chars('0');
    }

    if (!enough_inventory_3(mylly, found_match, amount_data)) {
      return 0;
    }
  }

  if (found_match) {
    output << buffer;
    if (append_amount_to_output) {
      output << amount_data;
    }
    output << '\n';

    output.write_if_buffer_holds_more_than(16384);
  } else if (stream_for_not_enough_inventory.is_open()) {
    stream_for_not_enough_inventory << buffer << '\n';
    stream_for_not_enough_inventory.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
enough_inventory(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output) {
  if (identifier_column < 0) {
    return enough_inventory_2(buffer, buffer, output);
  }

  const_IWSubstring string_lly;

  if (!buffer.word(identifier_column, string_lly)) {
    cerr << "Cannot extract column " << identifier_column << " from record\n";
    return 0;
  }

  return enough_inventory_2(buffer, string_lly, output);
}

static int
enough_inventory_tdt(const IW_TDT& tdt, IWString_and_File_Descriptor& output) {
  IWString id;
  if (!tdt.dataitem_value(identifier_tag, id)) {
    cerr << "Cannot extract identifier from TDT\n";
    return 0;
  }

  id.truncate_at_first(' ');

  IWString notused;

  int found_match;
  if (!enough_inventory_3(id, found_match, notused)) {
    return 0;
  }

  if (found_match) {
    output << tdt;
    output.write_if_buffer_holds_more_than(16384);
  } else if (stream_for_not_enough_inventory.is_open()) {
    stream_for_not_enough_inventory << tdt;
    stream_for_not_enough_inventory.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
enough_inventory_tdt(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    records_read++;

    if (!enough_inventory_tdt(tdt, output)) {  // fatal error
      return 0;
    }
  }

  return output.good();
}

static int
enough_inventory(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  input.set_translate_tabs(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    records_read++;

    if (!enough_inventory(buffer, output)) {  // fatal error
      return 0;
    }
  }

  return output.good();
}

static int
enough_inventory(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (input_is_tdt) {
    return enough_inventory_tdt(input, output);
  }

  input.set_dos(1);

  return enough_inventory(input, output);
}

static int
enough_inventory(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:N:c:m:C:B:akR:wu:yZ");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('c') && !cl.option_present('m') && !cl.option_present('u')) {
    cerr << "Must specify amount cutoff in mg via the -c, -m or -u options\n";
    usage(2);
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', identifier_column) || identifier_column < 1) {
      cerr << "Invalid value for identifier column (-C option)\n";
      usage(6);
    }

    if (verbose) {
      cerr << "Identifiers in column " << identifier_column << '\n';
    }

    identifier_column--;
  } else {
    identifier_column = 0;
  }

  if (cl.option_present('a')) {
    append_amount_to_output = 1;

    if (verbose) {
      cerr << "Will append the amount found to the output\n";
    }
  }

  if (cl.option_present('w')) {
    append_molecular_weight = 1;
    append_amount_to_output = 1;

    if (verbose) {
      cerr << "Will append the molecular weight to the output\n";
    }
  }

  if (cl.option_present('k')) {
    break_after_finding_a_match = 0;

    if (verbose) {
      cerr << "Will search all databases - even after a match is found\n";
    }
  }

  if (cl.option_present('c') && cl.option_present('u')) {
    cerr << "Sorry, cannot use both the -c and -u options\n";
    usage(2);
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', solid_cutoff) || solid_cutoff < 0.0) {
      cerr << "Invalid solid cutoff specification (-c option)\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will only select records with " << solid_cutoff
           << " or more mg. of inventory\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', liquid_cutoff) || liquid_cutoff < 0.0) {
      cerr << "Invalid liquid cutoff specification (-m option)\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will only select records with " << liquid_cutoff
           << " or more liquid units of inventory\n";
    }
  }

  if (cl.option_present('u')) {
    if (!cl.value('u', solid_micromoles_needed) || solid_micromoles_needed < 0.0) {
      cerr << "Invalid micromolar cutoff specification (-u option)\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will only select records with " << solid_micromoles_needed
           << " or more micromoles of solid inventory\n";
    }
  }

  if (cl.option_present('y')) {
    chop_solid_amounts_to_milligrams = 0;
    if (verbose) {
      cerr << "Solid amounts NOT truncated to whole milligrams\n";
    }
  }

  if (cl.option_present('Z')) {
    auto_strip_leading_zeros = 0;

    if (verbose) {
      cerr << "Will NOT strip leading zero's from identifiers\n";
    }
  }

  if (cl.option_present('R')) {
    const_IWSubstring r = cl.string_value('R');

    if (!iwre2::RE2Reset(id_rx, r)) {
      cerr << "Invalid identifier regular expression '" << r << "'\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Identifiers must match '" << id_rx->pattern() << "'\n";
    }
  }

  number_databases = cl.option_count('d');

  if (0 == number_databases) {
    cerr << "Must specify inventory database(s) via the -d option\n";
    usage(5);
  }

  int names_specified = cl.option_count('N');

  if (0 == names_specified) {  // no names
    ;
  } else if (names_specified != number_databases) {
    cerr << "Must specify as many database names via the -N option as databases via the "
            "-d option\n";
    usage(5);
  }

  database = new Inventory_Database*[number_databases];

  int number_liquid_databases = 0;
  int number_solid_databases = 0;
  int number_echo_databases = 0;

  if (cl.option_present('d')) {
    int i = 0;
    IWString dbname;
    while (cl.value('d', dbname, i)) {
      if (dbname.starts_with("L:")) {
        number_liquid_databases++;
        database[i] = new Liquid_Inventory_Database;
      } else if (dbname.starts_with("S:")) {
        number_solid_databases++;
        database[i] = new Solid_Inventory_Database;
      } else if (dbname.starts_with("E:")) {
        number_echo_databases++;
        database[i] = new Echo_Inventory_Database;
      } else {
        cerr << "Unrecognised database type '" << dbname << "'\n";
        usage(4);
      }

      dbname.remove_leading_chars(2);

      IWString dtype;

      if (names_specified) {
        cl.value('N', dtype, i);
      }

      if (!database[i]->build(dbname, dtype)) {
        cerr << "Bad news, cannot initialise database '" << dbname << "'\n";
        return 6;
      }

      i++;
    }
  }

  if (number_solid_databases && !(cl.option_present('c') || cl.option_present('u'))) {
    cerr << "Must specify a cutoff in mg via the -c option\n";
    usage(5);
  }

  if (number_liquid_databases && !cl.option_present('m')) {
    cerr << "Must specify a cutoff in liquid via the -m option\n";
    usage(5);
  }

  if (0 == number_solid_databases && cl.option_present('c')) {
    cerr << "Solid inventory amount cutoff (-c) but no solid inventory databases\n";
  }

  if (0 == number_solid_databases && cl.option_present('u')) {
    cerr << "Solid inventory micromolar cutoff (-u) but no solid inventory databases\n";
  }

  if (0 == number_liquid_databases && cl.option_present('m')) {
    cerr << "Liquid inventory amount cutoff (-m) but no liquid inventory databases\n";
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('B')) {
    IWString fname = cl.string_value('B');

    stream_for_not_enough_inventory.open(fname.null_terminated_chars());

    if (!stream_for_not_enough_inventory.good()) {
      cerr << "Sorry, cannot open file for not enough inventory '" << fname << "'\n";
      return 6;
    }

    if (verbose) {
      cerr << "Molecules with not enough inventory written to '" << fname << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!enough_inventory(cl[i], output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << records_read << " records\n";
    if (records_with_no_inventory_data) {
      cerr << records_with_no_inventory_data << " with no inventory data\n";
    }
    if (records_with_not_enough_inventory) {
      cerr << records_with_not_enough_inventory << " with less than required amount\n";
    }
  }

  for (int i = 0; i < number_databases; i++) {
    if (verbose) {
      database[i]->report(cerr);
    }

    delete database[i];
  }

  delete[] database;

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = enough_inventory(argc, argv);

  return rc;
}
