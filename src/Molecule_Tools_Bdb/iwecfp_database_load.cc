// Create database of extended connectivity fingerprints
// Profile a collection of molecules, and store the EC fingerprints
// encountered. That database is then used by iwecfp_database_lookup
// for making assessments of synthetic precedent.

#include <iostream>
#include <limits>
#include <memory>

#ifdef __WIN32__
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

#include <assert.h>

#include "db_cxx.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "iwecfp_database.h"

using std::cerr;

using iwecfp_database::Bit_Produced;
using iwecfp_database::Count_Radius;
using iwecfp_database::DBKey;
using iwecfp_database::Fingerprint_Characteristics;
using iwecfp_database::Set_of_Bits;

static Chemical_Standardisation chemical_standardisation;

static int verbose = 0;

static int molecules_read = 0;

static int only_set_bits_for_max_radius_shell = 0;

static Accumulator_Int<int> nbits_acc;

static int only_start_at_chiral_atoms = 0;

/*
  There are two ways we can update the global fingerprint.
  We can count all the bits present in each molecule, or we can
  just record presence of the bit in yet another molecule
*/

static int update_global_fingerprint_presence_only = 1;

static Db database(NULL, DB_CXX_NO_EXCEPTIONS);

static Report_Progress report_progress;

static int remove_leading_zeros = 0;

// static int db_cache_size = 33554432;

static int store_flag = DB_NOOVERWRITE;

static int duplicate_keys_not_stored = 0;

static int databases_contain_examples = 1;

static int need_to_swap_bytes = 0;  // we don't actually do this

static int bit_collisions = 0;
static int intra_molecular_bit_collisions = 0;

static int free_memory = 1;

class Count_Example
{
 private:
  int _number_examples;

  const int _radius;

  IWString _first_example_smiles;
  IWString _first_example_name;

  //  To help with collisions

  const int _atom_constant_centre_atom;

  //  private functions

  int _append_new_data_to_existing(Dbt& fromdb, IWString& to_store) const;
  int _store_bit_no_examples(const DBKey& dbkey, Db& bit_database) const;
  int _do_store(Db& bit_database, Dbt& key, Dbt& zdata, int store_flag) const;

 public:
  Count_Example(Molecule&, atom_number_t centre_atom, int r, int ac);

  void extra()
  {
    _number_examples++;
  }

  int radius() const
  {
    return _radius;
  }

  int atom_constant_centre_atom() const
  {
    return _atom_constant_centre_atom;
  }

  void set_number_examples(int s)
  {
    _number_examples = s;
  }

  void set_first_example_smiles(const const_IWSubstring& s)
  {
    _first_example_smiles = s;
  }

  int number_examples() const
  {
    return _number_examples;
  }

  //  const IWString & first_example () const {return _first_example;}

  int store_bit(const DBKey&, Db&);

  int do_write(const DBKey&, IWString_and_File_Descriptor&) const;
};

Count_Example::Count_Example(Molecule& m, atom_number_t centre_atom, int r, int ac)
    : _radius(r), _atom_constant_centre_atom(ac)
{
  _number_examples = 1;

  if (!databases_contain_examples) {
    return;
  }

  _first_example_smiles = m.smiles();  // store non isotope form in selimsteg database

  // if (0 == r)    // because isotope zero doesn't work
  //   r = 1;

  // radius 0 does not do anything.
  // That's OK, doing this off by one thing becomes to confusing.
  // Radius 0 should always be obvious
  m.set_isotope(centre_atom, r);

  _first_example_name << m.smiles() << ' ' << m.name();

  m.set_isotope(centre_atom, 0);

  return;
}

int
Count_Example::_store_bit_no_examples(const DBKey& dbkey, Db& bit_database) const
{
  Dbt key((void*)(&dbkey), sizeof(dbkey));  // loss of const OK

  int count_to_store = _number_examples;

  if (DB_NOOVERWRITE == store_flag) {
    Dbt fromdb;
    if (0 == bit_database.get(NULL, &key, &fromdb, 0)) {
      if (sizeof(Count_Radius) != fromdb.get_size()) {
        cerr << "Count_Example::_store_bit_no_examples:bad data in db, size "
             << fromdb.get_size() << '\n';
        return 0;
      }
      const Count_Radius* cr = reinterpret_cast<const Count_Radius*>(fromdb.get_data());
      count_to_store += cr->_count;
    }
  }

  Count_Radius cr;
  cr._count = count_to_store;
  cr._radius = _radius;

  Dbt zdata(&cr, sizeof(cr));

  return _do_store(bit_database, key, zdata, store_flag);
}

/*
  We are attempting to do a store, and we've found that ther is already
  data in the database. We increment the number
*/

int
Count_Example::_append_new_data_to_existing(Dbt& fromdb, IWString& to_store) const
{
  const_IWSubstring s((const char*)fromdb.get_data(), fromdb.get_size());

  const_IWSubstring string_count, other_stuff;

  s.split(string_count, ' ', other_stuff);

  int existing_count;

  if (!string_count.numeric_value(existing_count) || existing_count < 0) {
    cerr << "Count_Example::_append_new_data_to_existing:invalid stored data '" << s
         << "'\n";
    return 0;
  }

  to_store << (_number_examples + existing_count) << ' ' << other_stuff;

  return 1;
}

int
Count_Example::_do_store(Db& bit_database, Dbt& key, Dbt& zdata, int store_flag) const
{
  int rc = bit_database.put(NULL, &key, &zdata, store_flag);

  if (0 == rc) {
    return 1;
  }

  if (DB_KEYEXIST == rc) {
    duplicate_keys_not_stored++;
    return 1;
  }

  bit_database.err(rc, "Cannot store bit");
  return 0;
}

// #define DEBUG_STORE_BIT

int
Count_Example::store_bit(const DBKey& dbkey, Db& bit_database)
{
  assert(_radius == dbkey._radius);

  if (!databases_contain_examples) {
    return _store_bit_no_examples(dbkey, bit_database);
  }

  Dbt key((void*)(&dbkey), sizeof(dbkey));  // loss of const OK

  IWString to_store;

  // If we are not overwriting, we need to fetch anything already stored and
  // append our new data

  if (DB_NOOVERWRITE == store_flag) {
    Dbt fromdb;
    int rc = bit_database.get(NULL, &key, &fromdb, 0);
    if (0 == rc) {
      _append_new_data_to_existing(fromdb, to_store);
    }
  }

  if (0 == to_store.length()) {
    to_store << _number_examples << ' ' << _radius << ' ' << _first_example_name;
  }

  Dbt zdata;
  zdata.set_data(const_cast<char*>(to_store.rawchars()));
  zdata.set_size(to_store.size());

#ifdef DEBUG_STORE_BIT
  cerr << "Bit " << dbkey._bit << " storing '" << to_store << "'\n";

  cerr << "Storing ";
  debug_print_key_components(dbkey, cerr);
#endif

  return _do_store(bit_database, key, zdata, store_flag);
}

int
Count_Example::do_write(const DBKey& dbkey, IWString_and_File_Descriptor& output) const
{
  output << dbkey._bit << ' ' << dbkey._acca << ' ' << static_cast<int>(dbkey._radius);

  if (databases_contain_examples) {
    output << ' ' << _first_example_smiles << ' ' << _first_example_name;
  }

  output << '\n';

  return 1;
}

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << DB_VERSION_STRING << '\n';
  cerr << "Compute the Extended Connectivity fingerprints for molecules and store\n";
  cerr << "Warning, does not do updates, use 'iwbdb_cat.sh -p' for updates\n";
  cerr << "  -r <len>       min shell width for writing a fingerprint\n";
  cerr << "  -R <length>    set the maximum step for the connected shell\n";
  cerr << "  -d <dbname>    database(s) to use\n";
//cerr << "  -e <env>       specify database envirnment to use\n";
  cerr << "  -d <dbname>    database to populate when profiling\n";
  cerr << "  -p <number>    report progress every <number> molecules processed\n";
  cerr << "  -P <atype>     atom type to use (suggest SFX)\n";
  cerr << "  -e             record count rather than presence in each molecule\n";
//cerr << "  -c             building a chiral database, only start shells at chiral atoms\n";
  cerr << "  -y <type>      access type (btree, hash, recno, queue)\n";
  cerr << "  -z             remove leading zeros from identifiers (saves space)\n";
  cerr << "  -M ...         miscellaneous options, enter '-M help' for info\n";
  cerr << "  -i <type>      input type\n";
  //  (void) display_standard_aromaticity_options (cerr);
  //  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  //  (void) display_standard_sparse_fingerprint_options (cerr, 'F');
  cerr << "  -E ...         standard element options\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
preprocess_molecule(Molecule& m)
{
  m.remove_all(1);

  m.reduce_to_largest_fragment();  // always reduce to largest fragment

  if (!only_start_at_chiral_atoms) {
    m.remove_all_chiral_centres();  // not useful here
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (remove_leading_zeros) {
    const_IWSubstring s = m.name();
    if (s.starts_with('0')) {
      s.remove_leading_chars('0');
      m.set_name(s);
    }
  }

  return 1;
}

// Accumulates bits in RAM and then writes to be BerkeleyDB database.
class Global_Fingerprint
{
 private:
  typedef std::unordered_map<DBKey, Count_Example*, iwecfp_database::IWdbkeyHash>
      bce_hash;

  bce_hash _bits;

  //  private functions

  int
  _increment_count(const DBKey&, Molecule&, atom_number_t, int acca);

 public:
  Global_Fingerprint();
  ~Global_Fingerprint();

  int update(const Set_of_Bits<Bit_Produced>& sob, Molecule& m);

  int do_store(Db&, const int max_shell_radius) const;

  int do_write(const char* fname) const;
  int do_write(IWString_and_File_Descriptor& output) const;

  int size() const
  {
    return _bits.size();
  }
};

Global_Fingerprint::Global_Fingerprint()
{
}

Global_Fingerprint::~Global_Fingerprint()
{
  if (free_memory) {
    for (auto i = _bits.begin(); i != _bits.end(); ++i) {
      delete (*i).second;
    }
  }

  return;
}

int
Global_Fingerprint::_increment_count(const DBKey& dbkey, Molecule& m,
                                     atom_number_t centre_atom, int acca)
{
  DBKey tmp(dbkey);

  bce_hash::iterator f = _bits.find(tmp);

  int radius = dbkey._radius;

  if (f == _bits.end()) {
    Count_Example* ce = new Count_Example(m, centre_atom, radius, acca);
    _bits[dbkey] = ce;

    return 1;
  }

  Count_Example* ce = (*f).second;

  if (radius == ce->radius() && acca == ce->atom_constant_centre_atom()) {
    ce->extra();
    return 1;
  }

  // Potential collision. If the new bit is for a smaller radius, prefer it

  cerr << "Potential collision on bit " << dbkey._bit << " stored " << ce->radius()
       << " new " << radius << '\n';

  if (radius >= ce->radius())  // nope, new radius is no shorter than what's already there
  {
    bit_collisions++;
    return 0;
  }

  Count_Example* newce = new Count_Example(m, centre_atom, radius, acca);

  (*f).second = newce;

  delete ce;

  return 1;
}

int
Global_Fingerprint::update(const Set_of_Bits<Bit_Produced>& sob, Molecule& m)
{
  for (auto i : sob) {
    const DBKey& b = i.first;

    const Bit_Produced* bp = i.second;

    _increment_count(b, m, bp->centre_atom(), bp->atom_constant_centre_atom());
  }

  return 1;
}

int
Global_Fingerprint::do_store(Db& fingerprint_db, const int max_shell_radius) const
{
  if (verbose) {
    cerr << "Storing " << _bits.size() << " bits to radius " << max_shell_radius << '\n';
  }

  int* bits_at_radius;
  if (verbose) {
    bits_at_radius = new_int(max_shell_radius + 1);
  } else {
    bits_at_radius = nullptr;
  }

  for (bce_hash::const_iterator i = _bits.begin(); i != _bits.end(); ++i) {
    const DBKey& b = (*i).first;

    //  cerr << "Storing bit " << b._b << " at radius " << static_cast<int>(b._radius) <<
    //  '\n';

    Count_Example* ce = (*i).second;

    if (!ce->store_bit(b, fingerprint_db)) {
      cerr << "Cannot store bit " << b._bit << '\n';
      return 0;
    }

    if (verbose) {
      bits_at_radius[b._radius]++;
    }
  }

  if (verbose) {
    for (int i = 0; i <= max_shell_radius; i++) {
      cerr << "Stored " << bits_at_radius[i] << " bits for radius " << i << '\n';
    }

    delete[] bits_at_radius;
  }

  return 1;
}

int
Global_Fingerprint::do_write(const char* fname) const
{
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "Global_Fingerprint::do_write:cannot open '" << fname << "'\n";
    return 0;
  }

  return do_write(output);
}

int
Global_Fingerprint::do_write(IWString_and_File_Descriptor& output) const
{
  for (auto i : _bits) {
    i.second->do_write(i.first, output);

    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();

  return 1;
}

static int
iwecfp(Molecule& m, Fingerprint_Characteristics& fc,
       Global_Fingerprint& global_fingerprint)
{
  Set_of_Bits<Bit_Produced> sob;

  if (!iwecfp_database::compute_fingerprints(m, fc, sob)) {  // ignore errors here
    return 1;
  }

#ifdef ECHO_FINGERPRINT
  cerr << "Fingerprint contains " << sob.size() << " bits\n";
  for (auto i : sob) {
    cerr << " bot " << i.first._bit << " rad " << static_cast<int>(i.first._radius)
         << " acca " << i.first._acca << '\n';
  }
#endif

  global_fingerprint.update(sob, m);

  return 1;
}

static int
iwecfp(data_source_and_type<Molecule>& input, Fingerprint_Characteristics& fc,
       Global_Fingerprint& global_fingerprint)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (report_progress()) {
      cerr << "Processed " << molecules_read << " molecules, global fingerprint contains "
           << global_fingerprint.size() << " items\n";
    }

    if (!preprocess_molecule(*m)) {
      cerr << "Skipping non organic or too large '" << m->name() << "'\n";
      continue;
    }

    (void)m->compute_aromaticity_if_needed();

    if (!iwecfp(*m, fc, global_fingerprint)) {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
iwecfp(const char* fname, FileType input_type, Fingerprint_Characteristics& fc,
       Global_Fingerprint& global_fingerprint)
{
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return iwecfp(input, fc, global_fingerprint);
}

/*static int
store_atom_type(Db & database,
                const Fingerprint_Characteristics & fc)
{
  IWString string_atype;
  if (! fc.string_atom_type(string_atype))    // cannot happen
    return 0;

  Dbt dkey((void *) (ATYPE_KEY), ::strlen(ATYPE_KEY));
  Dbt d;

  int rc = database.get(NULL, &dkey, &d, 0);

  if (DB_NOTFOUND == rc)   // the easy case, store a new value
  {
    d.set_data((void *) string_atype.rawchars());
    d.set_size(string_atype.length());

    int rc = database.put(NULL, &dkey, &d, 0);

    if (0 == rc)   // great
      return 1;

    database.err(rc, "Cannot store atom type key");
    return 0;
  }

// Must make sure stored value is the same as what we are using

  const_IWSubstring tmp(reinterpret_cast<const char *>(d.get_data()), d.get_size());

  if (tmp == string_atype)
    return 1;

  cerr << "INconsistent atom types in use, stored " << tmp << "', using '" << string_atype
<< "', cannot continue\n"; return 0;
}*/

static int
store_check_whatever(Db& database, const char* skey, const IWString& to_store)
{
  Dbt dkey((void*)(skey), ::strlen(skey));
  Dbt d;

  int rc = database.get(NULL, &dkey, &d, 0);

  if (DB_NOTFOUND == rc)  // the easy case, store a new value
  {
    d.set_data((void*)to_store.rawchars());
    d.set_size(to_store.length());

    int rc = database.put(NULL, &dkey, &d, 0);

    if (0 == rc) {  // great
      return 1;
    }

    cerr << "Cannot store '" << skey << "' " << db_strerror(rc) << '\n';
    return 0;
  }

  // Must make sure stored value is the same as what we are using

  const_IWSubstring tmp(reinterpret_cast<const char*>(d.get_data()), d.get_size());

  if (tmp == to_store) {
    return 1;
  }

  cerr << "INconsistent " << skey << " in use, stored " << tmp << "', using '" << to_store
       << "', cannot continue\n";
  return 0;
}

static int
store_atom_type(Db& database, const Fingerprint_Characteristics& fc)
{
  IWString string_atype;
  if (!fc.string_atom_type(string_atype)) {  // cannot happen
    return 0;
  }

  return store_check_whatever(database, ATYPE_KEY, string_atype);
}

static int
store_max_shell_radius(Db& database, const int r)
{
  IWString s;
  s << r;

  return store_check_whatever(database, RADIUS_KEY, s);
}

static void
display_dash_M_options(std::ostream& os)
{
  os << " -M noex      no example structures in DB, just counts and radius\n";
  os << " -M nfreem    do NOT bother freeing memory (may help speed)\n";

  exit(0);
}

static int
iwecfp(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vE:A:g:i:r:R:P:mld:ec:p:zy:oM:W:");

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
      cerr << "Cannot initialise chemical standardisation (-g)\n";
      usage(14);
    }
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  set_global_aromaticity_type(Daylight);
  set_input_aromatic_structures(1);

  if (iw_little_endian()) {
    need_to_swap_bytes = 1;
  }

  if (cl.option_present('M')) {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++)) {
      if ("noex" == m) {
        databases_contain_examples = 0;
        if (verbose) {
          cerr << "Will load just counts and radii into database\n";
        }
      } else if ("nfreem" == m) {
        free_memory = 0;
        if (verbose) {
          cerr << "Will NOT de-allocate the global_fingerprint structures\n";
        }
      } else if ("help" == m) {
        display_dash_M_options(cerr);
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_M_options(cerr);
      }
    }
  }

  Fingerprint_Characteristics fc;

  if (!fc.build(cl, 0, verbose)) {
    cerr << "Cannot initialise fingerprint characteristics\n";
    usage(1);
  }

  if (cl.option_present('l')) {
    only_set_bits_for_max_radius_shell = 1;
    if (verbose) {
      cerr << "Only bits for the largest radius will be set\n";
    }
  }

  if (cl.option_present('e')) {
    update_global_fingerprint_presence_only = 0;
    if (verbose) {
      cerr << "Will record total counts when storing\n";
    }
  }

  if (cl.option_present('p')) {
    if (!report_progress.initialise(cl, 'p', verbose)) {
      cerr << "Cannot initialise progress reporting (-p)\n";
      return 1;
    }
  }

  if (cl.option_present('o')) {
    store_flag = 0;

    if (verbose) {
      cerr << "Will overwrite existing data\n";
    }
  }

  if (cl.option_present('z')) {
    remove_leading_zeros = 1;

    if (verbose) {
      cerr << "Leading zeros stripped from identifiers\n";
    }
  }

  if (cl.option_present('c')) {
    only_start_at_chiral_atoms = 1;

    if (verbose) {
      cerr << "Building a chirality database, shells only start at chiral atoms\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == ::strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot determine input type(s)\n";
    return 7;
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify database via the -d option\n";
    usage(3);
  } else if (cl.option_count('d') > 1) {
    cerr << "Only one database possible when writing\n";
    usage(3);
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  DBTYPE dbtype = DB_BTREE;

  if (cl.option_present('y')) {
    const_IWSubstring t = cl.string_value('y');

    if ("btree" == t) {
      ;
    } else if ("hash" == t) {
      dbtype = DB_HASH;
    } else if ("recno" == t) {
      dbtype = DB_RECNO;
    } else if ("queue" == t) {
      dbtype = DB_QUEUE;
    } else {
      cerr << "Unrecognised access method (-y) '" << t << "'\n";
    }
  }

  if (cl.option_present('d')) {
    const char* dbname = cl.option_value('d');
    int oflags = DB_CREATE;
    int mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;
    int rc = database.open(NULL, dbname, NULL, dbtype, oflags, mode);
    if (0 != rc) {
      cerr << "Cannot open database '" << dbname << "'\n";
      database.err(rc, "");
      exit(23);
    }

    if (verbose) {
      cerr << "Storing into '" << dbname << "'\n";
    }
  }

  if (!store_atom_type(database, fc)) {
    cerr << "Cannot store atom type, or inconsistent atom types\n";
    return 2;
  }

  if (!store_max_shell_radius(database, fc.max_shell_radius())) {
    cerr << "Cannot store max radius, or inconsistent max radius\n";
    return 2;
  }

  Global_Fingerprint global_fingerprint;

  int rc = 0;

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!iwecfp(cl[i], input_type, fc, global_fingerprint)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";

    if (nbits_acc.n() > 0) {
      cerr << "Fingerprints had between " << nbits_acc.minval() << " and "
           << nbits_acc.maxval() << " ave "
           << static_cast<float>(nbits_acc.average_if_available_minval_if_not())
           << " bits set\n";
    }
  }

  global_fingerprint.do_store(database, fc.max_shell_radius());

  database.close(0);

  if (bit_collisions) {
    cerr << bit_collisions << " bit collisions\n";
  }

  if (intra_molecular_bit_collisions) {
    cerr << intra_molecular_bit_collisions << " intra_molecular_bit_collisions\n";
  }

  if (duplicate_keys_not_stored) {
    cerr << duplicate_keys_not_stored << " duplicate_keys_not_stored\n";
  }

  if (cl.option_present('W')) {
    const char* w = cl.option_value('W');
    if (verbose) {
      cerr << "Writing db to '" << w << "'\n";
    }

    global_fingerprint.do_write(w);
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = iwecfp(argc, argv);
  return rc;
}
