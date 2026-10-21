/*
  Produces unique smiles for ring systems within a molecule
  Can store these ring systems in a Berkeley database, or can lookup
  ring systems in a Berkeley database
*/

#include <math.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/string_with_commas.h"

#include "db_cxx.h"

#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools_Bdb/smi2rings.pb.h"
#else
#include "smi2rings.pb.h"
#endif

#include "Molecule_Tools_Bdb/ring_database_builder.h"
#include "Molecule_Tools_Bdb/ring_database_lookup.h"
#include "Molecule_Tools_Bdb/ring_extraction.h"
#include "Molecule_Tools_Bdb/ring_preprocessing.h"

using std::cerr;

// Shared chemistry configuration is owned by reusable library objects. Only
// command-line input/output settings and job statistics remain in this wrapper.
static smi2rings_bdb::RingExtractor ring_extractor;
static smi2rings_bdb::RingPreprocessor preprocessor;
static smi2rings_bdb::RingDatabaseBuildOptions build_options;

struct LookupOptions {
  int _min_urs_needed = 0;
  int _report_all_ring_uniqueness_measurements_looked_up = 1;
  IWString _urs_string = " URS:";
  int _append_largest_ring_size_if_no_match = 0;
  int _max_exemplar_count = 20;
};

static LookupOptions lookup_options;

static int verbose = 0;
static int lookup_in_database = 0;
static int store_in_database = 0;
static resizable_array_p<Db> databases;
static std::unique_ptr<smi2rings_bdb::RingDatabaseBuilder> database_builder;

static int write_rings = 1;
static IWString smiles_tag("$SMI<");
static int write_parent_structure = 0;
static int use_parent_name_for_rings = 0;
static int write_non_ring_atoms = 0;
static int smiles_and_textproto = 0;
static int full_proto_output = 0;
static IWString_and_File_Descriptor stream_for_rejected_molecules;
static IWString_and_File_Descriptor stream_for_missing_ring_systems;
static Report_Progress report_progress;

static uint64_t molecules_read = 0;
static int discarded_for_urs_too_low = 0;
static int molecules_with_unique_ring_systems_found = 0;
static int ring_systems_not_found_in_database = 0;
static extending_resizable_array<int> ring_system_size_array;
static Accumulator_Int<uint64_t> ring_system_ring_stats;
static Accumulator_Int<uint64_t> ring_system_atom_stats;
static Accumulator_Int<uint64_t> acc_exemplar;
static extending_resizable_array<uint64_t> exemplar_count;

// CLI output ranking and annotation, separate from library lookup results.
class Ring_System_Uniqueness {
 private:
  int _uid;
  int _uniqueness_measure;

  //  If we are looking up in a database with names, we can return
  //  the name of the molecule that first saw that ring in the database

  IWString _name_from_database;

  // If we are doing lookups across multiple databases, we keep track
  // of the number of databases in which a ring is found.
  int _ndb_found;

 public:
  Ring_System_Uniqueness();

  void
  set_uid(int s) {
    _uid = s;
  }

  int
  uid() const {
    return _uid;
  }

  void
  set_uniqueness_measure(int s) {
    _uniqueness_measure = s;
  }

  int
  uniqueness_measure() const {
    return _uniqueness_measure;
  }

  void
  set_name_from_database(const const_IWSubstring& s) {
    _name_from_database = s;
  }

  const IWString&
  name_from_database() const {
    return _name_from_database;
  }

  void
  set_ndb_found(int s) {
    _ndb_found = s;
  }

  int
  ndb_found() const {
    return _ndb_found;
  }
};

Ring_System_Uniqueness::Ring_System_Uniqueness() {
  _uid = -1;
  _uniqueness_measure = std::numeric_limits<int>::max();
  _ndb_found = 0;
}

class Ring_System_Uniqueness_Comparitor {
 private:
 public:
  int operator()(const Ring_System_Uniqueness&, const Ring_System_Uniqueness&) const;
};

int
Ring_System_Uniqueness_Comparitor::operator()(const Ring_System_Uniqueness& rs1,
                                              const Ring_System_Uniqueness& rs2) const {
  int n1 = rs1.uniqueness_measure();
  int n2 = rs2.uniqueness_measure();

  if (n1 < n2) {
    return -1;
  }

  if (n1 > n2) {
    return 1;
  }

  return 0;
}

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
  cerr << R"(Identifies rings and ring systems. Stores or retrieves by unique smiles of the ring.
Data is stored in a BerkeleyDb database with key unique smiles and value a Smi2Rings.Ring textproto
describing the ring system, including parent structure and number of instances.
Usage:
smi2rings_bdb -d dbname -d STORE|LOOKUP <options> file.smi
  -j ...         options for adding atoms to rings, enter '-j help' for info
  -N add         add non-ring connection to positive aromatic nitrogen
  -N neut        convert positive aromatic nitrogen to neutral form
  -y             put isotopes on ring atoms to indicate ring system membership
  -z <size>      discard rings or ring systems with a ring of size > <size>
  -R <max>       max ring system size (number rings) to consider
  -d <dbname>    specify database for lookup or storage
  -d STORE       store ring systems found in the database
  -n             suppress output if writing to a database
  -d LOOKUP      lookup each ring system in the database
  -U <string>    string to append with the lowest ring count (default 'URS')
  -p w           write the parent structure - all atoms
  -p n           use the name of the parent for the name of the rings
  -s <min>       write non-ring fragments with at least <min> atoms
  -f <cutoff>    function as filter (lookup), discard molecules below <cutoff> URS
                 -f file=<fname> to write rejected molecules to <fname>
  -S <query>     one or more queries for functional groups. Typically single atom ring substituents
  -a             only report rarest in lookups
  -e             discard rings where the ring is the entire molecule
  -c             remove all chiral centres
  -Y ...         obscure options, enter '-Y help' for info
  -E ...         standard element options, enter '-E help' for info
  -A ...         standard aromaticity options
  -g ...         chemical standardisation options
  -v             verbose output
Typical usage
  smi2rings_bdb -d rings.bdb -d STORE -j ring -j iso -j env -j double -N add -n -c collection.smi
# To filter to molecules where the rarest ring has at least 10 instances
  smi2rings_bdb -d rings.bdb -d LOOKUP -n -f 10 -v -c -j double -N add -j iso query.smi
)";
  // clang-format on

  exit(rc);
}

struct PerMoleculeData : public smi2rings_bdb::RingExtractionData {
  std::unique_ptr<Ring_System_Uniqueness[]> _rsu;

  PerMoleculeData(Molecule& m, bool lookup)
      : RingExtractionData(m, ring_extractor.options()),
        _rsu(lookup ? std::make_unique<Ring_System_Uniqueness[]>(m.nrings()) : nullptr) {
  }
};

static int
all_fragments_too_small(Molecule& m, int min_atoms_per_fragment) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> fragment_membership = m.fragment_membership();

  int* to_remove = new_int(matoms);
  std::unique_ptr<int[]> free_to_remove(to_remove);

  int atoms_removed = 0;

  for (int i = 0; i < matoms; ++i) {
    const int frag = fragment_membership[i];
    if (m.atoms_in_fragment(frag) < min_atoms_per_fragment) {
      to_remove[i] = 1;
      ++atoms_removed;
    }
  }

  if (atoms_removed == 0) {
    return 0;
  }

  if (atoms_removed == matoms) {
    return 1;
  }

  m.remove_atoms(to_remove);

  return 0;
}

static int
do_write_non_ring_atoms(Molecule& m, int* in_system,
                        IWString_and_File_Descriptor& output) {
  int matoms = m.natoms();
  int atoms_in_subset = 0;
  for (int i = 0; i < matoms; i++) {
    if (in_system[i]) {  // was a ring atom, don't want those
      in_system[i] = 0;
    } else {
      in_system[i] = 1;
      atoms_in_subset++;
    }
  }

  if (atoms_in_subset < write_non_ring_atoms) {
    return 1;
  }

  Molecule subset;
  m.create_subset(subset, in_system, 1);

  if (subset.number_fragments() > 1) {
    if (all_fragments_too_small(subset, write_non_ring_atoms)) {
      return 1;
    }
  }

  if (ring_extractor.options()._remove_invalid_chiral_centres) {
    lillymol::do_remove_invalid_chiral_centres(subset);
  }

  output << subset.smiles() << ' ' << m.name() << " non ring\n";

  return output.good();
}

static void
WriteSmiles(Molecule& m, IWString_and_File_Descriptor& output) {
  if (ring_extractor.options()._write_3d_smiles_if_needed &&
      m.highest_coordinate_dimensionality() == 3) {
    lillymol::set_include_coordinates_with_smiles(1);
  }

  output << m.smiles();

  lillymol::set_include_coordinates_with_smiles(0);
}

static int
write_ring(Molecule& m, const IWString& unique_smiles,
           IWString_and_File_Descriptor& output) {
  output << unique_smiles;
  if (use_parent_name_for_rings) {
    output << ' ' << m.name() << 'R';
  }
  output << '\n';
  return output.good();
}

static Ring_System_Uniqueness_Comparitor ring_system_uniqueness_comparitor;

// Molecule `m` has been disarded for having a lowest ring count value
// too low. Increment global counter and maybe write.
static int
HandleDiscardedForMinUrs(Molecule& m, int lowest_ring_count) {
  discarded_for_urs_too_low++;
  if (stream_for_rejected_molecules.active()) {
    stream_for_rejected_molecules << m.smiles() << ' ' << m.name()
                                  << lookup_options._urs_string << lowest_ring_count
                                  << '\n';
    stream_for_rejected_molecules.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
MaybeWriteNotFoundRings(Molecule& m, const IWString& unique_smiles) {
  if (!stream_for_missing_ring_systems.is_open()) {
    return 0;
  }

  stream_for_missing_ring_systems << m.smiles() << ' ' << m.name() << '\n';
  stream_for_missing_ring_systems << unique_smiles << '\n';
  stream_for_missing_ring_systems.write_if_buffer_holds_more_than(32768);

  return 1;
}

// Return one for a hit, zero for a miss and minus one for a database error.
// Missing-ring output and CLI counters remain outside the database library.
static int
do_database_lookup(Molecule& m, PerMoleculeData& mdata, int uid,
                   resizable_array_p<Db>& databases, const IWString& unique_smiles) {
  smi2rings_bdb::RingLookupResult result;
  const smi2rings_bdb::LookupStatus status = smi2rings_bdb::LookupRingInDatabases(
      std::span<Db* const>(databases.rawdata(), databases.size()), unique_smiles, result);
  if (status == smi2rings_bdb::LookupStatus::kError) {
    return -1;
  }
  mdata._rsu[uid].set_uniqueness_measure(result._count);
  mdata._rsu[uid].set_name_from_database(IWString(result._exemplar));
  mdata._rsu[uid].set_ndb_found(result._databases_matched);
  if (status == smi2rings_bdb::LookupStatus::kFound) {
    return 1;
  }
  if (databases.size() == 1) {
    ++ring_systems_not_found_in_database;
  }
  MaybeWriteNotFoundRings(m, unique_smiles);
  return 0;
}

// These operations are ordered in a specific way, only for doing lookups.
// The most specific lookup is done first.

static int
DoLookups(Molecule& m, int ring_number, PerMoleculeData& mdata,
          resizable_array_p<Db>& databases, IWString_and_File_Descriptor& output) {
  for (smi2rings_bdb::RingKeyType type : ring_extractor.KeyTypes()) {
    IWString key;
    if (!ring_extractor.CreateKey(m, mdata, ring_number, type, key)) {
      return -1;
    }
    if (write_rings && type != smi2rings_bdb::RingKeyType::kBare) {
      write_ring(m, key, output);
    }
    if (const int rc = do_database_lookup(m, mdata, ring_number, databases, key);
        rc != 0) {
      return rc;
    }
    if (type != smi2rings_bdb::RingKeyType::kBare) {
      m.transform_to_non_isotopic_form();
    }
  }
  return 1;
}

static int
DoStore(Molecule& m, int ring_number, PerMoleculeData& mdata,
        IWString_and_File_Descriptor& output) {
  for (smi2rings_bdb::RingKeyType type : ring_extractor.KeyTypes(true)) {
    IWString key;
    if (!ring_extractor.CreateKey(m, mdata, ring_number, type, key, true)) {
      return -1;
    }
    if (write_rings && type != smi2rings_bdb::RingKeyType::kBare) {
      write_ring(m, key, output);
    }
    if (!mdata.AlreadySeen(key) && !database_builder->Store(key, m.name())) {
      return -1;
    }
    if (type != smi2rings_bdb::RingKeyType::kBare) {
      m.transform_to_non_isotopic_form();
    }
  }
  return 1;
}

static int
WriteProto(Molecule& m, const PerMoleculeData& mdata, int number_ring_systems,
           IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  Smi2Rings::Results results;
  int nprocess = lookup_options._report_all_ring_uniqueness_measurements_looked_up
                     ? number_ring_systems
                     : 1;
  for (int i = 0; i < nprocess; ++i) {
    Smi2Rings::Result* r = results.add_result();
    const auto& u = mdata._rsu[i];
    r->set_urs(u.uniqueness_measure());
    if (u.uniqueness_measure() > 0) {
      const IWString& s = u.name_from_database();
      r->set_ex(s.data(), s.length());
    }
  }

  if (full_proto_output) {
    const IWString& s = m.smiles();
    results.set_smiles(s.data(), s.length());
    const IWString& n = m.name();
    results.set_id(n.data(), n.length());
  } else {
    output << m.smiles() << ' ' << m.name() << ' ';
  }

  std::string buffer;
  if (!printer.PrintToString(results, &buffer)) {
    cerr << "PrintToString failed\n";
    return 0;
  }

  output << buffer << '\n';

  return 1;
}

// Update various file scope statistics

static void
UpdateGlobalStatistics(const PerMoleculeData& mdata, int ring_system_number) {
  if (verbose == 0) {
    return;
  }

  int rings_in_system = mdata._ring_system_size[ring_system_number];

  ring_system_size_array[rings_in_system]++;
  ring_system_ring_stats.extra(rings_in_system);

  int atoms_in_system =
      std::count(mdata._in_system, mdata._in_system + mdata.natoms(), ring_system_number);

  ring_system_atom_stats.extra(atoms_in_system);
}

static int
LargestRing(Molecule& m) {
  const Ring* r = m.ringi(m.nrings() - 1);

  return r->number_elements();
}

template <bool LookupMode>
static int
smi2rings(Molecule& m, PerMoleculeData& mdata, IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  const int n = ring_extractor.IdentifySystems(m, mdata);

  for (int i = 0; i < n; i++) {
    if constexpr (LookupMode) {
      mdata._rsu[i].set_uid(i);
    }

    UpdateGlobalStatistics(mdata, i);

    if (!ring_extractor.PrepareSystem(m, mdata, i)) {
      continue;
    }

    if constexpr (LookupMode) {
      if (DoLookups(m, i, mdata, databases, output) < 0) {
        return 0;
      }
    } else {
      if (DoStore(m, i, mdata, output) < 0) {
        return 0;
      }
    }
  }

  if (write_non_ring_atoms) {
    do_write_non_ring_atoms(m, mdata._in_system, output);
  }

  if constexpr (!LookupMode) {
    return output.good();
  }

  // If we were doing a lookup, check the results

  if (ring_extractor.options()._label_attachment_points ||
      ring_extractor.options()._label_attachment_points_with_environment ||
      ring_extractor.options()._always_process_unsubstituted_ring) {
    m.transform_to_non_isotopic_form();
  }

  if (n > 1) {
    iwqsort(mdata._rsu.get(), n, ring_system_uniqueness_comparitor);
  }

#ifdef ECHO_SORTED_LIST
  if (i > 1) {
    for (int i = 0; i < n; i++) {
      cerr << i << ' ' << rsu[i].uniqueness_measure() << '\n';
    }
  }
#endif

  int lowest_ring_count = mdata._rsu[0].uniqueness_measure();

  // none of the ring systems were looked up
  if (std::numeric_limits<int>::max() == lowest_ring_count) {
    return output.good();
  }

  if (verbose > 1) {
    cerr << m.name() << " of " << n << " ring systems, least common ring system count "
         << lowest_ring_count << '\n';
  }

  acc_exemplar.extra(lowest_ring_count);
  if (lowest_ring_count <= lookup_options._max_exemplar_count) {
    ++exemplar_count[lowest_ring_count];
  }

  if (0 == lowest_ring_count) {
    molecules_with_unique_ring_systems_found++;
  }

  if (lowest_ring_count < lookup_options._min_urs_needed) {
    return HandleDiscardedForMinUrs(m, lowest_ring_count);
  }

  if (smiles_and_textproto || full_proto_output) {
    return WriteProto(m, mdata, n, output);
  }

  IWString tmp(m.name());

  tmp << lookup_options._urs_string << mdata._rsu[0].uniqueness_measure();
  if (0 == mdata._rsu[0].uniqueness_measure()) {
    tmp << " .";
    if (lookup_options._append_largest_ring_size_if_no_match) {
      tmp << ' ' << LargestRing(m);
    }
  } else {
    tmp << ' ' << mdata._rsu[0].name_from_database();
  }

  if (lookup_options._report_all_ring_uniqueness_measurements_looked_up) {
    for (int i = 1; i < n; i++) {
      tmp << ' ' << mdata._rsu[i].uniqueness_measure();

      if (i != (n - 1)) {
        tmp << ',';
      }
    }
  }

#ifdef DEBUG_WRITE_RING
  cerr << ring_extractor.options()._label_ring_atoms << " label_ring_atoms, n = " << n
       << '\n';
  for (int i = 0; i < matoms; ++i) {
    cerr << i << " in_system " << mdata._in_system[i] << ' '
         << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

  mdata.MaybeRevertInitialIsotopes(m);

  if (ring_extractor.options()._label_ring_atoms) {
    int uid = mdata._rsu[0].uid();
    for (int i = 0; i < matoms; ++i) {
      if (mdata._in_system[i] == uid) {
        m.set_isotope(i, 1);
      }
    }
    IWString name(m.name());
    name << lookup_options._urs_string << mdata._rsu[0].uniqueness_measure();
    m.set_name(name);
  }

  WriteSmiles(m, output);
  output << ' ' << tmp << '\n';

  return output.good();
}

template <bool LookupMode>
static int
smi2rings(Molecule& m, IWString_and_File_Descriptor& output) {
  if (!report_progress()) {
  } else if constexpr (LookupMode) {
    cerr << "Read " << iwstring::with_commas(molecules_read) << ", "
         << iwstring::with_commas(exemplar_count[0]) << " unique ring systems found\n";
  } else {
    cerr << "Read " << iwstring::with_commas(molecules_read) << " molecules, found "
         << iwstring::with_commas(database_builder->new_ring_systems())
         << " new ring systems\n";
  }

  if (write_parent_structure || write_non_ring_atoms) {
    WriteSmiles(m, output);
    output << ' ' << m.name() << '\n';
  }

  if (m.nrings() == 0) {  // that's OK, we just have nothing to do
    return 1;
  }

  PerMoleculeData mdata(m, LookupMode);
  if (!ring_extractor.Prepare(m, mdata)) {
    return 0;
  }

  return smi2rings<LookupMode>(m, mdata, output);
}

template <bool LookupMode>
static int
smi2rings(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule()) && output.good()) {
    std::unique_ptr<Molecule> free_m(m);

    preprocessor.Process(*m);

    molecules_read++;

    if (!smi2rings<LookupMode>(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(16384);
  }

  return output.good();
}

template <bool LookupMode>
static int
smi2rings(const char* fname, FileType input_type, IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "smi2rings: cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return smi2rings<LookupMode>(input, output);
}

static void
display_dash_j_options(std::ostream& os) {
  // clang-format off
  os << R"(
 -j double      add doubly bonded neighbours to rings - always do this.
 -j all         add all singly bonded neighbours to rings.
 -j fg          add only singly bonded functional groups to rings.
 -j unsat       add all atoms bonded to an unsaturated neighbour.
 -j spiro       include spiro fusions in ring systems.
 -j safg        include all single atom functional groups, CH3, NH2, OH, F etc.
 -j iso         isotopically label the attachment points.
 -j isocon      isotope number is number of outside-system attachments.
 -j env         isotope number reflects environment of attachment.
 -j isoatt      isotope number is the highest isotope of attached atoms.
                use '-j env=UST:achry' to specify an atom type for the env.
 -j ring        always process the bare ring - in addition to isotopic forms.
When building a database (-d STORE) use all of '-j ring', '-j iso' and '-j env'.
When doing lookups use one or more of them.
)";
  // clang-format on

  exit(0);
}

static void
DisplayDashYOptions(std::ostream& output) {
  // clang-format off
  cerr << R"( -Y smi3d        write 3D smiles when 3D information present.
 -Y smiproto     write smiles, id followed by textproto output.
 -Y proto        write as Smi2Rings::Results textproto.
 -Y writemissing=<fname> during lookup write any missing rings to <fname>.smi.
 -Y hash         when storing, hash all results and write database at end (recommended).
 -Y rpt=<n>      when storing, report progress every <n> molecules processed.
 -Y applarge     if a molecule has a ring not in the database, also write the size of the largest ring.
                 Quite possible that the ring is larger than what is in the database.
                 Adds an extra token to the output 'smiles id URS: 0 12' for example.
 -Y scaffold     Store and lookup by scaffold instead of ring systems.
)";
  // clang-format on

  usage(1);
}

template <bool LookupMode>
static int
ProcessInputs(const Command_Line& cl, FileType input_type,
              IWString_and_File_Descriptor& output) {
  for (int i = 0; i < cl.number_elements(); ++i) {
    if (!smi2rings<LookupMode>(cl[i], input_type, output)) {
      return i + 1;
    }
  }
  return 0;
}

static int
smi2rings(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:d:j:D:g:nN:s:R:f:z:S:p:U:uycCaeY:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  if (cl.option_present('U')) {
    cl.value('U', lookup_options._urs_string);

    if (verbose) {
      cerr << "URS string set to '" << lookup_options._urs_string << "'\n";
    }

    if (lookup_options._urs_string.empty()) {
      cerr << "Warning, zero length Unique Ring System string\n";
    }
  }

  if (cl.option_present('g')) {
    if (!preprocessor.standardisation().construct_from_command_line(cl, verbose > 1,
                                                                    'g')) {
      usage(6);
    }
  }

  if (cl.option_present('j')) {
    int i = 0;
    const_IWSubstring j;
    while (cl.value('j', j, i++)) {
      if ("all" == j) {
        ring_extractor.options().set_add_single_bonded_neighbours_to_rings(2);
        if (verbose) {
          cerr << "Will include singly bonded atoms outside the ring system\n";
        }
      } else if ("iso" == j) {
        ring_extractor.options().set_label_attachment_points(1);
        if (verbose) {
          cerr << "Isotopic labels applied at attachment points\n";
        }
      } else if ("isocon" == j) {
        ring_extractor.options().set_label_attachment_points(2);
        if (verbose) {
          cerr << "The isotopic labels applied will be the number of connections\n";
        }
      } else if ("spiro" == j) {
        ring_extractor.options().set_bridge_across_spiro_groups(1);
        if (verbose) {
          cerr << "Will bridge across spiro fusions for ring systems\n";
        }
      } else if ("fg" == j) {
        ring_extractor.options().set_add_single_bonded_neighbours_to_rings(1);
        if (verbose) {
          cerr << "Will add singly bonded, identified functional groups to the ring "
                  "system\n";
        }
      } else if ("double" == j) {
        ring_extractor.options().set_add_double_bonded_neighbors_to_rings(1);
        if (verbose) {
          cerr << "Multiple bonds outside the ring system will be included\n";
        }
      } else if ("unsat" == j) {
        ring_extractor.options().set_continue_chain_neighbours_through_multiple_bonds(1);
        if (verbose) {
          cerr << "Will continue chain neighbours through multiple bonds\n";
        }
      } else if ("glaxo" == j) {
        ring_extractor.options().set_label_attachment_points(1);
        ring_extractor.options().set_add_double_bonded_neighbors_to_rings(1);
      } else if ("safg" == j) {
        ring_extractor.options().set_add_all_single_atom_ring_substituents(1);
        if (verbose) {
          cerr << "Will add all single atom functional groups\n";
        }
      } else if ("env" == j) {
        ring_extractor.options().set_label_attachment_points_with_environment(1);
        if (verbose) {
          cerr << "isotopic labels will reflect environment\n";
        }
      } else if (j == "isoatt") {
        ring_extractor.options().set_label_attachment_points_with_existing_isotope(1);
        if (verbose) {
          cerr << "Isotopic labels are the isotope of the attached atom\n";
        }
      } else if (j.starts_with("env=")) {
        j.remove_leading_chars(4);

        if (!ring_extractor.atom_typing().build(j)) {
          cerr << "Cannot determing atom typing for environment '" << j << "'\n";
          return 3;
        }
        ring_extractor.options().set_label_attachment_points_with_environment(1);
        if (verbose) {
          cerr << "Isotopic labels will reflect environment with type '" << j << "'\n";
        }
      } else if (j == "ring") {
        ring_extractor.options().set_always_process_unsubstituted_ring(1);
        if (verbose) {
          cerr << "Will always examine the non-isotopically labelled ring\n";
        }
      } else if ("help" == j) {
        display_dash_j_options(cerr);
      } else {
        cerr << "Unrecognised -j qualifier '" << j << "'\n";
        display_dash_j_options(cerr);
      }
    }
  }

  if (ring_extractor.options()._label_attachment_points_with_environment) {
  } else if (ring_extractor.options()._label_attachment_points) {
  } else if (ring_extractor.options()._always_process_unsubstituted_ring) {
  } else {
    cerr << "Default action, store labelled attachment points\n";
    ring_extractor.options().set_label_attachment_points(1);
  }

  if (cl.option_present('S')) {
    if (!process_queries(cl, ring_extractor.functional_group_queries(), verbose, 'S')) {
      cerr << "Cannot parse functional group queries\n";
      usage(4);
    }

    ring_extractor.options().set_add_single_bonded_neighbours_to_rings(1);

    if (verbose) {
      cerr << "Defined " << ring_extractor.functional_group_queries().number_elements()
           << " queries for singly bonded functional groups\n";
    }

    for (int i = 0; i < ring_extractor.functional_group_queries().number_elements();
         i++) {
      ring_extractor.functional_group_queries()[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('a')) {
    lookup_options._report_all_ring_uniqueness_measurements_looked_up = 0;
    if (verbose) {
      cerr << "Will only report the rarest ring system during lookups\n";
    }
  }

  if (cl.option_present('z')) {
    int max_ring_size_to_consider;
    if (!cl.value('z', max_ring_size_to_consider) || max_ring_size_to_consider < 3) {
      cerr << "The max ring size to consider option (-z) must be followed by a whole "
              "number > 2\n";
      usage(41);
    }

    ring_extractor.options().set_max_ring_size_to_consider(max_ring_size_to_consider);

    if (verbose) {
      cerr << "Will not consider rings larger than "
           << ring_extractor.options()._max_ring_size_to_consider << '\n';
    }
  }

  if (cl.option_present('e')) {
    ring_extractor.options().set_discard_molecules_where_ring_is_whole_molecule(1);

    if (verbose) {
      cerr << "Discard molecules where the ring is the whole molecule\n";
    }
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify database(s) to process via the -d option\n";
    usage(1);
  }

  if (cl.option_present('d')) {
    int oflags = 0;
    DBTYPE dbtype = DB_UNKNOWN;

    resizable_array<IWString> dbnames;

    IWString d;
    for (int i = 0; cl.value('d', d, i); ++i) {
      if ("STORE" == d) {
        oflags = DB_CREATE;
        dbtype = DB_BTREE;
        store_in_database = 1;
      } else if (d.starts_with("LOOK")) {
        oflags = DB_RDONLY;
        dbtype = DB_UNKNOWN;
        lookup_in_database = 1;
      } else if (d == "NAME") {
        // obsolete option now ignored.
      } else {
        dbnames << d;
      }
    }

    if (lookup_in_database && store_in_database) {
      cerr << "Cannot specify both store and lookup operations\n";
      usage(4);
    }

    if (lookup_in_database && !cl.option_present('n')) {
      cerr << "When doing database lookups, must add the -n option\n";
      usage(17);
    }

    if (0 == store_in_database && 0 == lookup_in_database) {
      cerr << "Must specify either '-d STORE' or '-d LOOKUP'\n";
      usage(17);
    }

    databases.reserve(dbnames.number_elements());
    for (IWString& dbname : dbnames) {
      std::unique_ptr<Db> db(new Db(NULL, DB_CXX_NO_EXCEPTIONS));
      if (int rc =
              db->open(NULL, dbname.null_terminated_chars(), 0, dbtype, oflags, 0644);
          rc != 0) {
        cerr << "Cannot open database '" << dbname << "'\n";
        db->err(rc, "");
        return 1;
      }
      if (verbose) {
        cerr << "Opened database '" << dbname << "'\n";
      }
      databases << db.release();
    }
  }

  if (cl.option_present('y')) {
    if (!lookup_in_database) {
      cerr << "The label atoms by ring system option (-y) only makes sense when doing "
              "lookups\n";
      usage(5);
    }

    ring_extractor.options().set_label_ring_atoms(1);

    if (verbose) {
      cerr << "Atoms labelled by ring system membership\n";
    }
  }

  if (cl.option_present('u')) {
    cerr << "The -u option is now ignored\n";
  }

  if (cl.option_present('R')) {
    int max_ring_system_size_to_consider;
    if (!cl.value('R', max_ring_system_size_to_consider) ||
        max_ring_system_size_to_consider < 1) {
      cerr << "The -R option must be followed by a whole positive number\n";
      usage(39);
    }

    ring_extractor.options().set_max_ring_system_size_to_consider(
        max_ring_system_size_to_consider);

    if (verbose) {
      cerr << "Ring systems with more than "
           << ring_extractor.options()._max_ring_system_size_to_consider
           << " rings will be ignored\n";
    }
  }

  if (cl.option_present('n')) {
    write_rings = 0;
    if (verbose) {
      cerr << "Default output to stdout suppressed\n";
    }
  }

  if (cl.option_present('N')) {
    int i = 0;
    const_IWSubstring n;
    while (cl.value('N', n, i++)) {
      if ("add" == n) {
        ring_extractor.options().set_add_connected_atom_to_positive_aromatic_nitrogen(1);
        if (verbose) {
          cerr << "Will add connected atoms to positive aromatic nitrogens\n";
        }
      } else if ("neut" == n) {
        ring_extractor.options().set_convert_positive_aromatic_nitrogen_to_neutral(1);
        if (verbose) {
          cerr << "Will convert positive aromatic nitrogen to neutral form\n";
        }
      } else {
        cerr << "Unrecognised -N qualifier '" << n << "'\n";
        usage(5);
      }
    }
  }

  if (cl.option_present('p')) {
    int i = 0;
    const_IWSubstring p;

    while (cl.value('p', p, i++)) {
      if ('w' == p) {
        write_parent_structure = 1;

        if (verbose) {
          cerr << "Will write the parent structure\n";
        }
      } else if ('n' == p) {
        use_parent_name_for_rings = 1;

        if (verbose) {
          cerr << "Each ring inherits the name of its parent structure\n";
        }
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(5);
      }
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', write_non_ring_atoms) || write_non_ring_atoms < 1) {
      cerr << "The minimum number of atoms in non-ring fragments to write (-s) must be a "
              "whole positive number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will write non ring fragments with " << write_non_ring_atoms
           << " or more atoms\n";
    }
  }

  if (cl.option_present('f')) {
    lookup_in_database = 1;

    IWString fname;
    const_IWSubstring f;
    for (int i = 0; cl.value('f', f, i); ++i) {
      if (f.starts_with("file=") || f.starts_with("FILE=")) {
        f.remove_leading_chars(5);
        fname = f;
      } else if (!f.numeric_value(lookup_options._min_urs_needed) ||
                 lookup_options._min_urs_needed < 1) {
        cerr << "The filter option requires a whole non-negative cutoff (-f)\n";
        usage(2);
      }
    }

    if (verbose) {
      cerr << "Will discard molecules with URS values below "
           << lookup_options._min_urs_needed << "\n";
    }

    if (!fname.empty()) {
      fname.EnsureEndsWith(".smi");
      if (!stream_for_rejected_molecules.open(fname)) {
        cerr << "Cannot open stream for discarded molecules '" << fname << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Discarded molecules written to '" << fname << "'\n";
      }
    }
  }

  if (cl.option_present('c')) {
    preprocessor.options().set_remove_chirality(1);
    if (verbose) {
      cerr << "Chirality will be discarded\n";
    }
  } else if (cl.option_present('C')) {
    ring_extractor.options().set_remove_invalid_chiral_centres(1);
    if (verbose) {
      cerr << "Will remove any invalid chiral centres in rings produced\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "smi3d") {
        ring_extractor.options().set_write_3d_smiles_if_needed(1);
        if (verbose) {
          cerr << "Will write 3D smiles if coordinates present\n";
        }
      } else if (y == "proto") {
        full_proto_output = 1;
        if (verbose) {
          cerr << "Will write Smi2Rings::Results textproto output\n";
        }
      } else if (y == "smiproto") {
        smiles_and_textproto = 1;
        if (verbose) {
          cerr << "Will write smiles + Smi2Rings::Results textproto\n";
        }
      } else if (y.starts_with("writemissing=")) {
        y.remove_leading_chars(13);
        IWString fname(y);
        if (!fname.ends_with(".smi")) {
          fname << ".smi";
        }
        if (!stream_for_missing_ring_systems.open(fname)) {
          cerr << "Cannot open stream for missing rings '" << fname << "'\n";
          return 1;
        }
        if (verbose) {
          cerr << "Missing rings written to '" << fname << "'\n";
        }
      } else if (y == "hash") {
        build_options.set_cache(1);
        if (verbose) {
          cerr << "Will cache rings and flush at the end\n";
        }
      } else if (y == "applarge") {
        lookup_options._append_largest_ring_size_if_no_match = 1;
        if (verbose) {
          cerr << "For molecules with missing rings, append size of largest ring\n";
        }
      } else if (y.starts_with("rpt=")) {
        uint32_t rpt;
        if (!y.numeric_value(rpt) || rpt == 0) {
          cerr << "Invalid report progress directive '" << y << "'\n";
          return 1;
        }
        report_progress.set_report_every(rpt);
        if (verbose) {
          cerr << "Will report progress every " << rpt << " molecules processed\n";
        }
      } else if (y == "scaffold") {
        ring_extractor.options().set_use_scaffolds(1);
        if (verbose) {
          cerr << "Will store and lookup scaffolds\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  if (store_in_database) {
    if (databases.empty()) {
      cerr << "No database specified for storage\n";
      return 1;
    }
    database_builder = std::make_unique<smi2rings_bdb::RingDatabaseBuilder>(
        *databases[0], build_options);
  }
  int rc;
  if (lookup_in_database) {
    rc = ProcessInputs<true>(cl, input_type, output);
  } else if (store_in_database) {
    rc = ProcessInputs<false>(cl, input_type, output);
  } else {
    cerr << "Specify a database build or lookup mode\n";
    return 1;
  }
  if (database_builder) {
    if (build_options._cache && verbose) {
      cerr << "Writing " << database_builder->pending_ring_systems()
           << " ring systems to the database\n";
    }
    if (!database_builder->Finish()) {
      rc = 1;
    }
  }

  output.flush();

  for (auto* database : databases) {
    database->close(0);
  }

  if (verbose) {
    cerr << "Produced rings for " << molecules_read << " molecules\n";

    for (int i = 0; i < ring_system_size_array.number_elements(); i++) {
      if (ring_system_size_array[i]) {
        cerr << ring_system_size_array[i] << " ring systems had " << i << " rings\n";
      }
    }

    cerr << "Ring systems had between " << ring_system_ring_stats.minval() << " and "
         << ring_system_ring_stats.maxval() << " rings";
    if (ring_system_ring_stats.n() > 1) {
      cerr << " average " << ring_system_ring_stats.average();
    }
    cerr << '\n';

    cerr << "Ring systems had between " << ring_system_atom_stats.minval() << " and "
         << ring_system_atom_stats.maxval() << " atoms";
    if (ring_system_atom_stats.n() > 1) {
      cerr << " average " << ring_system_atom_stats.average();
    }
    cerr << '\n';

    if (store_in_database) {
      cerr << "Stored " << database_builder->new_ring_systems() << " new ring systems\n";
    } else {
      cerr << "Found " << molecules_with_unique_ring_systems_found
           << " molecules with unique ring systems, "
           << ring_systems_not_found_in_database << " ring systems not in database\n";
      cerr << "Exemplar count values btw " << acc_exemplar.minval() << " and "
           << acc_exemplar.maxval() << " average "
           << static_cast<float>(acc_exemplar.average()) << "\n";
      cerr << "Note examplar counts are truncated at "
           << lookup_options._max_exemplar_count << '\n';
      for (int i = 0; i < exemplar_count.number_elements(); ++i) {
        if (exemplar_count[i]) {
          cerr << exemplar_count[i] << " rings had exemplar count " << i << '\n';
        }
      }
    }

    if (lookup_options._min_urs_needed > 0) {
      cerr << discarded_for_urs_too_low << " molecules discarded for URS below "
           << lookup_options._min_urs_needed << '\n';
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = smi2rings(argc, argv);

  return rc;
}
