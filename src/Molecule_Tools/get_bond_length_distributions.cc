// Scan a set of 3D molecules and discern the bond length distributions,
// and write a TFDataRecord file of bond_length_distribution::BondLengthDistribution
// protos with the distributions.

#include <iostream>

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/bond_length_distribution.pb.h"
#else
#include "bond_length_distribution.pb.h"
#endif

namespace bond_length_distribution {

using std::cerr;

// By convention the Usage function tells how to use the tool.
void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << R"(Extract bond length distributions from a set of 3D molecules.
Writes a TFDataRecord file of bond_length_distribution::BondLengthDistribution protos.
 -X <distance>          maximum bond length distance to consider (def 2)
 -S <fname>             name of output file
 -v                     verbose output
)";

  ::exit(rc);
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    float _min_distance;
    uint32_t _bonds_shorter_than_min_distance;

    float _max_distance;
    uint32_t _bonds_longer_than_max_distance;

    absl::flat_hash_map<uint32_t, bond_length_distribution::BondLengthDistribution> _bonds;

    uint64_t _number_aromatic_distances = 0;
    uint64_t _number_single_distances = 0;
    uint64_t _number_double_distances = 0;
    uint64_t _number_triple_distances = 0;

  // Private functions

    uint32_t Bucket(float dist) const;
    int Write(iw_tf_data_record::TFDataWriter& writer) const;

    int WriteDistribution(const bond_length_distribution::BondLengthDistribution& proto,
                IWString_and_File_Descriptor& output) const;
    int WriteDistribution(const BondLengthDistribution& proto, const IWString& stem) const;

    void AddBond(const Bond& b, uint32_t bucket, bond_length_distribution::BondLengthDistribution& proto);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    int Write(IWString& fname);
    int WriteDistributions(const IWString& stem) const;

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _min_distance = 1.0;
  _bonds_shorter_than_min_distance = 0;
  _bonds_longer_than_max_distance = 0;
  _max_distance = 2.0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name via the -S option\n";
    Usage(1);
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _bonds.size() << " different connection types\n";
  output << _number_aromatic_distances << " number aromatic bonds\n";
  output << _number_single_distances << " number single bonds\n";
  output << _number_double_distances << " number double bonds\n";
  output << _number_triple_distances << " number triple bonds\n";
  // Other information about what has happened.

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

// This hash function only works for atomic numbers less than 10
uint32_t
FormKey(atomic_number_t z1, atomic_number_t z2) {
  assert(z1 < 10);
  assert(z2 < 10);

  if (z1 < z2) {
    std::swap(z1, z2);
  }

  return 10 * z1 + z2;
}

// convert a distance into an array index.
uint32_t
Options::Bucket(float dist) const {
  assert(dist >= _min_distance);

  dist -= _min_distance;
  return static_cast<uint32_t> (dist * 100.0f);
}

void
Options::AddBond(const Bond& b, uint32_t bucket, bond_length_distribution::BondLengthDistribution& proto) {
  if (b.is_aromatic()) {
    ++(*proto.mutable_aromatic())[bucket];
    ++_number_aromatic_distances;
    return;
  }

  if (b.is_single_bond()) {
    ++(*proto.mutable_single())[bucket];
    ++_number_single_distances;
    return;
  }

  if (b.is_double_bond()) {
    ++(*proto.mutable_double_())[bucket];
    ++_number_double_distances;
    return;
  }

  if (b.is_triple_bond()) {
    ++(*proto.mutable_triple())[bucket];
    ++_number_triple_distances;
    return;
  }
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  m.compute_aromaticity_if_needed();

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    const Atom& at1 = m[a1];
    const Atom& at2 = m[a2];

    atomic_number_t z1 = at1.atomic_number();
    atomic_number_t z2 = at2.atomic_number();

    // We don't need to adjust the atoms since we are just dealing with the distance.
    if (z1 < z2) {
      std::swap(a1, a2);
    }

    float dist = m.distance_between_atoms(a1, a2);
    if (dist < _min_distance) {
      ++_bonds_shorter_than_min_distance;
      continue;
    }
    if (dist > _max_distance) {
      ++_bonds_longer_than_max_distance;
      continue;
    }

    uint32_t bucket = Bucket(dist);

    uint32_t key = FormKey(z1, z2);

    auto iter = _bonds.find(key);
    if (iter == _bonds.end()) {
      bond_length_distribution::BondLengthDistribution proto;
      proto.set_atomic_number1(z1);
      proto.set_atomic_number2(z2);
      proto.mutable_single()->Resize(101, 0);
      proto.mutable_double_()->Resize(101, 0);
      proto.mutable_triple()->Resize(101, 0);
      proto.mutable_aromatic()->Resize(101, 0);
      AddBond(*b, bucket, proto);
      _bonds[key] = std::move(proto);
    } else {
      AddBond(*b, bucket, iter->second);
    }
  }

  return 1;
}

int
Options::Write(IWString& fname) {
  iw_tf_data_record::TFDataWriter writer;
  if (! writer.Open(fname)) {
    cerr << "Options::Write:cannot open '" << fname << "'\n";
    return 0;
  }

  return Write(writer);
}

int
Options::Write(iw_tf_data_record::TFDataWriter& writer) const {
  if (_verbose) {
    cerr << "Writing " << _bonds.size() << " bond length distributions\n";
  }

  for (const auto& [key, value] : _bonds) {
    writer.WriteSerializedProto<BondLengthDistribution>(value);
    if (_verbose) {
      cerr << value.ShortDebugString() << '\n';
    }
  }

  return 1;
}

IWString
FileName(const IWString& stem,
         atomic_number_t z1,
         atomic_number_t z2) {
  IWString result(stem);
  result << '.' << z1 << '.' << z2 << ".txt";

  return result;
}

int
Options::WriteDistributions(const IWString& stem) const {
  for (const auto& [key, proto] : _bonds) {
    WriteDistribution(proto, stem);
  }

  return 1;
}

int
Options::WriteDistribution(const BondLengthDistribution& proto, const IWString& stem) const {
  IWString fname = FileName(stem, proto.atomic_number1(), proto.atomic_number2());
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteDistributions:cannot open '" << fname << "'\n";
    return 0;
  }

  if (! WriteDistribution(proto, output)) {
    cerr << "Options::WriteDistributions:cannot write '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Options::WriteDistribution(const bond_length_distribution::BondLengthDistribution& proto,
                IWString_and_File_Descriptor& output) const {
  static constexpr char kSep = ' ';

  output.reserve(32768);

  output << "Dist" << kSep << "Aromatic" << kSep << "NAromatic" << kSep <<
            "Single" << kSep << "NSingle" << kSep <<
            "Double" << kSep << "NDouble" << kSep <<
            "Triple" << kSep << "NTriple" << '\n';

  uint64_t narom = 0;
  uint64_t nsingle = 0;
  uint64_t ndouble = 0;
  uint64_t ntriple = 0;

  for (int i = 0; i <  proto.single().size(); ++i) {
    narom += proto.aromatic(i);
    nsingle += proto.single(i);
    ndouble += proto.double_(i);
    ntriple += proto.triple(i);
  }

  for (int i = 0; i <  proto.single().size(); ++i) {
    float d = 1.0f + 0.01 * i;
    output << d;
    
    if (proto.aromatic(i) == 0) {
      output << kSep << '0' << kSep << '0';
    } else {
      float f = iwmisc::Fraction<float>(proto.aromatic(i), narom);
      output << kSep << proto.aromatic(i) << kSep << f;
    }

    if (proto.single(i) == 0) {
      output << kSep << '0' << kSep << '0';
    } else {
      float f = iwmisc::Fraction<float>(proto.single(i), nsingle);
      output << kSep << proto.single(i) << kSep << f;
    }

    if (proto.double_(i) == 0) {
      output << kSep << '0' << kSep << '0';
    } else {
      float f = iwmisc::Fraction<float>(proto.double_(i), ndouble);
      output << kSep << proto.double_(i) << kSep << f;
    }

    if (proto.triple(i) == 0) {
      output << kSep << '0' << kSep << '0';
    } else {
      float f = iwmisc::Fraction<float>(proto.triple(i), ntriple);
      output << kSep << proto.triple(i) << kSep << f;
    }

    output << '\n';
  }

  return 1;
}

int
GetBondLengthDistributions(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
GetBondLengthDistributions(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! GetBondLengthDistributions(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
GetBondLengthDistributions(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "GetBondLengthDistributions:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return GetBondLengthDistributions(options, input, output);
}

int
GetBondLengthDistributions(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:X:S:P:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  // Not used here.
  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GetBondLengthDistributions(options, fname, input_type, output)) {
      cerr << "GetBondLengthDistributions::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  IWString fname = cl.string_value('S');

  if (! options.Write(fname)) {
    cerr << "Cannot write '" << fname << "'\n";
    return 1;
  }

  if (cl.option_present('P')) {
    IWString stem;
    cl.value('P', stem);
    options.WriteDistributions(stem);
  }

  return 0;
}

}  // namespace bond_length_distribution

int
main(int argc, char ** argv) {

  int rc = bond_length_distribution::GetBondLengthDistributions(argc, argv);

  return rc;
}
