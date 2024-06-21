// Read a series of DicerFragment protos and divide them into types.
// The diced fragments must have an isotopic label.
// The protos are writen to separate files, depending on what kind
// of molecule the fragment is.
// If a fragment has just one isotope, we form a name like
//  <stem>1_<n>.textproto
// where <n> is the number of atoms in the fragment.
// If a fragment has two isotopes, we form a name like
// <stem>2_<n>.textproto
// where <n> is the distance between the isotopic atoms.
// If a fragment has three isotopes, we form a name like
// <stem>3_<i>.<j>.<k>.textproto
// where i,j,k are the (ordered) bonds between the attachment points,
// in other words, the sides of the triangle defined by the three atoms.

#include <iostream>
#include <memory>

#include "absl/strings/string_view.h"
#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_to_topological_types {

using std::cerr;
using iw_tf_data_record::TFDataWriter;
using dicer_data::DicerFragment;

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
cerr << R"(Processes DicerFragment textproto's and creates separate files based on the fragment.
 -S <stem>              stem for generated files
 -z                     ignore smiles interpretation errors
 -T                     write TFDataRecord serialized protos
 -3                     write 3 connected fragments - by defaint only 1 and 2 connected are written
 -v                     verbose output
)";

  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    IWString _stem;

    IWString _suffix;

    // Not sure how these will be used.
    int _write_3_connected_fragments;

    // As fragments are encountered, we create some number of files and
    // write the textprotos to its corresponding file.
    // I was not able to figure out how to do this without the unique_ptr
    // because IWString_and_File_Descriptor has no move or copy operators.
    IW_STL_Hash_Map<IWString, std::unique_ptr<IWString_and_File_Descriptor>> _output;

    int _write_tfdata_records;

    // If writing tfdata records, these are our output streams.
    IW_STL_Hash_Map<IWString, std::unique_ptr<TFDataWriter>> _tfdata_output;

    int _molecules_read = 0;

    int _ignore_smiles_interpretation_errors;
    int _smiles_errors_encountered;

  // private functions
    int Substituent(const const_IWSubstring& buffer,
                    const dicer_data::DicerFragment& proto, Molecule& m);
    int Linker2(const const_IWSubstring& buffer,
                    const dicer_data::DicerFragment& proto, Molecule& m);
    int Linker3(const const_IWSubstring& buffer,
                    const dicer_data::DicerFragment& proto, Molecule& m);
    int Process(const const_IWSubstring& buffer,
                const dicer_data::DicerFragment& proto);
    int Write(const const_IWSubstring& buffer, const dicer_data::DicerFragment& proto, IWString& fname);
    int WriteTFData(const dicer_data::DicerFragment& proto,
               IWString& fname);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    int Process(const const_IWSubstring& buffer);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _molecules_read = 0;
  _write_3_connected_fragments = 0;
  _write_tfdata_records = 0;
  _ignore_smiles_interpretation_errors = 0;
  _smiles_errors_encountered = 0;
  _suffix = ".textproto";
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify the output file name stem via the -S option\n";
    Usage(1);
  }

  cl.value('S', _stem);
  if (_verbose) {
    cerr << "Files will be created with name stem " << _stem << "'\n";
  }

  if (cl.option_present('z')) {
    _ignore_smiles_interpretation_errors = 1;
    if (_verbose) {
      cerr << "Will ignore smiles interpretation errors\n";
    }
  }

  if (cl.option_present('3')) {
    _write_3_connected_fragments = 1;
    if (_verbose) {
      cerr << "Will write three connected fragments\n";
    }
  }

  if (cl.option_present('T')) {
    _write_tfdata_records = 1;
    if (_verbose) {
      cerr << "Will write TFDataRecord serialized protos\n";
    }
    _suffix = ".tfdata";
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_write_tfdata_records) {
    output << _tfdata_output.size();
  } else {
    output << _output.size();
  }
  output << " files created\n";

  output << _smiles_errors_encountered << " smiles errors encountered\n";

  return 1;
}

int
Options::Process(const const_IWSubstring& buffer) {
  ++_molecules_read;

  dicer_data::DicerFragment proto;
  const absl::string_view tmp(buffer.data(), buffer.length());
  if (!google::protobuf::TextFormat::ParseFromString(tmp, &proto)) {
    cerr << "Options::Process:cannot parse text proto " << buffer << '\n';
    return 0;
  }

  return Process(buffer, proto);
}

int
Options::Process(const const_IWSubstring& buffer,
                 const dicer_data::DicerFragment& proto) {
  Molecule m;
  if (! m.build_from_smiles(proto.smi()) || m.empty()) {
    cerr << "Invalid smiles " << proto.ShortDebugString() << '\n';
    ++_smiles_errors_encountered;
    return _ignore_smiles_interpretation_errors;
  }

  const int niso = m.number_isotopic_atoms();
  if (niso == 0) {
    cerr << "Molecule with no isotopic atoms " << proto.smi() << '\n';
    return 0;
  }

  switch (niso) {
    case 1:
      return Substituent(buffer, proto, m);
    case 2:
      return Linker2(buffer, proto, m);
    case 3:
      return Linker3(buffer, proto, m);
    default:
      cerr << "Molecules with " << niso << " attachments not processed " << proto.smi() << '\n';
      return 1;
  }
}

int
Options::Write(const const_IWSubstring& buffer,
               const dicer_data::DicerFragment& proto,
               IWString& fname) {

  if (_write_tfdata_records) {
    return WriteTFData(proto, fname);
  }

  auto iter = _output.find(fname);
  if (iter == _output.end()) {
    std::unique_ptr<IWString_and_File_Descriptor> tmp = std::make_unique<IWString_and_File_Descriptor>();
    if (! tmp->open(fname)) {
      cerr << "Write::cannot open '" << fname << '\n';
      return 0;
    }
    *tmp << buffer << '\n';
    _output.emplace(std::make_pair(fname, std::move(tmp)));
  } else {
    *(iter->second) << buffer << '\n';
    iter->second->write_if_buffer_holds_more_than(4192);
  }

  return 1;
}

int
Options::WriteTFData(const dicer_data::DicerFragment& proto,
               IWString& fname) {

  assert(_write_tfdata_records);

  auto iter = _tfdata_output.find(fname);
  if (iter == _tfdata_output.end()) {
    std::unique_ptr<TFDataWriter> tmp = std::make_unique<TFDataWriter>();
    if (! tmp->Open(fname)) {
      cerr << "WriteTFData::cannot open '" << fname << '\n';
      return 0;
    }
    tmp->WriteSerializedProto<dicer_data::DicerFragment>(proto);
    _tfdata_output.emplace(std::make_pair(fname, std::move(tmp)));
  } else {
    iter->second->WriteSerializedProto<dicer_data::DicerFragment>(proto);
  }

  return 1;
}
               

int
Options::Substituent(const const_IWSubstring& buffer,
                     const dicer_data::DicerFragment& proto,
                     Molecule& m) {
  const int matoms = m.natoms();
  IWString fname(_stem);
  fname << "_1_" << matoms << _suffix;

  return Write(buffer, proto, fname);
}

int
DistanceBetweenIsotopicAtoms(Molecule& m) {
  atom_number_t iso1 = kInvalidAtomNumber;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i) == 0) {
      continue;
    }

    if (iso1 == kInvalidAtomNumber) {
      iso1 = i;
    } else {
      return m.bonds_between(iso1, i);
    }
  }

  cerr << "No pair of isotopes in " << m.smiles() << '\n';
  return -1;
}

int
Options::Linker2(const const_IWSubstring& buffer,
                 const dicer_data::DicerFragment& proto,
                 Molecule& m) {
  const int d = DistanceBetweenIsotopicAtoms(m);
  if (d < 0) {
    return 0;
  }

  IWString fname(_stem);
  fname << "_2_" << d << _suffix;
  return Write(buffer, proto, fname);
}

int
OrderedDistances(int& d1, int& d2, int& d3) {
  if (d1 < d2 && d2 < d3) {
    return 1;
  }

  if (d1 > d2) {
    std::swap(d1, d2);
  }

  if (d2 > d3) {
    std::swap(d2, d3);
  }

  if (d1 > d2) {
    std::swap(d1, d2);
  }

  return 1;
}

int
GetDistances(Molecule& m, int& d1, int& d2, int& d3) {
  atom_number_t iso1 = kInvalidAtomNumber;
  atom_number_t iso2 = kInvalidAtomNumber;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i) == 0) {
      continue;
    }
    if (iso1 == kInvalidAtomNumber) {
      iso1 = i;
    } else if (iso2 == kInvalidAtomNumber) {
      iso2 = i;
    } else {
      d1 = m.bonds_between(iso1, iso2);
      d2 = m.bonds_between(iso1, i);
      d3 = m.bonds_between(iso2, i);

      return OrderedDistances(d1, d2, d3);
    }
  }

  cerr << "Did not find 3 isotopes in " << m.smiles() << '\n';

  return -1;
}
             

int
Options::Linker3(const const_IWSubstring& buffer,
                 const dicer_data::DicerFragment& proto,
                 Molecule& m) {
  if (! _write_3_connected_fragments) {
    return 1;
  }

  int d1, d2, d3;
  if (! GetDistances(m, d1, d2, d3)) {
    return 0;
  }

  IWString fname(_stem);
  fname << "_3_" << d1 << '.' << d2 << '.' << d3 << _suffix;

  return Write(buffer, proto, fname);
}

int
DicerToTopologicalTypes(Options& options,
                iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! options.Process(buffer)) {
      cerr << "Error processing " << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
DicerToTopologicalTypes(Options& options,
             const char * fname) {

  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DicerToTopologicalTypes:cannot open '" << fname << "'\n";
    return 0;
  }

  return DicerToTopologicalTypes(options, input);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:S:zT");

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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }


  for (const char * fname : cl) {
    if (! DicerToTopologicalTypes(options, fname)) {
      cerr << "DicerToTopologicalTypes::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_to_topological_types

int
main(int argc, char ** argv) {

  int rc = dicer_to_topological_types::Main(argc, argv);

  return rc;
}

