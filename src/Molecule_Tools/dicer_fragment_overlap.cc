// Examines the output of dicer in order to identify the
// degree of overlap between fragments.
// Dicer must have been run with the `-I ini` option so
// that the isotopes are the initial atom numbers in the
// parent.

#include <iostream>
#include <memory>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_fragment_overlap {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Examines dicer proto output with the `-I ini` option to discern overlap between fragments\n";
  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    int _molecules_read = 0;

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(const const_IWSubstring& buffer,
                IWString_and_File_Descriptor& output);
    int Process(dicer_data::DicedMolecule& proto,
                 IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }
};

Options::Options() {
  _verbose = 0;
  _molecules_read = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  return 1;
}

int
Options::Process(const const_IWSubstring& buffer,
                 IWString_and_File_Descriptor& output) {
  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());
  dicer_data::DicedMolecule proto;
  if (! google::protobuf::TextFormat::Parse(&input, &proto)) {
    cerr << "Options::Process:cannot parse '" << buffer << "'\n";
    return 0;
  }

  return Process(proto, output);
}

// For every isotope, including zero, in `m` set
// the corresponding entry in `isotope`.
void
SetIsotopes(const Molecule& m, int * isotope) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);
    isotope[iso] = 1;
  }
}

// Return the number of isotopes present in `m` that are
// already set in `isotope`.
int
IsotopeOverlap(const Molecule& m, const int* isotope) {
  int result = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);
    if (isotope[iso] > 0) {
      ++result;
    }
  }

  return result;
}

// Determine the fragment overlap in `proto`.
// We add FragmentOverlap messages to `proto` and
// write to `output`.
int
Options::Process(dicer_data::DicedMolecule& proto,
                 IWString_and_File_Descriptor& output) {
  Molecule parent;
  if (! parent.build_from_smiles(proto.smiles())) {
    cerr << "Options::Process:cannot parse parent smiles '" << proto.smiles() << "'\n";
    return 0;
  }

  parent.set_name(proto.name());

  const int nfrags = proto.fragment_size();

  if (nfrags < 2) {
    // do something
    return 1;
  }

  std::unique_ptr<Molecule[]> fragments = std::make_unique<Molecule[]>(nfrags);
  for (int i = 0; i < nfrags; ++i) {
    const dicer_data::DicerFragment& frag = proto.fragment(i);
    if (! fragments[i].build_from_smiles(frag.smi())) {
      cerr << "Options::Process:cannot parse fragment smiles '" << frag.smi() << "'\n";
      return 0;
    }
  }

  const int matoms = parent.natoms();
  std::unique_ptr<int[]> isotope = std::make_unique<int[]>(matoms);

  for (int i = 0; i < nfrags; ++i) {
    std::fill_n(isotope.get(), matoms, 0);
    const uint32_t f1 = proto.fragment(i).id();
    SetIsotopes(fragments[i], isotope.get());
    for (int j = i + 1; j < nfrags; ++j) {
      int overlap = IsotopeOverlap(fragments[j], isotope.get());
      dicer_data::FragmentOverlap* ov = proto.add_overlap();
      const uint32_t f2 = proto.fragment(j).id();
      if (f1 < f2) {
        ov->set_f1(f1);
        ov->set_f2(f2);
      } else {
        ov->set_f1(f2);
        ov->set_f2(f1);
      }

      ov->set_ov(overlap);
    }
  }

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);
  std::string buffer;
  if (! printer.PrintToString(proto, &buffer)) {
    cerr << "Options::Process:cannot write '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  output << buffer;
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  return 1;
}

int
DicerFragmentOverlapSingle(Options& options,
                const const_IWSubstring& buffer,
                IWString_and_File_Descriptor& output) {
  return options.Process(buffer, output);
}

int
DicerFragmentOverlap(Options& options,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! DicerFragmentOverlapSingle(options, buffer, output)) {
      cerr << "DicerFragmentOverlap:error '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
DicerFragmentOverlap(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "DicerFragmentOverlap::cannot open '" << fname << "'\n";
    return 0;
  }

  return DicerFragmentOverlap(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "v");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! DicerFragmentOverlap(options, fname, output)) {
      cerr << "DicerFragmentOverlap::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_fragment_overlap

int
main(int argc, char ** argv) {

  int rc = dicer_fragment_overlap::Main(argc, argv);

  return rc;
}
