// Read a series of 3D conformers with isotopic labels indicating
// attachment points.
// Create 3 files.
//   A file containing distances between attachment points.
//   A file containing bond angles
//   A file containing torsions.
// Each file is sorted.

#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace conf2indices {

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
// clang-format off
  cerr << R"(Given 3D isotopically labelled molecules, generates files containing sorted distances.
 -S <fname>  file name stem for output files (mandatory).
 -v          verbose output.
)";
// clang-format on

  ::exit(rc);
}

// All 3 files consist of these structs.

struct DistId {
  float dist;
  uint32_t id;
};

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

    // Upon first call, we count the number of molecules in the input.
    uint64_t _molecules_in_input = 0;

    uint64_t _molecules_read = 0;

    int _ignore_molecules_with_no_attachment_points = 0;

    IWString_and_File_Descriptor _distances_stream, _angles_stream, _torsions_stream;

    DistId* _distances;
    DistId* _angles;
    DistId* _torsions;
    

  // Private functions

    int DoOpen(const IWString& stem, const char* extension,
                IWString_and_File_Descriptor& output);
    int TwoAttachmentPoints(Molecule& m, const Set_of_Atoms& iso);
    int ThreeAttachmentPoints(Molecule& m, const Set_of_Atoms& iso);
    int OneAttachmentPoints(Molecule& m, atom_number_t zatom);

  public:
    Options();
    ~Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    int ignore_molecules_with_no_attachment_points() const {
      return _ignore_molecules_with_no_attachment_points;
    }

    // When the input file is opened, we know how many molecules are in the input
    // stream, so the DistId arrays can be allocated.
    int Resize(uint64_t s);

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int Process(Molecule& mol);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;

    // Once data is accumulated, write.
    int Write();
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;

  _distances = nullptr;
  _angles = nullptr;
  _torsions = nullptr;
}

Options::~Options() {
  if (_distances != nullptr) {
    delete [] _distances;
  }
  if (_angles != nullptr) {
    delete [] _angles;
  }
  if (_torsions != nullptr) {
    delete [] _torsions;
  }
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
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
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  if (cl.option_present('S')) {
    IWString stem = cl.string_value('S');
    if (! DoOpen(stem, "distances", _distances_stream)) {
      return 0;
    }
    if (! DoOpen(stem, "angles", _angles_stream)) {
      return 0;
    }
    if (! DoOpen(stem, "torsions", _torsions_stream)) {
      return 0;
    }
  }

  if (cl.option_present('z')) {
    _ignore_molecules_with_no_attachment_points = 1;
    if (_verbose) {
      cerr << "Will ignore molecules containing non attachment points\n";
    }
  }

  return 1;
}

int
Options::DoOpen(const IWString& stem, const char* extension,
                IWString_and_File_Descriptor& output) {
  IWString fname;
  fname << stem << ',' << extension;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "Options::DoOpen:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Options::Resize(uint64_t s) {
  _distances = new DistId[s];
  _angles = new DistId[s];
  _torsions = new DistId[s];

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
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
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return 1;
}

// Note that we do not differentiate anything about the connectivity
// of the attachment points - maybe in a later iteration.
int
Options::TwoAttachmentPoints(Molecule& m, const Set_of_Atoms& iso) {
  float d = m.distance_between_atoms(iso[0], iso[1]);
  if (d < 1.0f) {
    return 0;
  }

  DistId& di = _distances[_molecules_read - 1];
  di.dist = d;
  di.id = _molecules_read - 1;

  return 1;
}

// Not implemented.
int
Options::ThreeAttachmentPoints(Molecule& m, const Set_of_Atoms& iso) {
  return 1;
}

// Nothing to do with these.
int
Options::OneAttachmentPoints(Molecule& m, atom_number_t zatom) {
  return 1;
}

int
Options::Process(Molecule& m) {
  ++_molecules_read;

  Set_of_Atoms iso;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i) > 0) {
      iso << i;
    }
  }

  const int niso = iso.number_elements();
  if (niso == 2) {
    return TwoAttachmentPoints(m, iso);
  }
  if (niso == 3) {
    return ThreeAttachmentPoints(m, iso);
  }
  if (niso == 1) {
    return OneAttachmentPoints(m, iso[0]);
  }

  if (_ignore_molecules_with_no_attachment_points) {
    return 1;
  }

  cerr << m.smiles() << ' ' << m.name() << " no isotopes\n";
  return 0;
}

int
Options::Write() {
  if (_molecules_read == 0) {
    cerr << "Options::Write:no molecules\n";
    return 0;
  }

  std::sort(_distances, _distances + _molecules_read, [] (const DistId& di1, const DistId& di2) {
    return di1.dist < di2.dist;
  });

  int fd = _distances_stream.fd();

  size_t to_write = _molecules_read * sizeof(DistId);
  size_t written = IW_FD_WRITE(fd, _distances, to_write);
  if (written != to_write) {
    cerr << "Options::Write:wrote " << written << " of " << to_write << " bytes\n";
    return 0;
  }

  return 1;
}

int
ApplicationName(Options& options,
                data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! options.Process(*m)) {
      return 0;
    }
  }

  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  uint64_t n = input.molecules_remaining();
  if (n == 0) {
    cerr << "ApplicationName:no molecules in input\n";
    return 0;
  }

  if (! options.Resize(n)) {
    cerr << "ApplicationName:cannot initialise for " << n << " molecules\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ApplicationName(options, input);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vcg:lS:z:");

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
  
  if (cl.size() > 1) {
    cerr << "Cannot process multiple files\n";
    return 1;
  }

  for (const char * fname : cl) {
    if (! ApplicationName(options, fname, input_type)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  if (verbose) {
    cerr << "Writing...\n";
  }

  options.Write();

  return 0;
}

}  // namespace conf2indices

int
main(int argc, char ** argv) {

  int rc = conf2indices::Main(argc, argv);

  return rc;
}
