// Adds centroid atoms to 3D molecules according to substructure query matches.

#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Lib/target.h"

namespace meaningful_name {

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
  cerr << R"(Adds centroid atoms to 3D molecules according to substructure query matches
 -C <ele>       use <ele> as the centroid atom, default 'Ce'
 -R arom        place a centroid atom at the centre of each aromatic ring
 -R aromsys     place a centroid atom at the centre of each aromatic ring system : NOT IMPLEMENTED. DO NOT USE.
 -s <smarts>    smarts describing matched atoms.
 -q <query>     query describing matched atoms.
 -z i           ignore molecules not matching any of the queries.
 -S <fname>     ouput file.
 -i ...         input file specification.
 -o ...         output file specification. Default .sdf. Use '-i smi -o smi3d' for 3d smiles.
 -v             verbose output
)";
// clang-format on

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

    uint64_t _molecules_read = 0;

    resizable_array_p<Substructure_Query> _query;

    const Element* _centroid_element;

    int _aromatic_ring_centroid;
    int _aromatic_ring_system_centroid;

    int _ignore_molecules_not_matching_query;
    uint64_t _molecules_not_matching_queries;

  // Private functions.
    int AromaticRingCentroid(Molecule& m, Molecule_Output_Object& output);
    int AromaticRingSystemCentroid(Molecule& m, Molecule_Output_Object& output);
    int PlaceCentroid(Molecule& m, const Set_of_Atoms& embedding);

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

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, Molecule_Output_Object& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;

  _centroid_element = nullptr;
  _aromatic_ring_centroid = 0;
  _aromatic_ring_system_centroid = 0;

  _ignore_molecules_not_matching_query = 0;
  _molecules_not_matching_queries = 0;
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

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('C')) {
    IWString ele = cl.string_value('C');
    _centroid_element = get_element_from_symbol_no_case_conversion(ele);
    if (_centroid_element == nullptr) {
      cerr << "Invalid centroid element specification '" << ele << "'\n";
      return 1;
    }
    if (_verbose) {
      cerr << "Centroid element " << _centroid_element->symbol() << '\n';
    }
  } else {
    _centroid_element = get_element_from_atomic_number(58);
  }

  if (cl.option_present('R') &&
      (cl.option_present('s') || cl.option_present('q'))) {
    cerr << "The ring centroid option (-R) is not compatible with the -s/-q options\n";
    Usage(1);
  }

  if (cl.option_present('R')) {
    IWString r;
    for (int i = 0; cl.value('R', r, i); ++i) {
      if (r == "arom") {
        _aromatic_ring_centroid = 1;
      } else if (r == "aromsys") {
        _aromatic_ring_system_centroid = 1;
      } else {
        cerr << "Unrecognised -R qualifier '" << r << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matching_query = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching queries\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 1;
      }
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(smarts)) {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _query << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _query, _verbose, 'q')) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  for (Substructure_Query* q : _query) {
    q->set_find_unique_embeddings_only(1);
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_ignore_molecules_not_matching_query) {
    cerr << _molecules_not_matching_queries << " molecules did not match any of the queries\n";
  }

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

int
Options::Process(Molecule& m,
                 Molecule_Output_Object& output) {
  ++_molecules_read;

  m.compute_aromaticity_if_needed();

  if (_aromatic_ring_centroid) {
    return AromaticRingCentroid(m, output);
  }

  if (_aromatic_ring_system_centroid) {
    return AromaticRingSystemCentroid(m, output);
  }

  Molecule_to_Match target(&m);

  int rc = 0;

  Molecule mcopy(m);

  for (Substructure_Query* q : _query) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      rc += PlaceCentroid(mcopy, *e);
    }
  }

  if (rc == 0) {
    ++_molecules_not_matching_queries;
    return _ignore_molecules_not_matching_query;
  }

  return output.write(mcopy);
}

int
Options::AromaticRingCentroid(Molecule& m,
                              Molecule_Output_Object& output) {
  int rc = 0;

  Molecule mcopy(m);

  for (const Ring* r : m.sssr_rings()) {
    if (! r->is_aromatic()) {
      continue;
    }

    rc += PlaceCentroid(mcopy, *r);
  }

  if (rc == 0) {
    ++_molecules_not_matching_queries;
    return _ignore_molecules_not_matching_query;
  }

  return output.write(mcopy);
}

int
Options::AromaticRingSystemCentroid(Molecule& m,
                              Molecule_Output_Object& output) {
  cerr << "Options:;AromaticRingCentroid:not implemented, see Ian\n";
  return 0;
}

int
Options::PlaceCentroid(Molecule& m,
                       const Set_of_Atoms& embedding) {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  // cerr << "PlaceCentroid " << embedding << '\n';
  for (atom_number_t a : embedding) {
    const Atom& atom = m[a];
    x += atom.x();
    y += atom.y();
    z += atom.z();
  }

  const int matoms = m.natoms();

  m.add(_centroid_element);
  const double ring_size = static_cast<double>(embedding.size());
  m.setxyz(matoms, x / ring_size, y / ring_size, z / ring_size);

  return 1;
}

int
ApplicationName(Options& options,
                Molecule& m,
                Molecule_Output_Object& output) {
  return options.Process(m, output);
}

int
ApplicationName(Options& options,
                data_source_and_type<Molecule>& input,
                Molecule_Output_Object& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! ApplicationName(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             FileType input_type,
             Molecule_Output_Object& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ApplicationName(options, input, output);
}

int
ApplicationName(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:C:R:S:i:o:s:q:z:g:");

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

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  Molecule_Output_Object output;
  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SDF);
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    return 1;
  }

  if (cl.option_present('S')) {
    IWString s = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, s)) {
      cerr << "Cannot overwrite input file(s) with stem '" << s << "'\n";
      return 0;
    }

    if (! output.new_stem(s)) {
      cerr << "Cannot open stream for output '" << s << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Output written to '" << s << "'\n";
    }
  }

  for (const char * fname : cl) {
    if (! ApplicationName(options, fname, input_type, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace meaningful_name

int
main(int argc, char ** argv) {

  int rc = meaningful_name::ApplicationName(argc, argv);

  return rc;
}
