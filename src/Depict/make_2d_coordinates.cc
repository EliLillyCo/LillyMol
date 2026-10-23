// Generate 2D coordinates for depiction.

#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"

#include "Depict/coords_2d.h"
#include "Depict/wedge_2d.h"

namespace make_2d_coordinates {

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
// clang-format off
  cerr << R"(Generates 2D coordinates for depiction.
Only .sdf output can hold the coordinates, so that is the default output type.
Molfiles are written with full V2000 counts lines, since the point of the
coordinates is that other software can read them, and LillyMol's default terse
counts line is rejected by some toolkits. Use -t to get LillyMol's terse form.
  make_2d_coordinates -S out file.smi
 -b <length>    bond length in the generated layout, default 1.5, the molfile convention
 -p quick       lay out quickly, accepting clashes that would otherwise be resolved
 -p standard    the default effort
 -p best        work hardest to resolve clashes, slowest
 -e             spread substituents evenly around an atom
 -k             skip the force field refinement that follows construction
 -n             do not centre the layout on the origin
 -t             write LillyMol's terse molfile counts line instead of full V2000
 -w             do not assign wedge bonds. Chirality is then written only as the
                MDL atom parity field, which most other software ignores
 -z             do not lay cis/trans double bonds out to match the molecule.
                The layout engine then draws every double bond trans, so Z
                double bonds are silently written as E
 -h             remove explicit hydrogens first. Usually wanted, they are laid
                out as atoms in their own right otherwise
 -c             discard molecules for which no layout could be generated,
                rather than treating that as a fatal error
 -l             reduce to the largest fragment
 -g ...         chemical standardisation options
 -S <fname>     output file name stem
 -i ...         input file specification
 -o ...         output file specification. Default .sdf
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_explicit_hydrogens = 0;

    Chemical_Standardisation _chemical_standardisation;

    depict::Coords2DOptions _coords_2d;

    uint64_t _molecules_read = 0;

    // Layout failure is rare enough that by default it is fatal. With -c the
    // molecule is skipped instead and counted here.
    int _ignore_layout_failures = 0;
    uint64_t _layout_failures = 0;

    // Layouts in which two atoms that are not bonded came out closer than
    // _crowded_threshold. These are written out - they are drawable, just
    // cramped - and only counted, for the verbose report.
    uint64_t _crowded_layouts = 0;
    float _crowded_threshold = 0.0f;

    // Wedge bonds are what carries chirality to other software, so they are on
    // by default. -w turns them off.
    int _assign_wedge_bonds = 1;

    // Double bonds whose E/Z the molecule stated, and how many of those the
    // layout drew that way. Only counted in verbose mode, since that is the only
    // mode that asks for a Coords2DResult.
    uint64_t _cis_trans_bonds = 0;
    uint64_t _cis_trans_honoured = 0;

    uint64_t _chiral_centres = 0;
    uint64_t _centres_wedged = 0;
    uint64_t _centres_not_stereogenic = 0;
    uint64_t _centres_unresolved = 0;

  public:
    int Initialise(Command_Line& cl);

    int
    verbose() const {
      return _verbose;
    }

    int Preprocess(Molecule& m);

    int Process(Molecule& m, Molecule_Output_Object& output);

    int Report(std::ostream& output) const;
};

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to the largest fragment\n";
    }
  }

  if (cl.option_present('h')) {
    _remove_explicit_hydrogens = 1;
    if (_verbose) {
      cerr << "Will remove explicit hydrogens\n";
    }
  }

  if (cl.option_present('c')) {
    _ignore_layout_failures = 1;
    if (_verbose) {
      cerr << "Will skip molecules that cannot be laid out\n";
    }
  }

  if (cl.option_present('b')) {
    float b;
    if (!cl.value('b', b) || b <= 0.0f) {
      cerr << "The bond length (-b) must be a positive number\n";
      return 0;
    }
    _coords_2d.bond_length = b;
    if (_verbose) {
      cerr << "Bond length set to " << b << '\n';
    }
  }

  if (cl.option_present('p')) {
    const IWString p = cl.string_value('p');
    if (p == "quick") {
      _coords_2d.precision = depict::Precision::kQuick;
    } else if (p == "standard") {
      _coords_2d.precision = depict::Precision::kStandard;
    } else if (p == "best") {
      _coords_2d.precision = depict::Precision::kBest;
    } else {
      cerr << "Unrecognised precision (-p) '" << p << "', must be one of "
           << "quick, standard or best\n";
      return 0;
    }
  }

  if (cl.option_present('e')) {
    _coords_2d.even_angles = true;
  }

  if (cl.option_present('k')) {
    _coords_2d.skip_minimization = true;
  }

  if (cl.option_present('n')) {
    _coords_2d.centre = false;
  }

  if (cl.option_present('w')) {
    _assign_wedge_bonds = 0;
  }

  if (cl.option_present('z')) {
    _coords_2d.honour_cis_trans = false;
    if (_verbose) {
      cerr << "Will not lay cis/trans double bonds out to match the molecule\n";
    }
  }

  // Half a bond length. Closer than this and the two atoms' labels overlap in
  // any reasonable rendering.
  _crowded_threshold = 0.5f * _coords_2d.bond_length;

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

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  // After standardisation, which can itself add or remove hydrogens.
  if (_remove_explicit_hydrogens) {
    m.remove_all(1);
  }

  if (m.empty()) {
    return 0;
  }

  return 1;
}

int
Options::Process(Molecule& m, Molecule_Output_Object& output) {
  ++_molecules_read;

  // Measuring how crowded the layout is costs time quadratic in the atom
  // count, so only ask for it when there is a verbose report to put it in.
  depict::Coords2DResult result;
  const int rc = _verbose ? depict::Generate2DCoordinates(m, _coords_2d, result)
                          : depict::Generate2DCoordinates(m, _coords_2d);

  if (!rc) {
    ++_layout_failures;
    cerr << "Cannot generate 2D coordinates for '" << m.name() << "'\n";
    return _ignore_layout_failures;
  }

  if (result.closest_nonbonded_approach >= 0.0f &&
      result.closest_nonbonded_approach < _crowded_threshold) {
    ++_crowded_layouts;
  }

  _cis_trans_bonds += result.cis_trans_bonds;
  _cis_trans_honoured += result.cis_trans_honoured;

  // Must follow layout: the direction a wedge points is decided from the
  // coordinates.
  if (_assign_wedge_bonds) {
    depict::WedgeResult wedges;
    depict::AssignWedgeBonds(m, wedges);
    _chiral_centres += wedges.chiral_centres;
    _centres_wedged += wedges.wedged;
    _centres_not_stereogenic += wedges.not_stereogenic;
    _centres_unresolved += wedges.unresolved;
  }

  return output.write(m);
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_layout_failures > 0) {
    output << _layout_failures << " molecules could not be laid out\n";
  }
  if (_crowded_layouts > 0) {
    output << _crowded_layouts << " layouts have non bonded atoms closer than "
           << _crowded_threshold << ", and will look crowded\n";
  }
  if (_cis_trans_bonds > 0) {
    output << _cis_trans_bonds << " cis/trans double bonds, " << _cis_trans_honoured
           << " drawn as the molecule says\n";
    if (_cis_trans_honoured < _cis_trans_bonds) {
      output << (_cis_trans_bonds - _cis_trans_honoured)
             << " double bonds are drawn as the other isomer, and will be read back\n"
                "as such\n";
    }
  }
  if (_chiral_centres > 0) {
    output << _chiral_centres << " chiral centres, " << _centres_wedged
           << " given a wedge bond\n";
    if (_centres_not_stereogenic > 0) {
      output << _centres_not_stereogenic
             << " marked chiral but not stereogenic, nothing to draw\n";
    }
    if (_centres_unresolved > 0) {
      output << _centres_unresolved
             << " stereocentres could not be wedged, their stereochemistry will\n"
                "not reach software that reads wedge bonds\n";
    }
  }

  return 1;
}

int
MakeCoordinates(Options& options, data_source_and_type<Molecule>& input,
                Molecule_Output_Object& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!options.Preprocess(*m)) {
      continue;
    }

    if (!options.Process(*m, output)) {
      return 0;
    }
  }

  return 1;
}

int
MakeCoordinates(Options& options, const char* fname, FileType input_type,
                Molecule_Output_Object& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "MakeCoordinates:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return MakeCoordinates(options, input, output);
}

int
MakeCoordinates(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:K:i:o:S:g:lhcb:p:ekntwz");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  // LillyMol writes an abbreviated molfile counts line by default, which some
  // toolkits refuse to read. Coordinates that only LillyMol can read are of no
  // use to a depiction, so ask for the full V2000 records. Done before the -o
  // options are processed so that anything given there still takes effect.
  if (!cl.option_present('t')) {
    set_write_isis_standard(1);
    set_write_mdl_charges_as_m_chg(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (!cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  Molecule_Output_Object output;
  if (!cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SDF);
  } else if (!output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    return 1;
  }

  IWString s = cl.string_value('S');
  if (output.would_overwrite_input_files(cl, s)) {
    cerr << "Cannot overwrite input file(s) with stem '" << s << "'\n";
    return 1;
  }

  if (!output.new_stem(s)) {
    cerr << "Cannot open stream for output '" << s << "'\n";
    return 1;
  }

  if (verbose) {
    cerr << "Output written to '" << s << "'\n";
  }

  for (const char* fname : cl) {
    if (!MakeCoordinates(options, fname, input_type, output)) {
      cerr << "MakeCoordinates::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace make_2d_coordinates

int
main(int argc, char** argv) {
  int rc = make_2d_coordinates::MakeCoordinates(argc, argv);

  return rc;
}
