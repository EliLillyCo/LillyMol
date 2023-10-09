// Score molecules with respect to their position relative to a grid.
// The molecules must all have been pre-aligned.

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/insight_grid.h"

namespace grid_proximity {

using std::cerr;

IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");

void
Usage(int rc) {
  cerr << "Computes a proximity score of molecules wrt a precomputed grid\n";
  cerr << " -G <fname>     file containing InsightII grid\n";
  cerr << " -t <dist>      count density within <dist> of matched atoms in the ligands\n";
  cerr << " -q <qry>       query specifying ligand atoms to score\n";
  cerr << " -s <smarts>    smarts specifying ligand atoms to score\n";
  cerr << " -z i           ignore molecules not matching the query\n";
  cerr << " -M <fname>     write the grid as a molecule\n";
  cerr << " -m <float>     multiple grid raw values by <float> and place isotopes\n";
  cerr << " -i ...         input specification(s)\n";
  cerr << " -A ...         aromaticity options\n";
  cerr << " -E ...         element related options\n";
  ::exit(rc);
}

class GridProximityConfig {
  private:
    int _verbose = 0;

    FileType _input_type = FILE_TYPE_INVALID;

    int _reduce_to_largest_fragment = 0;

    Chemical_Standardisation _chemical_standardisation;

    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    // As read from the -G option.
    insight_grid::InsightGrid _grid;

    // The -t option.
    float _distance;

    resizable_array_p<Substructure_Query> _queries;

    int _no_query_match = 0;

    int _ignore_no_query_match = 0;

    // the number of molecules that are close to a grid point and
    // therefore generate a non zero result.
    int _molecules_with_non_zero_values;

    // We keep track of the maximum value across all molecules.
    Accumulator<double> _acc_max;

    // output file token separator.
    static constexpr char kSep = ' ';

    Molecule_Output_Object _stream_for_grid_as_molecule;

    IWString _missing;

    // We can generate fingerrints.
    int _generate_fingerprint = 0;
    // If fixed size, they will be the same number of points as the grid.
    int _generate_fixed_size_fingerprint = 0;
    int _generate_sparse_fingerprint = 0;
    // The tag for any fingerprint generated.
    IWString _tag;

  // private functions.

    int IdentifyAtoms(Molecule& m, int * matched);
    int Process(Molecule& m, const int * atoms,
                IWString_and_File_Descriptor& output);
    int WriteGridAsMolecule(Molecule& m);

    int GenerateFingerprint(Molecule& m,
                        const int* in_grid,
                        IWString_and_File_Descriptor& output);

  public:
    GridProximityConfig();
    ~GridProximityConfig() {
      if (_stream_for_grid_as_molecule.active()) {
        _stream_for_grid_as_molecule.do_close();
      }
    }

    int Initialise(Command_Line& cl);

    int WriteHeader(IWString_and_File_Descriptor& output) const;

    int Preprocess(Molecule& m);

    int MaybeDiscernInputType(const char * fname);

    FileType input_type() const {
      return _input_type;
    }

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }
};

GridProximityConfig::GridProximityConfig() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _distance = 0.0f;
  _input_type = FILE_TYPE_INVALID;
  _molecules_read = 0;
  _molecules_with_non_zero_values = 0;

  _no_query_match = 0;
  _ignore_no_query_match = 0;

  _missing = '.';
}

int
GridProximityConfig::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(1);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (! cl.option_present('G')) {
    cerr << "Must specify name of grid file via the -G option\n";
    Usage(1);
  } else {
    const char* fname = cl.option_value('G');
    if (! _grid.Build(fname)) {
      cerr << "GridProximityConfig::Initialise:cannot read grid '" << fname << "'\n";
      return 0;
    }
  }

  if (_verbose) {
    cerr << "Grid contains " << _grid.npoints() << " points\n";
  }

  if (cl.option_present('t')) {
    if (! cl.value('t', _distance) || _distance <= 0.0f) {
      cerr << "GridProximityConfig::Initialise:invalid distance (-t)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will examine grid points within " << _distance << " of matched ligand atoms\n";
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Options::Initialise:cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _queries << qry.release();
    }
  }

  if (_queries.empty()) {
    cerr << "GridProximityConfig::Initialise:no query atoms specified\n";
    return 0;
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_no_query_match = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching the queries\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    if (! _tag.ends_with('<')) {
      _tag << '<';
    }
    if (_tag.starts_with("FP")) {
      _generate_fixed_size_fingerprint = 1;
    } else if (_tag.starts_with("NC")) {
      _generate_sparse_fingerprint = 1;
    } else {
      cerr << "Unrecognised fingerprint type '" << _tag << "'\n";
      return 0;
    }
    _generate_fingerprint = 1;
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.option_present('M')) {
    const const_IWSubstring fname = cl.string_value('M');
    _stream_for_grid_as_molecule.add_output_type(FILE_TYPE_SDF);
    if (! _stream_for_grid_as_molecule.new_stem(fname)) {
      cerr << "Cannot open molecule grid file '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will write the grid as molecules to '" << fname << "'\n";
    }

    if (cl.option_present('m')) {
      float m;
      if (! cl.value('m', m) || m < 0.0f) {
        cerr << "The -m option must be >0\n";
        return 0;
      }
      _grid.set_isotope_multiplier(m);
      if (_verbose) {
        cerr << "Will place isotopes on grid points, multipler " << m << '\n';
      }
    }
  }

  return 1;
}

int
GridProximityConfig::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "ID" << kSep
         << "NonZero" << kSep
         << "NValues" << kSep
         << "Minval" << kSep
         << "Mean" << kSep
         << "Maxval" << kSep
         << "Sum" << kSep
         << '\n';

  return 1;
}

int
GridProximityConfig::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _no_query_match << " molecules did not match the query\n";
  output << _molecules_with_non_zero_values << " molecules had non zero values\n";
  output << "Across " << _acc_max.n() << " calculations, max field interaction btw "
         << _acc_max.minval() << " and " << _acc_max.maxval()
         << " mean " << _acc_max.average() << '\n';
  return 1;
}

int
GridProximityConfig::MaybeDiscernInputType(const char * fname) {
  if (_input_type == FILE_TYPE_INVALID) {
    _input_type = discern_file_type_from_name(fname);
  }
  return 1;
}

int
GridProximityConfig::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
GridProximityConfig::Process(Molecule& m,
                             IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  const int matoms = m.natoms();
  std::unique_ptr<int[]> atoms(new_int(matoms));
  if (! IdentifyAtoms(m, atoms.get())) {
    cerr << "GridProximityConfig::Process:cannot identify atoms to process\n";
   ++_no_query_match;
   if (_ignore_no_query_match) {
     return 1;
   }
    return 0;
  }

  return Process(m, atoms.get(), output);
}

int
GridProximityConfig::Process(Molecule& m,
                             const int * atoms,
                             IWString_and_File_Descriptor& output) {
  const int n = _grid.npoints();
  std::unique_ptr<int[]> in_grid(new_int(n));

  const int matoms = m.natoms();
  int atoms_with_nonzero_values = 0;
  for (int i = 0; i < matoms; ++i) {
    if (atoms[i] == 0) {
      continue;
    }
    const Atom * a = m.atomi(i);
    if (_grid.IndicesInRange(a->x(), a->y(), a->z(), _distance, in_grid.get())) {
      ++atoms_with_nonzero_values;
    }
  }

  if (atoms_with_nonzero_values) {
    ++_molecules_with_non_zero_values;
  }

  if (_generate_fingerprint) {
    return GenerateFingerprint(m, in_grid.get(), output);
  }

  Accumulator<double> acc;
  _grid.Score(in_grid.get(), acc);

  output << m.name() << kSep << atoms_with_nonzero_values << kSep
         << acc.n() << kSep;
  if (acc.n() == 0) {
    //output << _missing << kSep << _missing << kSep << _missing << kSep << _missing << '\n';
    output << _missing << kSep;
    output << _missing << kSep;
    output << _missing << kSep;
    output << _missing << kSep;
    output << '\n';
  } else {
//  output << acc.minval() << kSep;
    output << static_cast<float>(acc.average()) << kSep;
    output << acc.maxval() << kSep;
    output << acc.sum() << '\n';
    _acc_max.extra(acc.maxval());
  }

  output.write_if_buffer_holds_more_than(4096);

  cerr << "WRiting molecule gric " << _stream_for_grid_as_molecule.active() << '\n';
  if (_stream_for_grid_as_molecule.active()) {
    WriteGridAsMolecule(m);
  }

  return 1;
}

int
GridProximityConfig::GenerateFingerprint(Molecule& m,
                        const int* in_grid,
                        IWString_and_File_Descriptor& output) {
  const uint32_t npoints = _grid.npoints();

  output << smiles_tag << m.smiles() << ">\n";
  output << identifier_tag << m.name() << ">\n";

  if (_generate_fixed_size_fingerprint) {
    fixed_bit_vector::FixedBitVector bits(npoints);
    for (uint32_t i = 0; i < npoints; ++i) {
      if (in_grid[i]) {
        bits.set_bit(i);
      }
    }
    const IWString fp = bits.DaylightAsciiRepresentationIncludingNsetInfo();
    output << _tag << fp << ">\n";
  } else if (_generate_sparse_fingerprint) {
    Sparse_Fingerprint_Creator sfc;
    for (uint32_t i = 0; i < npoints; ++i) {
      if (in_grid[i]) {
        sfc.hit_bit(i, 1);
      }
    }
    IWString fp;
    sfc.daylight_ascii_form_with_counts_encoded(fp);
    output << _tag << fp << ">\n";
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

#ifdef COUNTS_ATOMS_MULTIPLE_TIMES
int
GridProximityConfig::Process(Molecule& m,
                             const int * atoms,
                             IWString_and_File_Descriptor& output) {
  Accumulator<double> acc;

  const int matoms = m.natoms();
  int atoms_with_nonzero_values = 0;
  for (int i = 0; i < matoms; ++i) {
    if (atoms[i] == 0) {
      continue;
    }
    const Atom * a = m.atomi(i);
    if (_grid.Sum(a->x(), a->y(), a->z(), _distance, acc)) {
      ++atoms_with_nonzero_values;
    }
  }

  if (atoms_with_nonzero_values) {
    ++_molecules_with_non_zero_values;
  }

  output << m.name() << kSep << atoms_with_nonzero_values << kSep
         << acc.n() << kSep;
  if (acc.n() == 0) {
    output << _missing << kSep << _missing << kSep << _missing << kSep << _missing << '\n';
  } else {
    output << acc.minval() << kSep;
    output << static_cast<float>(acc.average()) << kSep;
    output << acc.maxval() << kSep;
    output << acc.sum() << '\n';
    _acc_max.extra(acc.maxval());
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}
#endif

// Set `matched` for the matched atoms in the first query that matches `m`.
// Return 1 if any query matches, 0 otherwise.
int
GridProximityConfig::IdentifyAtoms(Molecule& m, int * matched) {
  Molecule_to_Match target(&m);

  Substructure_Results sresults;
  for (Substructure_Query* q : _queries) {
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    sresults.each_embedding_set_vector(matched, 1);
    return 1;
  }

  return 0;
}

int
GridProximityConfig::WriteGridAsMolecule(Molecule& m) {
  // Chose a single letter element, Iodine probably the most compact
  const Element * iodine = get_element_from_atomic_number(53);
  Molecule grid_mol;
  _grid.MakeMolecule(iodine, _distance, m, grid_mol);

  grid_mol.add_molecule(&m);
  cerr << "Combined grid + molecule has " << grid_mol.natoms() << " atoms\n";
  return _stream_for_grid_as_molecule.write(grid_mol);
}

int
GridProximity(GridProximityConfig& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
GridProximity(GridProximityConfig& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! GridProximity(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
GridProximity(GridProximityConfig& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.good()) {
    cerr << "GridProximity:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return GridProximity(options, input, output);
}

int
GridProximity(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:l:G:g:q:s:i:t:z:M:m:J:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  set_default_iwstring_float_concatenation_precision(4);
  set_default_iwstring_double_concatenation_precision(4);

  GridProximityConfig options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  options.WriteHeader(output);

  for (const char * fname : cl) {
    if (! GridProximity(options, fname, output)) {
      cerr << "GridProximity::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace grid_proximity

int
main(int argc, char ** argv) {

  int rc = grid_proximity::GridProximity(argc, argv);

  return rc;
}
