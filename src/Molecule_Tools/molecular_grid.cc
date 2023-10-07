/*
  Creates a grid of dummy atoms that can encompas a set of input molecules
*/

#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static Report_Progress report_progress;

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
  cerr << "  -d <dist>     spacing between grid atoms\n";
  cerr << "  -x <dist>     ensure grid expands <dist> from any atom (default 2A)\n";
  cerr << "  -t            trim the grid to only those atoms with -x of a ligand atom\n";
  cerr << "  -j            for each ligand, compute total and min grid point occupancy\n";
  cerr << "  -r <number>   during the -t option, report progress every <number> ligands examined\n";
  cerr << "  -S <fname>    write grid to <fname>\n";
  cerr << "  -o <type>     output type (default sdf)\n";
  cerr << "  -e <ele>      element to use for the grid\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static void
place_isotopes_on_grid_atoms_within_range(const Molecule& ligand, Molecule& grid,
                                          const coord_t extra_distance)
{
  const auto ligand_matoms = ligand.natoms();
  const auto grid_natoms = grid.natoms();

  if (report_progress()) {
    cerr << report_progress.times_called() << " grid coverage determination processing '"
         << ligand.name() << "'\n";
  }

  for (auto i = 0; i < ligand_matoms; ++i) {
    const Atom* a = ligand.atomi(i);

    for (auto j = 0; j < grid_natoms; ++j) {
      if (a->closer_than(*(grid.atomi(j)), extra_distance)) {
        grid.increment_isotope(j, 1);
      }
    }
  }
}

static int
grid_intensity_of_nearby_atoms(const Molecule& ligand, const Molecule& grid,
                               const coord_t extra_distance, unsigned int& sum,
                               unsigned int& minval)
{
  const auto ligand_natoms = ligand.natoms();
  const auto grid_natoms = grid.natoms();

  sum = 0;
  minval = numeric_limits<isotope_t>::max();

  if (report_progress()) {
    cerr << report_progress.times_called() << " ligand score processing '"
         << ligand.name() << "'\n";
  }

  for (auto i = 0; i < ligand_natoms; ++i) {
    const Atom* a = ligand.atomi(i);

    for (auto j = 0; j < grid_natoms; ++j) {
      if (!a->closer_than(*(grid.atomi(j)), extra_distance)) {
        continue;
      }

      const isotope_t iso = grid.isotope(j);
      sum += iso;
      if (iso < static_cast<isotope_t>(minval)) {
        minval = iso;
      }
    }
  }

  return sum;
}

static int
molecular_grid(Molecule& m, float& xmin, float& xmax, float& ymin, float& ymax,
               float& zmin, float& zmax)
{
  if (1 == molecules_read) {
    m.spatial_extremeties(xmin, xmax, ymin, ymax, zmin, zmax);
    return 1;
  }

  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; i++) {
    const auto x = m.x(i);
    const auto y = m.y(i);
    const auto z = m.z(i);

    if (x < xmin) {
      xmin = x;
    } else if (x > xmax) {
      xmax = x;
    }

    if (y < ymin) {
      ymin = y;
    } else if (y > ymax) {
      ymax = y;
    }

    if (z < zmin) {
      zmin = z;
    } else if (z > zmax) {
      zmax = z;
    }
  }

  return 1;
}

static int
molecular_grid(data_source_and_type<Molecule>& input,
               resizable_array_p<Molecule>& molecules, float& xmin, float& xmax,
               float& ymin, float& ymax, float& zmin, float& zmax)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    preprocess(*m);

    molecules.add(m);

    if (!molecular_grid(*m, xmin, xmax, ymin, ymax, zmin, zmax)) {
      return 0;
    }
  }

  return 1;
}

static int
molecular_grid(const char* fname, FileType input_type,
               resizable_array_p<Molecule>& molecules, float& xmin, float& xmax,
               float& ymin, float& ymax, float& zmin, float& zmax)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return molecular_grid(input, molecules, xmin, xmax, ymin, ymax, zmin, zmax);
}

static int
molecular_grid(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lS:o:d:x:O:e:tr:j");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  coord_t inter_atom_grid_spacing = static_cast<coord_t>(1.0);

  if (cl.option_present('d')) {
    if (!cl.value('d', inter_atom_grid_spacing) || inter_atom_grid_spacing <= 0.0f) {
      cerr << "The inter atom grid spacing (-d) must be a positive value\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Grid produced with inter atom spacing " << inter_atom_grid_spacing << endl;
    }
  }

  coord_t extra_distance = static_cast<coord_t>(2.0);
  if (cl.option_present('x')) {
    if (!cl.value('x', extra_distance) || extra_distance <= 0.0f) {
      cerr << "The extra distance option (-x) must be a positive value\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will extend the grid to " << extra_distance
           << " beyond the most extreme atom\n";
    }
  }

  if (!cl.option_present('S')) {
    cerr << "Must specify output file name via the -S option\n";
    usage(2);
  }

  Molecule_Output_Object output;

  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      return 2;
    }
  } else {
    output.add_output_type(FILE_TYPE_SDF);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s)) {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 2;
    }

    if (!output.new_stem(s)) {
      cerr << "gack, cannot open output file '" << s << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Output written to '" << s << "'\n";
    }
  }

  coord_t xmin = -numeric_limits<coord_t>::max();
  coord_t xmax = numeric_limits<coord_t>::max();
  coord_t ymin = -numeric_limits<coord_t>::max();
  coord_t ymax = numeric_limits<coord_t>::max();
  coord_t zmin = -numeric_limits<coord_t>::max();
  coord_t zmax = numeric_limits<coord_t>::max();

  resizable_array_p<Molecule> molecules;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!molecular_grid(cl[i], input_type, molecules, xmin, xmax, ymin, ymax, zmin,
                        zmax)) {
      cerr << "Fatal error processing grid for '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << "Atoms range X: " << xmin << ',' << xmax << '\n';
    cerr << "            Y: " << ymin << ',' << ymax << '\n';
    cerr << "            Z: " << zmin << ',' << zmax << '\n';
  }

  if (extra_distance > 0.0f) {
    xmin -= extra_distance;
    xmax += extra_distance;
    ymin -= extra_distance;
    ymax += extra_distance;
    zmin -= extra_distance;
    zmax += extra_distance;
  }

  coord_t xorigin = 0.0;
  coord_t yorigin = 0.0;
  coord_t zorigin = 0.0;

  if (cl.option_present('O')) {
    const_IWSubstring o = cl.string_value('O');

    if (3 != o.nwords(',')) {
      cerr << "The origin specification must be of the form 'x,z,y', '" << o
           << "' invalid\n";
      return 2;
    }

    int i = 0;
    const_IWSubstring token;
    o.nextword(token, i, ',');
    if (0 == token.length() || !token.numeric_value(xorigin)) {
      cerr << "INvalid grid origin specification '" << o << "'\n";
      return 2;
    }
    o.nextword(token, i, ',');
    if (0 == token.length() || !token.numeric_value(yorigin)) {
      cerr << "INvalid grid origin y specification '" << o << "'\n";
      return 2;
    }
    o.nextword(token, i, ',');
    if (0 == token.length() || !token.numeric_value(zorigin)) {
      cerr << "INvalid grid origin z specification '" << o << "'\n";
      return 2;
    }
  }

  const auto xrange = xmax - xmin;
  const auto yrange = ymax - ymin;
  const auto zrange = zmax - zmin;

  if (xrange <= inter_atom_grid_spacing || yrange <= inter_atom_grid_spacing ||
      zrange <= inter_atom_grid_spacing) {
    cerr << "Inter atom grid spacing " << inter_atom_grid_spacing
         << " too small for range of atoms\n";
    return 1;
  }

  const Element* grid_element = get_element_from_atomic_number(20);

  if (cl.option_present('e')) {
    const_IWSubstring a = cl.string_value('e');

    grid_element = get_element_from_symbol_no_case_conversion(a);

    if (nullptr == grid_element) {
      cerr << "Cannot fetch element '" << a << "'\n";
      return 0;
    }
  }

  // First fill up the rectangular grid

  uint64_t npoints_estimate =
      static_cast<uint64_t>((xmax - xmin) / inter_atom_grid_spacing) *
      static_cast<uint64_t>((ymax - ymin) / inter_atom_grid_spacing) *
      static_cast<uint64_t>((zmax - zmin) / inter_atom_grid_spacing);
  if (verbose) {
    cerr << "Estimate " << npoints_estimate << " points in grid\n";
  }

  Molecule grid;

  for (auto x = xmin; x <= xmax; x += inter_atom_grid_spacing) {
    for (auto y = ymin; y <= ymax; y += inter_atom_grid_spacing) {
      for (auto z = zmin; z <= zmax; z += inter_atom_grid_spacing) {
        grid.add(grid_element);
        grid.setxyz(grid.natoms() - 1, x, y, z);
      }
    }
  }

  const auto grid_natoms = grid.natoms();

  if (verbose) {
    cerr << "Grid contains " << grid_natoms << " atoms\n";
  }

  if (cl.option_present('t')) {
    if (cl.option_present('r')) {
      if (!report_progress.initialise(cl, 'r', verbose)) {
        return 2;
      }
    }

    for (auto i = 0; i < molecules.number_elements(); ++i) {
      const Molecule* m = molecules[i];
      place_isotopes_on_grid_atoms_within_range(*m, grid, extra_distance);
    }

    int atoms_removed = 0;

    for (int i = grid.natoms() - 1; i >= 0; --i) {
      if (0 == grid.isotope(i)) {
        grid.remove_atom(i);
        atoms_removed++;
      }
    }

    if (verbose) {
      cerr << "Trimming to proximity to ligands removed " << atoms_removed
           << " atoms from the grid, final grid molecule contains " << grid.natoms()
           << " atoms\n";
    }
  }

  grid.set_name("GRID");

  output.write(grid);

  if (cl.option_present('j')) {
    cout << "ID NATOMS SUM MINVAL\n";
    for (auto i = 0; i < molecules.number_elements(); ++i) {
      const Molecule* m = molecules[i];

      unsigned int sum = 0;
      unsigned int minval = numeric_limits<int>::max();
      grid_intensity_of_nearby_atoms(*m, grid, extra_distance, sum, minval);
      cout << m->name() << ' ' << m->natoms() << ' ' << sum << ' ' << minval << '\n';
    }
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = molecular_grid(argc, argv);

  return rc;
}
