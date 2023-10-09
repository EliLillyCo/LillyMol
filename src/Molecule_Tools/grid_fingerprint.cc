/*
  Computes atom pair fingerprints between ligand atoms and a grid
*/

#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static Atom_Typing_Specification atom_typing_specification;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString tag;

static std::vector<float> distance_cutoff;
static float longest_distance_cutoff;

static int write_3d_smiles = 0;

static Accumulator_Int<int> nset;

static IWString tag_for_each_grid_atoms;

static resizable_array_p<Substructure_Query> queries;

static int bonds_from_matched_atoms = 0;

static int molecules_matching_queries = 0;

static int ignore_molecules_not_matching_queries = 0;

static int matched_atom_bit_multiplier = 1;

// One way of fingerprinting is to have as many bits
// as there are in the grid. Every time an atom is
// within range of a grid point, increment that bit.
// That way we get a fingerprint that is how often a
// ligand atom is near each grid point.
static int fingerprint_grid_occupancy = 0;

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
  cerr << "Generates atom pair fingerprints between ligand atoms and grid atoms\n";
  cerr << " -T <dist>     produce a bit for each atom pair within radius <dist>. Can be multiple -T options\n";
  cerr << " -P ...        atom type specification for ligand\n";
  cerr << " -F <tag>      tag to create\n";
  cerr << " -w            write 3D smiles (contains coordinates)\n";
  cerr << " -s <smarts>   only fingerprint atoms matching <smarts>\n";
  cerr << " -x <nbonds>   used in conjunction with the -s option, extend <nbonds> from matched atoms\n";
  cerr << " -w <n>        increase importance (count) of matched atoms relative to unmatched atoms (n > 1)\n";
  cerr << " -z            ignore molecules not matching queries\n";
  cerr << " -l            reduce to largest fragment\n";
  cerr << " -M ...        more options, enter '-M help' for info\n";
  cerr << " -i <type>     input specification\n";
  cerr << " -g ...        chemical standardisation options\n";
  cerr << " -E ...        standard element specifications\n";
  cerr << " -A ...        standard aromaticity specifications\n";
  cerr << " -v            verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  When computing the fingerprint of the grid points, we need a means of keeping track of
  what kinds of ligand atoms are near each grid atom
*/

class Ligand_Atoms_Near_Grid_Point {
 private:
  std::unordered_map<unsigned int, unsigned int> _b;

  unsigned int _max_count;

 public:
  Ligand_Atoms_Near_Grid_Point();

  void extra(unsigned int b, unsigned int c) {
    _b[b] += c;
  }

  void extra(unsigned int b) {
    _b[b]++;
  }

  unsigned int max_count();
};

Ligand_Atoms_Near_Grid_Point::Ligand_Atoms_Near_Grid_Point() {
  _max_count = 0;

  return;
}

unsigned int
Ligand_Atoms_Near_Grid_Point::max_count() {
  if (_max_count > 0) {
    return _max_count;
  }

  _max_count = 0;

  for (auto i : _b) {
    if (i.second > _max_count) {
      _max_count = i.second;
    }
  }

  return _max_count;
}

static void
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
identify_ligand_atoms_to_process(Molecule& ligand,
                                 resizable_array_p<Substructure_Query>& queries,
                                 const int bonds_from_matched_atoms,
                                 int* process_these_ligand_atoms) {
  const int nq = queries.number_elements();

  if (0 == nq) {
    std::fill_n(process_these_ligand_atoms, ligand.natoms(), 1);
    return 1;
  }

  if (matched_atom_bit_multiplier > 1) {
    std::fill_n(process_these_ligand_atoms, ligand.natoms(), 1);
  } else {
    std::fill_n(process_these_ligand_atoms, ligand.natoms(), 0);
  }

  Molecule_to_Match target(&ligand);

  int matches_found = 0;

  for (int i = 0; i < nq; ++i) {
    Substructure_Results sresults;
    const int nhits = queries[i]->substructure_search(target, sresults);
    if (0 == nhits) {
      continue;
    }

    sresults.each_embedding_set_vector(process_these_ligand_atoms,
                                       matched_atom_bit_multiplier);

    matches_found++;
  }

  if (0 == matches_found) {
    return 0;
  }

  if (bonds_from_matched_atoms > 0) {
    const int matoms = ligand.natoms();

    for (int i = 0; i < matoms; ++i) {
      if (0 == process_these_ligand_atoms[i]) {
        continue;
      }

      for (int j = 0; j < matoms; ++j) {
        if (process_these_ligand_atoms[j]) {
          continue;
        }

        if (ligand.bonds_between(i, j) <= bonds_from_matched_atoms) {
          process_these_ligand_atoms[j] = 1;
        }
      }
    }
  }

  molecules_matching_queries++;

  return matches_found;
}

/*
*/

static int
grid_fingerprint(Molecule& ligand, const Molecule& grid, const int* grid_atom_type,
                 IWString_and_File_Descriptor& output) {
  const auto ligand_natoms = ligand.natoms();

  // cerr << "Ligand " << ligand.name() << " contains " << ligand_natoms << " atoms\n";

  int* ligand_atom_type = new int[ligand_natoms];
  std::unique_ptr<int[]> free_ligand_atom_type(ligand_atom_type);
  if (!atom_typing_specification.assign_atom_types(ligand, ligand_atom_type)) {
    cerr << "Cannot assign atom types '" << ligand.name() << "'\n";
    return 0;
  }

  const auto grid_natoms = grid.natoms();

  Ligand_Atoms_Near_Grid_Point* langp = nullptr;

  if (tag_for_each_grid_atoms.length() > 0) {
    langp = new Ligand_Atoms_Near_Grid_Point[grid_natoms];
  }

  int* process_these_ligand_atoms = new int[ligand_natoms];
  std::unique_ptr<int[]> free_process_these_ligand_atoms(process_these_ligand_atoms);

  if (identify_ligand_atoms_to_process(ligand, queries, bonds_from_matched_atoms,
                                       process_these_ligand_atoms)) {
    ;
  } else if (ignore_molecules_not_matching_queries) {
    return 1;
  } else {
    cerr << "None of " << queries.number_elements() << " queries matched "
         << ligand.name() << '\n';
    return 0;
  }

  Sparse_Fingerprint_Creator sfc;

  for (int i = 0; i < ligand_natoms; ++i) {
    if (!process_these_ligand_atoms[i]) {
      continue;
    }

    const Atom* a = ligand.atomi(i);

    for (auto j = 0; j < grid_natoms; ++j) {
      const auto d = a->distance(*(grid.atomi(j)));

      if (d >= longest_distance_cutoff) {
        continue;
      }

      for (unsigned int k = 0; k < distance_cutoff.size(); ++k) {
        if (d > distance_cutoff[k]) {
          break;
        }

        //      const unsigned int b = ligand_atom_type[i] * (5712 + k + 1) +
        //      grid_atom_type[j];    // very interesting that I used K here. Is that a
        //      good idea or not? Means you must always use the same set of -T options.
        //      May need to revisit this...
        if (fingerprint_grid_occupancy) {
          sfc.hit_bit(j, 1);
        } else {
          const unsigned int b = ligand_atom_type[i] * 5712 + grid_atom_type[j];
          sfc.hit_bit(b, process_these_ligand_atoms[i]);
        }

        if (nullptr != langp) {
          langp[j].extra(ligand_atom_type[i]);
        }
      }
    }
  }

  if (verbose) {
    nset.extra(sfc.nbits());
  }

  if (write_3d_smiles) {
    set_append_coordinates_after_each_atom(1);
  }

  output << smiles_tag << ligand.smiles() << ">\n";
  output << identifier_tag << ligand.name() << ">\n";

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);
  output << tmp << '\n';

  output << "|\n";

  output.write_if_buffer_holds_more_than(8192);

  return output.good();
}

static int
grid_fingerprint(data_source_and_type<Molecule>& input, const Molecule& grid,
                 const int* grid_atom_type, IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!grid_fingerprint(*m, grid, grid_atom_type, output)) {
      return 0;
    }
  }

  return 1;
}

static Molecule*
read_grid(const char* fname, FileType input_type) {
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  Molecule* grid = input.next_molecule();

  return grid;
}

static int
grid_fingerprint(const char* fname, FileType input_type, const Molecule& grid,
                 const int* grid_atom_type, IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return grid_fingerprint(input, grid, grid_atom_type, output);
}

static void
DisplayDashMOptions(std::ostream& output) {
  output << " -M grid        fingerprint grid occupancy - regardless of atom type\n";

  ::exit(0);
}

static int
grid_fingerprint(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lT:P:s:x:zw:M:");

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
  } else {
    set_global_aromaticity_type(Daylight);
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

  if (cl.option_present('M')) {
    IWString m;
    for (int i = 0; cl.value('M', m, i); ++i) {
      if (m == "grid") {
        fingerprint_grid_occupancy = 1;
        if (verbose) {
          cerr << "Will fingerprint grid occupancy\n";
        }
      } else if (m == "help") {
        DisplayDashMOptions(cerr);
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        DisplayDashMOptions(cerr);
      }
    }
  }

  if (cl.option_present('F')) {
    cl.value('F', tag);

    if (!tag.ends_with('<')) {
      tag << '<';
    }

    if (verbose) {
      cerr << "Fingerprints written with tag '" << tag << "'\n";
    }
  } else {
    tag = "NCGD<";
  }

  if (!cl.option_present('T')) {
    cerr << "Must specify inter molecular distance cutoff via the -T option\n";
    usage(1);
  } else {
    float t;
    for (auto i = 0; cl.value('T', t, i); ++i) {
      if (t <= 0.0f) {
        cerr << "The distance cutoff (-T) must be positive value\n";
        return 2;
      }
      distance_cutoff.push_back(t);
    }

    if (verbose) {
      cerr << "Defined " << distance_cutoff.size() << " distance cutoffs\n";
    }

    std::sort(distance_cutoff.begin(), distance_cutoff.end());

    longest_distance_cutoff = distance_cutoff[distance_cutoff.size() - 1];

    if (verbose) {
      for (unsigned int i = 0; i < distance_cutoff.size(); ++i) {
        cerr << "Bits set for atoms within " << distance_cutoff[i] << " distance\n";
      }
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Cannot parse ligand atom type specification '" << p << "'\n";
      return 1;
    }
  } else {
    const_IWSubstring p("none");

    atom_typing_specification.build(p);

    if (verbose) {
      cerr << "Using default 'none' ligand atom type\n";
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      Substructure_Query* q = new Substructure_Query;
      if (!q->create_from_smarts(s)) {
        cerr << "Cannot interpret smarts '" << s << "'\n";
        delete q;
        return 2;
      }

      queries.add(q);
    }

    if (verbose) {
      cerr << "Defined " << queries.number_elements() << " matched atom queries\n";
    }

    if (cl.option_present('x')) {
      if (!cl.value('x', bonds_from_matched_atoms) || bonds_from_matched_atoms < 0) {
        cerr << "The extend from matched atoms option (-x) must be a non-negative "
                "number\n";
        usage(1);
      }

      if (verbose) {
        cerr << "Will extend " << bonds_from_matched_atoms
             << " bonds away from matched atoms\n";
      }
    }

    if (cl.option_present('z')) {
      ignore_molecules_not_matching_queries = 1;

      if (verbose) {
        cerr << "Will ignore molecules not matching any of the queries\n";
      }
    }

    if (cl.option_present('w')) {
      if (!cl.value('w', matched_atom_bit_multiplier) ||
          matched_atom_bit_multiplier < 2) {
        cerr << "The matched atom bit multiplier option (-w) must be a whole +ve number, "
                "greater than 1\n";
        usage(1);
      }

      if (verbose) {
        cerr << "Will create " << matched_atom_bit_multiplier
             << " count multiplier for matched atoms\n";
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SDF;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() < 2) {
    cerr << "Must specify grid file and ligands file\n";
    usage(2);
  }

  Molecule* grid = read_grid(cl[0], input_type);

  if (nullptr == grid) {
    cerr << "Cannot read grid molecule '" << cl[0] << "'\n";
    return 2;
  }

  const auto grid_natoms = grid->natoms();

  if (verbose) {
    cerr << "Grid contains " << grid_natoms << " atoms\n";
  }

  int* grid_atom_type = new int[grid_natoms];
  std::unique_ptr<int[]> free_grid_atom_type(grid_atom_type);

  for (auto i = 0; i < grid_natoms; ++i) {
    grid_atom_type[i] = i;
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 1; i < cl.number_elements(); i++) {
    if (!grid_fingerprint(cl[i], input_type, *grid, grid_atom_type, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    if (queries.number_elements() > 0) {
      cerr << molecules_matching_queries << " molecules matched one or more queries\n";
    }
    if (nset.n() > 0) {
      cerr << "Fingerprints had between " << nset.minval() << " and " << nset.maxval()
           << " bits set, ave " << static_cast<float>(nset.average()) << "\n";
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = grid_fingerprint(argc, argv);

  return rc;
}
