/*
  Extracts subsets of molecules
*/

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/target.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int invert_subset = 0;

static int max_matches_to_process = 0;

static int embedding_to_process = -1;

static int ignore_queries_not_hitting = 0;

static int write_molecules_not_matching_queries = 0;

static IWString append_to_unchanged_molecules;

static int molecules_read = 0;

static int write_each_hit_as_separate_molecule = 0;

static int isotope_at_break_points = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static extending_resizable_array<int> atoms_in_molecules_written;
static extending_resizable_array<int> atoms_removed;

static int create_subsets_by_bond = 0;

static int make_implicit_hydrogens_explicit = 0;

// This seems very strange, since we are using a Molecule_Output_Object
// for output. Seems wrong.
static int blank_line_to_cout_for_non_matches = 0;

// Aug 2004. Chavali Krishna wanted to preserve counterions

static int also_write_small_fragments = 0;

// Given a matched atom, include all other atoms that are within
// this distance.
static distance_t include_atoms_within = static_cast<distance_t>(0.0);

// We can extract all atoms that are within range of a given 
// point in space. This only works for molecules that are spatially
// aligned.

class SpatialPosition {
  private:
    Coordinates _xyz;

    float _radius;
  public:
    // The input must be x,y,z,radius
    int Build(const const_IWSubstring& buffer);

    // Return true of `coords` is within _radius of _xyz
    int Matches(const Coordinates& coords) const;

    // Examine the atoms in `embedding` and remove any that do not
    // match our position.
    // Returns true if any atoms removed.
    int RemoveNonMatchingAtoms(const Molecule& m, Set_of_Atoms& embedding) const;
};

int
SpatialPosition::Build(const const_IWSubstring& buffer) {
  int i = 0;
  const_IWSubstring token;
  for (int ndx = 0; buffer.nextword(token, i, ','); ++ndx) {
    float d;
    if (! token.numeric_value(d)) {
      cerr << "SpatialPosition::Build:invalid value '" << token << "'\n";
      return 0;
    }
    if (ndx == 0) {
      _xyz.set_x(d);
    } else if (ndx == 1) {
      _xyz.set_y(d);
    } else if (ndx == 2) {
      _xyz.set_z(d);
    } else if (ndx == 3) {
      if (d < 0.0f) {
        cerr << "SpatialPosition::Build:invalid radius '" << token << "'\n";
        return 0;
      }
      _radius = d;
    } else {
      cerr << "SpatialPosition::Build:invalid input '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
SpatialPosition::RemoveNonMatchingAtoms(const Molecule& m, Set_of_Atoms& embedding) const {
  int rc = 0;
  for (int i = embedding.number_elements() - 1; i >= 0; --i) {
    const Atom* a = m.atomi(embedding[i]);
    if (! Matches(*a)) {
      embedding.remove_item(i);
      ++rc;
    }
  }

  return rc;  // The number of items removed.
}

int
SpatialPosition::Matches(const Coordinates& coords) const {
  float d = _xyz.distance(coords);
  return d <= _radius;
}

static resizable_array_p<SpatialPosition> spatial_position;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << "usage: " << prog_name << " -i <input type> -o <output type> file1 file2...\n";
  cerr << "Creates subsets of molecules defined by query matches\n";
  cerr << "  -q <query>     specify substructure query\n";
  cerr << "  -s <smarts>    specify smarts (ss)\n";
  cerr << "  -n             write the molecule subset which does not match the query\n";
  cerr << "  -b             create the subsets by bond rather than by atom\n";
  cerr << "  -m <number>    only process the first <number> hits for any query\n";
  cerr << "  -m do=<nn>     process embedding number <nn> (starts with 0)\n";
  cerr << "  -m each        process each embedding separately\n";
  cerr << "  -z i           ignore molecules not matching any query\n";
  cerr << "  -z w           write molecules not matching any query\n";
  cerr << "  -z b           blank line to stdout for non-matching molecules\n";
  cerr << "  -z nmapp=xxx   append 'xxx' to the names of non-matching molecules\n";
  cerr << "  -f             also write small fragments with selected atoms\n";
  cerr << "  -T <dist>      with 3D input, include atoms within <dist> of any matched "
          "atom\n";
  cerr << "  -C x,y,z,rad   with 3D input, only search atoms within <rad> of <x,y,z>\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -M ...         various miscellaneous options, enter '-M help' for info\n";
  cerr << "  -E ...         element specifications, enter '-E help' for details\n";
  cerr << "  -A ...         standard aromaticity options, enter '-A help' for details\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int
identify_fragments_with_matched_atoms(Molecule& m, const Set_of_Atoms& embedding,
                                      int* hits_in_fragment) {
  if (nullptr == hits_in_fragment) {
    return 1;
  }

  int n = embedding.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = embedding[i];

    int f = m.fragment_membership(j);

    hits_in_fragment[f]++;
  }

  return 1;
}

static int
add_back_atoms_in_small_fragments(Molecule& m, int* remove_atom,
                                  const int* hits_in_fragment) {
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (0 == remove_atom[i]) {  // atom being kept
      continue;
    }

    int f = m.fragment_membership(i);

    if (0 == hits_in_fragment[f]) {  // atom is in fragment that was not matched by query
      remove_atom[i] = 0;
    }
  }

  return 1;
}

static int
do_write_subset(const Molecule& m, const int* include_bond, const int* hits_in_fragment,
                Molecule_Output_Object& output) {
  Molecule subset;

  m.create_subset_by_bond(subset, include_bond, 1);

  if (0 == subset.natoms()) {
    cerr << "Warning, empty molecule created from '" << m.name() << "'\n";
    return 0;
  }

  return output.write(subset);
}

static int
identify_bonds_to_include(Molecule& m, const Set_of_Atoms& e, int* include_bond,
                          int* fragment_used = nullptr) {
  int n = e.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++) {
    atom_number_t ai = e[i];

    const Atom* atom_ai = m.atomi(ai);

    for (int j = i + 1; j < n; j++) {
      atom_number_t aj = e[j];

      const Bond* b = atom_ai->bond_to_atom(aj);

      if (nullptr == b) {
        continue;
      }

      int bn = b->bond_number();

      assert(bn >= 0);

      if (include_bond[bn]) {
        continue;
      }

      include_bond[bn] = 1;

      if (nullptr != fragment_used) {
        fragment_used[m.fragment_membership(aj)] = 1;
      }

      rc++;
    }
  }

  return rc;
}

static int
write_bond_subsets(Molecule& m, int nhits, const Substructure_Results& sresults,
                   int* include_bond, const int* hits_in_fragment,
                   Molecule_Output_Object& output) {
  set_vector(include_bond, m.nedges(), 0);

  for (int i = 0; i < nhits; i++) {
    const Set_of_Atoms* e = sresults.embedding(i);

    m.convert_set_of_atoms_to_bond_numbers(*e, include_bond);
  }

  return do_write_subset(m, include_bond, hits_in_fragment, output);
}

/*
  What do we do if no queries match
*/

static int
handle_no_queries_matching(Molecule& m, Molecule_Output_Object& output) {
  if (ignore_queries_not_hitting) {
    ;
  } else {
    cerr << "No query hits to '" << m.name() << "'\n";
    return 0;
  }

  if (write_molecules_not_matching_queries) {
    if (append_to_unchanged_molecules.length()) {
      m.append_to_name(append_to_unchanged_molecules);
    }

    return output.write(m);
  }

  if (blank_line_to_cout_for_non_matches) {
    std::cout << '\n';
  }

  return output.good();
}

static int
molecule_subset_by_bond(Molecule& m, int* hits_in_fragment,
                        Molecule_Output_Object& output) {
  m.assign_bond_numbers_to_bonds();

  int nq = queries.number_elements();

  Molecule_to_Match target(&m);

  int queries_matching = 0;

  int* include_bond = new_int(m.nedges(), 0);
  std::unique_ptr<int[]> free_include_bond(include_bond);

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      if (verbose > 1) {
        cerr << "Zero hits to query " << i << " '" << queries[i]->comment()
             << "', only matched " << sresults.max_query_atoms_matched_in_search()
             << " query atoms\n";
      }

      continue;
    }

    queries_matching++;

    if (max_matches_to_process > 0 && nhits > max_matches_to_process) {
      nhits = max_matches_to_process;
    }

    if (write_each_hit_as_separate_molecule) {
      write_bond_subsets(m, nhits, sresults, include_bond, hits_in_fragment, output);
      continue;
    }

    if (embedding_to_process < 0) {
      ;
    } else if (embedding_to_process < nhits) {
      const Set_of_Atoms* e = sresults.embedding(embedding_to_process);
      identify_bonds_to_include(m, *e, include_bond);
      identify_fragments_with_matched_atoms(m, *e, hits_in_fragment);
      continue;
    } else {
      cerr << "Requested embedding " << embedding_to_process << " but only " << nhits
           << " matches\n";
      return 0;
    }

    for (int i = 0; i < nhits; i++) {
      const Set_of_Atoms* e = sresults.embedding(i);

      identify_bonds_to_include(m, *e, include_bond);
      identify_fragments_with_matched_atoms(m, *e, hits_in_fragment);
    }
  }

  if (verbose > 1) {
    cerr << queries_matching << " of " << nq << " queries hit '" << m.name() << "'\n";
  }

  if (0 == queries_matching) {
    return handle_no_queries_matching(m, output);
  }

  if (verbose > 1) {
    cerr << queries_matching << " of " << nq << " matched '" << m.name() << "'\n";
  }

  if (write_each_hit_as_separate_molecule) {  //  output already done
    return output.good();
  }

  return do_write_subset(m, include_bond, hits_in_fragment, output);
}

static int
do_isotope_at_break_points(Molecule& m, const int* remove_atom) {
  int nb = m.nedges();

  atom_number_t a1, a2;  // scope here for efficiency

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    a1 = b->a1();
    a2 = b->a2();

    if (remove_atom[a1] == remove_atom[a2]) {  // either both kept or both removed
      continue;
    }

    if (remove_atom[a1]) {
      m.set_isotope(a2, isotope_at_break_points);
    } else if (remove_atom[a2]) {
      m.set_isotope(a1, isotope_at_break_points);
    }
  }

  return 1;
}

// #define DEBUG_REMOVE_ATOMS

static int
do_output(Molecule& m, int* remove_atom, const int* hits_in_fragment,
          Molecule_Output_Object& output) {
  assert(nullptr != remove_atom);

#ifdef DEBUG_REMOVE_ATOMS
  for (int i = 0; i < m.natoms(); i++) {
    if (remove_atom[i]) {
      cerr << "Will remove atom " << i << ", '" << m.smarts_equivalent_for_atom(i)
           << "'\n";
    }
  }
#endif

  if (nullptr != hits_in_fragment) {
    add_back_atoms_in_small_fragments(m, remove_atom, hits_in_fragment);
  }

  int number_atoms_removed = count_non_zero_occurrences_in_array(remove_atom, m.natoms());

  if (isotope_at_break_points) {
    do_isotope_at_break_points(m, remove_atom);
  }

  if (number_atoms_removed) {
    m.remove_atoms(remove_atom);
  } else {
    cerr << "Warning, no atoms removed\n";
  }

  int matoms = m.natoms();

  atoms_in_molecules_written[matoms]++;
  atoms_removed[number_atoms_removed]++;

  return output.write(m);
}

static void
identify_atoms_to_be_removed(const Set_of_Atoms& e, int* remove_atom) {
  if (invert_subset) {
    e.set_vector(remove_atom, 1);
  } else {
    e.set_vector(remove_atom, 0);
  }

  return;
}

static int
identify_atoms_proximal_to_retained_atoms(const Molecule& m, int* remove_atom,
                                          const distance_t include_atoms_within) {
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms;
       ++i)  // ensure that only the number 1 is used in the remove_atom array
  {
    if (remove_atom[i]) {
      remove_atom[i] = 1;
    }
  }

  int rc = 0;

  for (auto i = 0; i < matoms; ++i) {
    if (remove_atom[i]) {  // start with atoms that were initially NOT being removed
      continue;
    }

    const Atom* a = m.atomi(i);

    for (auto j = 0; j < matoms; ++j) {
      if (2 == remove_atom[j] || 0 == remove_atom[j]) {
        continue;
      }

      if (j == i) {
        continue;
      }

      if (!a->closer_than(*(m.atomi(j)), include_atoms_within)) {
        continue;
      }

      remove_atom[j] = 2;  // we will switch this later
      rc++;
    }
  }

  if (rc) {
    for (auto i = 0; i < matoms; ++i) {
      if (2 == remove_atom[i]) {
        remove_atom[i] = 0;
      }
    }
  }

  return rc;
}

int
OkSpatialConstraints(const Molecule& m,
                     Set_of_Atoms& embedding,
                     const resizable_array_p<SpatialPosition>& position) {
  for (const SpatialPosition* p : position) {
    if (p->RemoveNonMatchingAtoms(m, embedding) == 0) {
      continue;
    }
    if (embedding.empty()) {
      return 0;
    }
  }

  return 1;
}

static int
do_write_each_hit_as_separate_molecule(Molecule& m,
                                       int nhits,  // may be less than in SRESULTS
                                       const Substructure_Results& sresults,
                                       int* remove_atom, int* hits_in_fragment,
                                       Molecule_Output_Object& output) {
  int matoms = m.natoms();

  for (int i = 0; i < nhits; i++) {
    Set_of_Atoms e = *sresults.embedding(i);
    if (! OkSpatialConstraints(m, e, spatial_position)) {
      continue;
    }

    Molecule mcopy(m);

    mcopy.set_name(m.name());

    if (invert_subset) {
      set_vector(remove_atom, matoms, 0);
    } else {
      set_vector(remove_atom, matoms, 1);
    }

    identify_atoms_to_be_removed(e, remove_atom);
    identify_fragments_with_matched_atoms(m, e, hits_in_fragment);

    if (include_atoms_within > 0.0f) {
      identify_atoms_proximal_to_retained_atoms(m, remove_atom, include_atoms_within);
    }

    do_output(mcopy, remove_atom, hits_in_fragment, output);
  }

  return output.good();
}

static int
molecule_subset(Molecule& m, int* hits_in_fragment, Molecule_Output_Object& output) {
  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  int nq = queries.number_elements();

  int queries_matching = 0;

  int matoms = m.natoms();

  int* remove_atom = new int[matoms];
  std::unique_ptr<int[]> free_remove_atom(remove_atom);

  if (write_each_hit_as_separate_molecule) {
    ;
  } else if (invert_subset) {
    set_vector(remove_atom, matoms, 0);
  } else {
    set_vector(remove_atom, matoms, 1);
  }

  Molecule_to_Match target(&m);

  for (Substructure_Hit_Statistics* qry : queries) {
    Substructure_Results sresults;

    int nhits = qry->substructure_search(target, sresults);

    if (verbose > 2) {
      cerr << nhits << " matches to query " << qry->comment() << '\n';
    }

    if (0 == nhits) {
      if (verbose > 1) {
        cerr << "Zero hits to query '" << qry->comment()
             << "', only matched " << sresults.max_query_atoms_matched_in_search()
             << " query atoms\n";
      }

      continue;
    }

    queries_matching++;

    if (max_matches_to_process > 0 && nhits > max_matches_to_process) {
      nhits = max_matches_to_process;
    }

    if (write_each_hit_as_separate_molecule) {
      Molecule mcopy(m);
      mcopy.set_name(m.name());
      do_write_each_hit_as_separate_molecule(mcopy, nhits, sresults, remove_atom,
                                             hits_in_fragment, output);
      continue;
    }

    // This seems wrong, all hits get processed in the loop below.
    // TODO:ianwatson investigate
    if (embedding_to_process < 0) {
      ;
    } else if (embedding_to_process < nhits) {
      const Set_of_Atoms* e = sresults.embedding(embedding_to_process);
      identify_atoms_to_be_removed(*e, remove_atom);
      identify_fragments_with_matched_atoms(m, *e, hits_in_fragment);
      continue;
    } else {
      cerr << "Requested embedding " << embedding_to_process << " but only " << nhits
           << " matches\n";
      return 0;
    }

    for (int i = 0; i < nhits; i++) {
      Set_of_Atoms e = *sresults.embedding(i);
      if (! OkSpatialConstraints(m, e, spatial_position)) {
        continue;
      }
      if (e.empty()) {
        cerr << "Embedding empty\n";
      }

      identify_atoms_to_be_removed(e, remove_atom);
      identify_fragments_with_matched_atoms(m, e, hits_in_fragment);
    }
  }

  if (verbose > 1) {
    cerr << queries_matching << " of " << nq << " queries hit '" << m.name() << "'\n";
  }

  if (0 == queries_matching) {
    return handle_no_queries_matching(m, output);
  }

  if (verbose > 1) {
    cerr << queries_matching << " of " << nq << " matched '" << m.name() << "'\n";
  }

  if (write_each_hit_as_separate_molecule) {  //  output already done
    return output.good();
  }

  // If all atoms are removed - 3D constraints not matched, return now.
  if (std::count(remove_atom, remove_atom + matoms, 1) == matoms) {
    return 1;
  }

  return do_output(m, remove_atom, hits_in_fragment, output);
}

static int
molecule_subset(data_source_and_type<Molecule>& input, Molecule_Output_Object& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    int* tmp;

    if (!also_write_small_fragments) {
      tmp = nullptr;
    } else if (1 == m->number_fragments()) {
      tmp = nullptr;
    } else {
      tmp = new_int(m->number_fragments());
    }

    std::unique_ptr<int[]> free_tmp(tmp);

    if (create_subsets_by_bond) {
      if (!molecule_subset_by_bond(*m, tmp, output)) {
        return 0;
      }
    } else if (!molecule_subset(*m, tmp, output)) {
      return 0;
    }
  }

  return 1;
}

static int
molecule_subset(const char* fname, FileType input_type, Molecule_Output_Object& output) {
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  return molecule_subset(input, output);
}

static int
display_standard_dash_M_options(std::ostream& os) {
  os << "  -M aeuao       Atom Environments match Unmatched Atoms Only\n";
  os << "  -M imp2exp     make implicit Hydrogens explicit\n";
  os << "  -M onlysubiso  isotopic atoms in query molecules indicate subsitution "
        "points\n";
  os << "  -M iso=nn      apply isotope <nn> to the points of breakage\n";

  return os.good();
}

static int
molecule_subset(int argc, char** argv) {
  Command_Line cl(argc, argv, "A:E:S:q:s:o:i:nm:z:vbM:f:T:C:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  if (!cl.option_present('q') && !cl.option_present('s')) {
    cerr << "Must specify one or more queries via the -q or -s option(s)\n";
    usage(8);
  }

  if (cl.option_present('M')) {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++)) {
      if ("aeuao" == m) {
        set_atom_environment_only_matches_unmatched_atoms(1);
      } else if ("imp2exp" == m) {
        make_implicit_hydrogens_explicit = 1;
      } else if ("onlysubiso" == m) {
        set_substituents_only_at_isotopic_atoms(1);
        if (verbose) {
          cerr << "When building queries from molecules, isotopic atoms indicate "
                  "substitution points\n";
        }
      } else if (m.starts_with("iso=")) {
        m.remove_leading_chars(4);
        if (!m.numeric_value(isotope_at_break_points) || isotope_at_break_points <= 0) {
          cerr << "The isotope at break points must be a whole +ve number\n";
          usage(3);
        }

        if (verbose) {
          cerr << "Will label break points with isotope " << isotope_at_break_points
               << '\n';
        }
      } else if (m == "coords") {
        set_append_coordinates_after_each_atom(1);
        if (verbose) {
          cerr << "Will append coordinates after each atom in smiles\n";
        }
      } else if ("help" == m) {
        display_standard_dash_M_options(cerr);
        return 0;
      } else {
        cerr << "Unrecognised -m qualifier '" << m << "'\n";
        usage(4);
      }
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose, 'q')) {
      cerr << "Cannot read queries (-q)\n";
      return 4;
    }
  }

  if (cl.option_present('s')) {
    int i = 0;
    const_IWSubstring s;
    while (cl.value('s', s, i++)) {
      Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;
      if (!q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 3;
      }

      queries.add(q);
    }
  }

  if (verbose) {
    cerr << "Created " << queries.number_elements() << " queries\n";
  }

  if (cl.option_present('T')) {
    if (!cl.value('T', include_atoms_within) || include_atoms_within <= 0.0f) {
      cerr << "The include atoms distance (-T) option must be a valid distance\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will include atoms within a spatial distance of " << include_atoms_within
           << " angstroms\n";
    }
  }

  if (cl.option_present('C')) {
    const_IWSubstring token;
    for (int i = 0; cl.value('C', token, i); ++i) {
      std::unique_ptr<SpatialPosition> pos = std::make_unique<SpatialPosition>();
      if (! pos->Build(token)) {
        cerr << "Invalid spatial position specification '-C " << token << "'\n";
        return 1;
      }
      spatial_position << pos.release();
    }
    if (verbose) {
      cerr << "Defined " << spatial_position.size() << " spatial position selections\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  Molecule_Output_Object output;

  if (!cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (!output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    usage(28);
  }

  if (cl.option_present('S')) {
    const_IWSubstring stem_for_output = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, stem_for_output)) {
      cerr << "Cannot overwrite input file(s), '" << stem_for_output << "'\n";
      return 4;
    }

    if (!output.new_stem(stem_for_output)) {
      cerr << "Cannot initialise output stem '" << stem_for_output << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Will use '" << stem_for_output << "' for output file(s)\n";
    }
  } else {
    output.new_stem("-");
  }

  if (cl.option_present('n')) {
    invert_subset = 1;
    if (verbose) {
      cerr << "Will write atoms NOT matching the query(s)\n";
    }
  }

  if (cl.option_present('f')) {
    also_write_small_fragments = 1;

    if (verbose) {
      cerr << "Small fragments will also be included\n";
    }
  }

  if (also_write_small_fragments && invert_subset) {
    cerr << "Sorry, the -n and -f options are mutually exclusive\n";
    usage(4);
  }

  if (cl.option_present('b')) {
    create_subsets_by_bond = 1;

    if (verbose) {
      cerr << "Subsets based on bonds will be used\n";
    }
  }

  if (create_subsets_by_bond && also_write_small_fragments) {
    cerr << "Sorry the -f and -b options are mutually exclusive\n";
    usage(4);
  }

  if (cl.option_present('m')) {
    const_IWSubstring m = cl.string_value('m');

    if ("each" == m) {
      write_each_hit_as_separate_molecule = 1;

      if (verbose) {
        cerr << "Each substructure match written as a separate molecule\n";
      }
    } else if (m.starts_with("do=")) {
      m.remove_leading_chars(3);

      if (!m.numeric_value(embedding_to_process) || embedding_to_process < 0) {
        cerr << "The embedding to process must be a whole non-negative number '" << m
             << "'\n";
        return 8;
      }

      if (verbose) {
        cerr << "Will process embedding " << embedding_to_process << '\n';
      }
    } else if (!m.numeric_value(max_matches_to_process) || max_matches_to_process < 1) {
      cerr << "The number of substructure hits to process must be a whole positive "
              "number\n";
      usage(5);
    } else if (verbose) {
      cerr << "Will process a maximum of " << max_matches_to_process
           << " embeddings per query match\n";
    }
  }

  if (max_matches_to_process > 0 && embedding_to_process >= max_matches_to_process) {
    cerr << "Fatal error. max_matches_to_process = " << max_matches_to_process << '\n';
    cerr << "But you asked to process embedding number " << embedding_to_process << '\n';
    cerr << "impossible\n";
    usage(37);
  }

  if (cl.option_present('z')) {
    int i = 0;
    const_IWSubstring z;
    while (cl.value('z', z, i++)) {
      if ('i' == z) {
        ignore_queries_not_hitting = 1;
        if (verbose) {
          cerr << "Will ignore any query which does not match\n";
        }
      } else if ('w' == z) {
        write_molecules_not_matching_queries = 1;
        ignore_queries_not_hitting = 1;
        if (verbose) {
          cerr << "Will write molecules not matched by any query\n";
        }
      } else if ('b' == z) {
        blank_line_to_cout_for_non_matches = 1;
        ignore_queries_not_hitting = 1;

        if (verbose) {
          cerr << "Blank line to cout for non matching molecules\n";
        }
      } else if (z.starts_with("nmapp=")) {
        append_to_unchanged_molecules = z;
        append_to_unchanged_molecules.remove_leading_chars(6);

        if (' ' != append_to_unchanged_molecules[0]) {
          append_to_unchanged_molecules.insert_before(0, ' ');
        }

        ignore_queries_not_hitting = 1;
        write_molecules_not_matching_queries = 1;
        if (verbose) {
          cerr << "Will append '" << append_to_unchanged_molecules
               << "' to the names of non-subsetted molecules\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(7);
      }
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!molecule_subset(cl[i], input_type, output)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < atoms_removed.number_elements(); i++) {
      if (atoms_removed[i]) {
        cerr << atoms_removed[i] << " molecules had " << i << " atoms removed\n";
      }
    }

    for (int i = 0; i < atoms_in_molecules_written.number_elements(); i++) {
      if (atoms_in_molecules_written[i]) {
        cerr << atoms_in_molecules_written[i] << " final molecules had " << i
             << " atoms\n";
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = molecule_subset(argc, argv);

  return rc;
}

/* arch-tag: 0c5a495c-cc67-44a7-ba4a-07c6eb29add2

*/
