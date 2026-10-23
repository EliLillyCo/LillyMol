/* We want to fingerprint just a substructure - useful for
  grouping molecules
*/

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "circular_fingerprint.h"
#include "iwecfp_lib.h"

using std::cerr;

// When dealing with queries that are down the bond types, we need to keep
// track of the matched atom numbers that define the bond.
struct PairOfMatchedAtoms {
  uint32_t a1;
  uint32_t a2;

  PairOfMatchedAtoms() {
    a1 = kInvalidAtomNumber;
    a2 = kInvalidAtomNumber;
  }

  PairOfMatchedAtoms(atom_number_t c1, atom_number_t c2) : a1(c1), a2(c2) {
  }
};

class Options {
 private:
  const char* prog_name = nullptr;

  int verbose = 0;

  Chemical_Standardisation chemical_standardisation;
  Element_Transformations element_transformations;
  int reduce_to_largest_fragment = 0;

  resizable_array_p<Substructure_Hit_Statistics> queries;

  IWString fingerprint_tag;

  IWString empty_fingerprint_with_closing_angle_bracket_and_newline;

  uint64_t molecules_processed = 0;

  int break_at_first_match = 0;

  int fingerprint_each_substructure_match = 0;

  Molecule_Output_Object stream_for_labelled_molecules;

  int extend_shell = 0;

  std::vector<std::optional<PairOfMatchedAtoms>> query_is_down_the_bond;

  IWString input_smiles_tag = "$SMI<";
  IWString output_smiles_tag = "SBSMI<";
  IWString tag_for_isotopically_labelled_parent;

  // when $SMI< is written we can mark the included atoms.
  int write_smiles_with_included_atoms_marked = 0;

  IWString identifier_tag = "PCN<";

  // True if the -f option is specified.
  int working_as_filter = 0;

  IWString natoms_tag;

  int produce_atom_pair_fingerprint = 0;

  int ec_fingerprint_radius = -1;
  iwecfp::Iwecfp ec_fingerprint_generator;

  // By default we stop processing.
  int ignore_molecules_not_matching_any_queries = 0;

  uint32_t molecules_not_matching_queries = 0;

  int join_disconnected_sections = 0;

  Atom_Typing_Specification atom_typing_specification;

  extending_resizable_array<int> disconnected_groups;

  int truncate_counted_fingerprints_to_one = 0;

  Circular_Fingerprint_Generator circular_fingerprint_generator;

  // Private functions.
  int write_sparse_fingerprint(Sparse_Fingerprint_Creator& sfc,
                               IWString_and_File_Descriptor& output);
  int write_natoms_if_needed(int natoms, IWString_and_File_Descriptor& output);
  int do_produce_atom_pair_fingerprint(Molecule& m, const int* include_these_atoms,
                                       const uint32_t* atype,
                                       IWString_and_File_Descriptor& output);
  int do_produce_ec_fingerprint(Molecule& m, const int* include_these_atoms,
                                const uint32_t* atype,
                                IWString_and_File_Descriptor& output);
  void WriteSmiles(Molecule& m, const int* include_these_atoms,
                   IWString_and_File_Descriptor& output);
  int WriteNoQueryMatches(Molecule& m,
                          IWString_and_File_Descriptor& output);
  int fingerprint_substructure2(Molecule& m, const int* include_these_atoms,
                                const uint32_t* atype,
                                IWString_and_File_Descriptor& output);
  int do_extend_shell(Molecule& m, int* include_these_atoms, int depth);
  int join_disconnected_fragments_2(Molecule& m, atom_number_t a1, atom_number_t astop,
                                    int* include_these_atoms, int flag);
  int join_disconnected_fragments(Molecule& m, const Set_of_Atoms& a1,
                                  const Set_of_Atoms& a2, int* include_these_atoms,
                                  int flag);
  int do_join_disconnected_sections(Molecule& m, int* include_these_atoms);
  int do_fingerprint_each_substructure_match(Molecule& m,
                                             IWString_and_File_Descriptor& output);
  int WriteLabelledForm(Molecule& m, const int* include_these_atoms,
                        Molecule_Output_Object& output);
  int fingerprint_substructure(Molecule& m, int* include_these_atoms,
                               int embeddings_processed,
                               IWString_and_File_Descriptor& output);
  void preprocess(Molecule& m);
  int fingerprint_substructure(Molecule& m, IWString_and_File_Descriptor& output);
  int fingerprint_substructure_filter(const const_IWSubstring& smiles,
                                      IWString_and_File_Descriptor& output);
  int fingerprint_substructure_filter(iwstring_data_source& input,
                                      IWString_and_File_Descriptor& output);
  int fingerprint_substructure_filter(IWString_and_File_Descriptor& output);
  int fingerprint_substructure(data_source_and_type<Molecule>& input,
                               IWString_and_File_Descriptor& output);
  int fingerprint_substructure(const char* fname, FileType input_type,
                               IWString_and_File_Descriptor& output);

 public:
  int Main(int argc, char** argv);
};

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
  cerr << R"(Fingerprints just a subset of the atoms in a molecule.
 -q ...        query specifications for identifying the subset.
 -s <smarts>   smarts to identify the subset.
 -x <number>   include atoms within <number> atoms of the subset.
 -d a1-a2      query is a down-the-bond type. Fingerprint all atoms down the bond
               defined by matched atoms a1->a2.
 -P ...        atom typing specification, enter '-P help' for info.
 -I <stem>     file name stem for isotopically labelled subset molecules.
 -Y ...        standard fingerprint options.
 -J <tag>      tag to use for fingerprints.
 -M            produce atom pair fingerprints.
 -C <radius>    produce EC fingerprints with maximum shell radius <radius>.
 -f            work as a filter.
 -T ...        standard element transformations -T I=Cl -T Br=Cl ....
 -O ...        Other options, enter '-O help' for info.
 -l            reduce to largest fragment.
 -i <type>     input specification.
 -g ...        chemical standardisation options.
 -E ...        standard element specifications.
 -A ...        standard aromaticity specifications.
 -v            verbose output.

To fingerprint all atoms within 4 bonds of an imidazole
fingerprint_substructure -s 'n1cncc1' -x 4 file.smi > file.gfp
or for a non-fused imidazole
fingerprint_substructure -s '[IWfss1n]1cncc1' -x 4 file.smi > file.gfp
)";
  // clang-format on

  ::exit(rc);
}

static void
DisplayDashOOptions(int rc) {
  cerr << R"(The following options are recognised.
 -O trunc01             truncate any counted fingerprints to 0/1.
 -O NATAG=<tag>         write the number of atoms in the subset the output
                         -O NATAG=NATOMS results in NATOMS<4> in the output.
 -O join=<n>            if there are disconnected sets of matched atoms closer than <n> bonds
                         apart, also include all the atoms between them in the fingerprint.
 -O INTAG=<tag>         smiles tag when reading GFP fingerprints.
 -O OUTTAG=<tag>        smiles tag when writing GFP fingerprints.
 -O ISO=<tag>           tag for isotopically labelled parent smiles in the output.
 -O isosmi              write the molecule's smiles with isotopes - only works when reading a structure
                         file, NOT when reading a pipeline, the -f option.
)";

  ::exit(rc);
}

static void
DisplayDashzOptions(int rc) {
  cerr << R"(The following options related to substructure matching are recognised.
 -z ...        options related to substructure matching. Enter -'z help' for info.
 -z i          ignore molecules not hitting any queries.
 -z e          fingerprint each substructure match. Generates multiple fingerprints.
 -z f          if multiple matches, take the first.
 -z d
)";

  ::exit(rc);
}

// Parse `d` as (\d+)-(\d+) and place the result in `destination`.
static int
AddQueryIsDownTheBond(const IWString& d,
                      std::vector<std::optional<PairOfMatchedAtoms>>& destination) {
  const_IWSubstring a1, a2;
  if (!d.split(a1, '-', a2) || a1.empty() || a2.empty()) {
    cerr << "AddQueryIsDownTheBond:invalid form, must be (\\d+)-(\\d+) '" << d << "'\n";
    return 0;
  }

  uint32_t i1, i2;
  if (!a1.numeric_value(i1) || !a2.numeric_value(i2) || i1 == i2) {
    cerr << "AddQueryIsDownTheBond:invalid numerics '" << d << "'\n";
    return 0;
  }

  destination.emplace_back(PairOfMatchedAtoms(i1, i2));

  return 1;
}

int
Options::write_sparse_fingerprint(Sparse_Fingerprint_Creator& sfc,
                                  IWString_and_File_Descriptor& output) {
  if (truncate_counted_fingerprints_to_one) {
    sfc.flatten_to_01();
  }

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(fingerprint_tag, tmp);

  output << tmp << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::write_natoms_if_needed(int natoms, IWString_and_File_Descriptor& output) {
  if (natoms_tag.length()) {
    output << natoms_tag << natoms << ">\n";
  }

  return output.good();
}

int
Options::WriteNoQueryMatches(Molecule& m,
                             IWString_and_File_Descriptor& output) {
  ++molecules_not_matching_queries;

  if (! ignore_molecules_not_matching_any_queries) {
    cerr << m.smiles() << " no query matches\n";
    return 0;
  }

  // Ignoring the non-match, write an empty fingerprint.

  if (! working_as_filter) {
    output << input_smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  output << empty_fingerprint_with_closing_angle_bracket_and_newline;
  write_natoms_if_needed(0, output);

  if (!working_as_filter) {
    output << "|\n";
  }

  return 1;
}

/*
  Same code from map.cc, although I don't think there is any
  reason it needs to be the same
*/

static uint32_t
form_atom_pair(uint32_t atype1, int distance, uint32_t atype2) {
  if (atype1 > atype2) {
    std::swap(atype1, atype2);
  }

  return 975535 * atype1 + distance * 7927 + atype2;
}

int
Options::do_produce_atom_pair_fingerprint(Molecule& m, const int* include_these_atoms,
                                          const uint32_t* atype,
                                          IWString_and_File_Descriptor& output) {
  int matoms = m.natoms();

  Sparse_Fingerprint_Creator sfp;

  for (int i = 0; i < matoms; i++) {
    if (!include_these_atoms[i]) {
      continue;
    }

    for (int j = i + 1; j < matoms; j++) {
      if (!include_these_atoms[j]) {
        continue;
      }

      int d = m.bonds_between(i, j);

      uint32_t b = form_atom_pair(atype[i], atype[j], d);

      sfp.hit_bit(b);
    }
  }

  if (!write_sparse_fingerprint(sfp, output)) {
    return 0;
  }

  const int atoms_in_subset = std::count_if(
      include_these_atoms, include_these_atoms + matoms, [](int s) { return s != 0; });
  return write_natoms_if_needed(atoms_in_subset, output);
}


int
Options::do_produce_ec_fingerprint(Molecule& m, const int* include_these_atoms,
                                   const uint32_t* atype,
                                   IWString_and_File_Descriptor& output) {
  if (!working_as_filter) {
    WriteSmiles(m, include_these_atoms, output);
  }

  Sparse_Fingerprint_Creator sfc;
  const iwecfp::FingerprintResult result =
      ec_fingerprint_generator.Fingerprint(m, atype, include_these_atoms, &sfc);
  if (result == iwecfp::FingerprintResult::kFatal) {
    return 0;
  }

  if (!write_sparse_fingerprint(sfc, output)) {
    return 0;
  }

  const int atoms_in_subset = std::count_if(
      include_these_atoms, include_these_atoms + m.natoms(), [](int s) { return s != 0; });
  write_natoms_if_needed(atoms_in_subset, output);

  if (!working_as_filter) {
    output << "|\n";
  }

  return output.good();
}

void
Options::WriteSmiles(Molecule& m, const int* include_these_atoms,
                     IWString_and_File_Descriptor& output) {
  if (!write_smiles_with_included_atoms_marked) {
    output << input_smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
    return;
  }

  std::unique_ptr<isotope_t[]> isotopes = m.GetIsotopes();
  m.transform_to_non_isotopic_form();
  m.set_isotopes(include_these_atoms);

  output << input_smiles_tag << m.smiles() << ">\n";
  output << identifier_tag << m.name() << ">\n";

  m.set_isotopes(isotopes.get());

  output.write_if_buffer_holds_more_than(4096);

  return;
}

int
Options::fingerprint_substructure2(Molecule& m, const int* include_these_atoms,
                                   const uint32_t* atype,
                                   IWString_and_File_Descriptor& output) {
  // cerr << "SM<ILES is " << m.smiles() << " contains " << m.natoms() << " atoms " <<
  // output_smiles_tag << m.smiles() << ">\n";
  if (!working_as_filter) {
    WriteSmiles(m, include_these_atoms, output);
  }

  const int atoms_in_subset =
      std::count_if(include_these_atoms, include_these_atoms + m.natoms(),
                    [](int s) { return s != 0; });

  IWMFingerprint fp;

  fp.construct_fingerprint(m, atype, include_these_atoms);

  IWString tmp;

  fp.daylight_ascii_representation_including_nset_info(tmp);

  output << fingerprint_tag << tmp << ">\n";

  write_natoms_if_needed(atoms_in_subset, output);

  if (!working_as_filter) {
    output << "|\n";
  }

  return output.good();
}

/*
  We mark any atoms we extend with -1
*/

int
Options::do_extend_shell(Molecule& m, int* include_these_atoms, int depth) {
  int matoms = m.natoms();

  int atoms_added = 0;

  for (int i = 0; i < matoms; i++) {
    if (include_these_atoms[i] > 0)
      ;
    else if (-99 == include_these_atoms[i]) {
      continue;
    } else if (0 == include_these_atoms[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++) {
      atom_number_t k = a->other(i, j);
      if (0 != include_these_atoms[k]) {
        continue;
      }

      include_these_atoms[k] = -99;

      atoms_added++;
    }
  }

  if (0 == atoms_added) {
    return 1;
  }

  for (int i = 0; i < matoms; i++) {
    if (-99 == include_these_atoms[i]) {
      include_these_atoms[i] = -(depth + 1);
    }
  }

  if (0 == depth) {
    return 1;
  }

  return do_extend_shell(m, include_these_atoms, depth - 1);
}

int
Options::join_disconnected_fragments_2(Molecule& m, atom_number_t a1, atom_number_t astop,
                                       int* include_these_atoms, int flag) {
  int d = m.bonds_between(a1, astop);

#ifdef DEBUG_JOIN_DISCONNECTED_FRAGMENTS
  cerr << "join_disconnected_fragments, atom " << a1 << " to " << astop << " " << d
       << " bonds\n";
#endif

  const Atom* a = m.atomi(a1);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(a1, i);

    if (j == astop) {  // got to the end
      return 1;
    }

    if (d - 1 != m.bonds_between(j, astop)) {
      continue;
    }

    if (0 == include_these_atoms[j]) {
      include_these_atoms[j] = flag;
    }

    join_disconnected_fragments_2(m, j, astop, include_these_atoms, flag);
  }

  return 1;
}

int
Options::join_disconnected_fragments(Molecule& m, const Set_of_Atoms& a1,
                                     const Set_of_Atoms& a2, int* include_these_atoms,
                                     int flag) {
  int n = a1.number_elements();

  assert(n > 0);

  for (int i = 0; i < n; i++) {
    atom_number_t j = a1[i];

    for (int k = 0; k < n; k++) {
      atom_number_t l = a2[k];

      join_disconnected_fragments_2(m, j, l, include_these_atoms, flag);
    }
  }

  return 1;
}

static int
grow_group(const Molecule& m, atom_number_t zatom, int* include_these_atoms, int flag) {
  include_these_atoms[zatom] = flag;

  int rc = 1;

  for (const Bond* b : m[zatom]) {
    atom_number_t j = b->other(zatom);

    if (1 != include_these_atoms[j]) {
      continue;
    }

    rc += grow_group(m, j, include_these_atoms, flag);
  }

  return rc;
}

/*
  We have possibly disconnected sections. Join them if they are
  close enough. We use the flag -2 to denote atoms flagged this way
*/

int
Options::do_join_disconnected_sections(Molecule& m, int* include_these_atoms) {
  int f = 2;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (include_these_atoms[i] <= 0) {
      continue;
    }

    if (include_these_atoms[i] > 1) {
      continue;
    }

    grow_group(m, i, include_these_atoms, f);
    f++;
  }

  // cerr << "Found " << (f - 2) << " disconnected groups\n";
  disconnected_groups[f - 2]++;

  if (3 == f) {  // just one grouping
    return 1;
  }

  // We have possibly separated groups.  Join them by finding shortest
  // paths between the disconnected sections

  // iw_write_array (include_these_atoms, matoms, "include_these_atoms", cerr);

  for (int f1 = 2; f1 < f; f1++) {
    for (int f2 = f1 + 1; f2 < f; f2++) {
      int closest_separation = matoms + 2;
      Set_of_Atoms a1;
      Set_of_Atoms a2;

      for (int i = 0; i < matoms; i++) {
        if (f1 != include_these_atoms[i]) {
          continue;
        }

        for (int j = i + 1; j < matoms; j++) {
          if (f2 != include_these_atoms[j]) {
            continue;
          }

          int dij = m.bonds_between(i, j);
#ifdef DEBUG_JOIN_DISCONNECTED_FRAGMENTS
          cerr << dij << " bonds between " << i << " and " << j << '\n';
#endif

          if (dij < closest_separation) {
            a1.resize_keep_storage(0);
            a2.resize_keep_storage(0);
            a1.add(i);
            a2.add(j);
            closest_separation = dij;
          } else if (dij == closest_separation) {
            a1.add(i);
            a2.add(j);
          }
        }
      }

      assert(closest_separation > 0);

#ifdef DEBUG_JOIN_DISCONNECTED_FRAGMENTS
      cerr << f1 << " to " << f2 << ", closest separation " << closest_separation << '\n';
#endif

      if (closest_separation <= 1) {  // should never be 0
        continue;
      }

      if (closest_separation > join_disconnected_sections) {
        continue;
      }

      assert(a1.number_elements() == a2.number_elements());
      assert(a1.number_elements() > 0);

      join_disconnected_fragments(m, a1, a2, include_these_atoms, -2);
    }
  }

  return 1;
}

// Given an embedding, and two matched atoms in it, fill `first_atom` and
// `next_atom`.
static bool
IdentifyMatchedAtoms(const Set_of_Atoms& embedding,
                     const PairOfMatchedAtoms& matched_atoms, atom_number_t& first_atom,
                     atom_number_t& next_atom) {
  if (embedding.size() < 2) {
    cerr << "IdentifyMatchedAtoms:Not enough matched atoms in " << embedding << '\n';
    return false;
  }

  if (!embedding.ok_index(matched_atoms.a1) || !embedding.ok_index(matched_atoms.a2)) {
    cerr << "IdentifyMatchedAtoms:invalid matched atom numbers. Requested "
         << matched_atoms.a1 << " and " << matched_atoms.a2 << " but embedding size "
         << embedding.size() << '\n';
    return false;
  }

  first_atom = embedding[matched_atoms.a1];
  next_atom = embedding[matched_atoms.a2];

  return true;
}

// `embedding` is assumed to specify a down the bond request.
// Identify all the atoms implied by that and mark them in `include_these_atoms`
// with the value `inc`.
static int MarkDownTheBond(Molecule& m, const Set_of_Atoms& embedding,
                           const PairOfMatchedAtoms& matched_atoms, int* visited,
                           int* include_these_atoms, int inc);

// The embedding is assumed to be a down the bond kind of match.
// Identify all the atoms implied by that and put them in `result`.
static int
MarkDownTheBond(Molecule& m, const Set_of_Atoms& embedding,
                const PairOfMatchedAtoms& matched_atoms, int* visited,
                Set_of_Atoms& result) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> include_these_atoms = std::make_unique<int[]>(matoms);

  const int rc =
      MarkDownTheBond(m, embedding, matched_atoms, visited, include_these_atoms.get(), 1);
  if (rc == 0) {
    return 0;
  }

  for (int i = 0; i < matoms; ++i) {
    if (include_these_atoms[i]) {
      result << i;
    }
  }

  return rc;
}

// Recursively identify all the unvisited atoms attached to `zatom`.
// For the atoms discovered, mark the corresponding entry in `include_these_atoms`
// with the value `inc`.
static int
MarkDownTheBond(Molecule& m, atom_number_t zatom, int* visited, int* include_these_atoms,
                int inc) {
  visited[zatom] = 1;
  include_these_atoms[zatom] = inc;

  int rc = 1;

  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (visited[o]) {
      continue;
    }
    rc += MarkDownTheBond(m, o, visited, include_these_atoms, inc);
  }

  return rc;
}

static int
MarkDownTheBond(Molecule& m, const Set_of_Atoms& embedding,
                const PairOfMatchedAtoms& matched_atoms, int* visited,
                int* include_these_atoms, int inc) {
  atom_number_t first_atom, next_atom;
  if (!IdentifyMatchedAtoms(embedding, matched_atoms, first_atom, next_atom)) {
    return 0;
  }

  std::fill_n(visited, m.natoms(), 0);
  visited[first_atom] = 1;

  return MarkDownTheBond(m, next_atom, visited, include_these_atoms, inc);
}

int
Options::do_fingerprint_each_substructure_match(Molecule& m,
                                                IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  m.recompute_distance_matrix();

  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);
  atom_typing_specification.assign_atom_types(m, atype.get());

  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms + matoms);
  int* all_atoms_matched = tmp.get() + matoms;
  std::fill_n(all_atoms_matched, matoms, 0);

  Sparse_Fingerprint_Creator sfc;

  std::unique_ptr<int[]> dtb;
  if (!query_is_down_the_bond.empty()) {
    dtb = std::make_unique<int[]>(matoms);
  }

  Molecule_to_Match target(&m);

  int queries_matching = 0;

  const int nq = queries.number_elements();

  for (int i = 0; i < nq; ++i) {
    Substructure_Results sresults;

    const int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    if (verbose > 2) {
      cerr << m.name() << " " << nhits << " hits to query " << i << '\n';
    }

    queries_matching++;

    //  mark all atoms within extend_shell of the matched atoms with the number 2. Matched
    //  atoms get 1

    std::fill_n(tmp.get(), matoms, 0);

    for (int j = 0; j < nhits; ++j) {
      const Set_of_Atoms* e = sresults.embedding(j);

      // We need a Set_of_Atoms containing the atoms around which we do shell expansion.
      // If this is a down the bond query, we must compute it. Otherwise it is `e`.
      Set_of_Atoms matched_atoms;

      if (i < static_cast<int>(query_is_down_the_bond.size()) &&
          query_is_down_the_bond[i]) {
        if (!MarkDownTheBond(m, *e, *query_is_down_the_bond[i], dtb.get(), matched_atoms)) {
          return 0;
        }
        matched_atoms.set_vector(tmp.get(), 1);  // note value of 1
      } else {
        e->set_vector(tmp.get(), 1);  // note value of 1
        matched_atoms = *e;
      }

      for (atom_number_t matched : matched_atoms) {
        for (int k = 0; k < matoms; ++k) {
          if (0 != tmp[k]) {  // already marked
            continue;
          }

          if (m.bonds_between(matched, k) <= extend_shell) {
            tmp[k] = 2;  // note value 2
          }
        }
      }
    }

    circular_fingerprint_generator.generate_fingerprint(m, atype.get(), tmp.get(), sfc);
    // Update the all atom match array
    for (int j = 0; j < matoms; ++j) {
      if (tmp[j] > 0) {
        all_atoms_matched[j] = 1;
      }
    }
  }

  if (queries_matching == 0) {
    return WriteNoQueryMatches(m, output);
  }

  if (!working_as_filter) {
    WriteSmiles(m, all_atoms_matched, output);
  }

  if (!write_sparse_fingerprint(sfc, output)) {
    return 0;
  }

  if (!working_as_filter) {
    output << "|\n";
  }

  return 1;
}

int
Options::WriteLabelledForm(Molecule& m, const int* include_these_atoms,
                           Molecule_Output_Object& output) {
  std::unique_ptr<isotope_t[]> isosave = m.GetIsotopes();

  m.transform_to_non_isotopic_form();
  m.set_isotopes(include_these_atoms);
  stream_for_labelled_molecules.write(m);

  m.set_isotopes(isosave.get());

  return 1;
}

int
Options::fingerprint_substructure(Molecule& m, int* include_these_atoms,
                                  int embeddings_processed,
                                  IWString_and_File_Descriptor& output) {
  if (extend_shell > 0) {
    do_extend_shell(m, include_these_atoms, extend_shell - 1);
  }

  if (join_disconnected_sections && embeddings_processed > 1) {
    do_join_disconnected_sections(m, include_these_atoms);
  }

  const int matoms = m.natoms();

  // Deliberately not an isotope_t because we use negative values.
  int* isotopes = new_int(matoms);
  std::unique_ptr<int[]> free_isotopes(isotopes);

  if (join_disconnected_sections || extend_shell) {
    for (int i = 0; i < matoms; i++) {
      if (include_these_atoms[i] < 0) {
        isotopes[i] = -include_these_atoms[i];
        include_these_atoms[i] = 1;
      } else if (include_these_atoms[i] > 0) {
        isotopes[i] = include_these_atoms[i];
        include_these_atoms[i] = 1;
      }
    }
  }

  if (stream_for_labelled_molecules.active()) {
    WriteLabelledForm(m, include_these_atoms, stream_for_labelled_molecules);
  }

  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);

  atom_typing_specification.assign_atom_types(m, atype.get());

  if (produce_atom_pair_fingerprint) {
    return do_produce_atom_pair_fingerprint(m, include_these_atoms, atype.get(), output);
  }

  if (ec_fingerprint_radius >= 0) {
    return do_produce_ec_fingerprint(m, include_these_atoms, atype.get(), output);
  }

  if (tag_for_isotopically_labelled_parent.length()) {
    m.set_isotopes(isotopes);
    output << tag_for_isotopically_labelled_parent << m.smiles() << ">\n";
  }
  // cerr << "Processing " << subset.smiles() << '\n';

  return fingerprint_substructure2(m, include_these_atoms, atype.get(), output);
}

#ifdef THIS_IS_HORRIBLE
Or maybe kind of interesting.THis generates molecules that can be partially
    aromatic.Oc(cc) cc would be a possible output from here.While interesting,
    not sure it is useful.

        int *xref = new int[matoms];
std::unique_ptr<int[]> free_xref(xref);
// All the fingerprint generators should be able to fingerprint a subset of a molecule.

Molecule subset;

if (!m.create_subset(subset, include_these_atoms, 1, xref)) {
  cerr << "fingerprint_substructure::cannot create subset\n";
  return 0;
}
// cerr << "Created subset " << subset.smiles() << '\n';

for (int i = 0; i < matoms; i++) {
  const atom_number_t j = xref[i];

  if (j < 0) {
    continue;
  }

  if (!m.is_aromatic(i)) {
    continue;
  }

  subset.set_permanent_aromatic(j, 1);
}

// Now bonds

for (int i = 0; i < matoms; i++) {
  const atom_number_t j = xref[i];

  if (j < 0) {
    continue;
  }

  if (!m.is_aromatic(i)) {
    continue;
  }

  const Atom* a = m.atomi(i);

  for (const Bond* b : *a) {
    const atom_number_t l = b->other(i);

    const atom_number_t n = xref[l];

    if (n < 0) {
      continue;
    }

    if (!m.in_same_aromatic_ring(i, l)) {
      continue;
    }

    Bond* bn = const_cast<Bond*>(subset.bond_between_atoms(j, n));

    bn->set_permanent_aromatic(1);  // just one of setting the bond to perm arom, or
                                    // setting to single should be enough

    subset.set_bond_type_between_atoms(j, n, SINGLE_BOND);
  }
}

if (nullptr != atype) {
  for (int i = 0; i < matoms; ++i) {
    if (xref[i] < 0) {
      continue;
    }

    assert(xref[i] <= i);  // this only works because we have removed atoms and the new
                           // atom number should be less than where we started

    atype[xref[i]] = atype[i];
  }
}

subset.set_name(m.name());

#endif

void
Options::preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  return;
}

int
Options::fingerprint_substructure(Molecule& m, IWString_and_File_Descriptor& output) {
  output.write_if_buffer_holds_more_than(8192);

  molecules_processed++;

  m.revert_all_directional_bonds_to_non_directional();  // until I get that working
                                                        // properly

  if (fingerprint_each_substructure_match) {
    return do_fingerprint_each_substructure_match(m, output);
  }

  const int nq = queries.number_elements();

  std::unique_ptr<int[]> include_these_atoms = std::make_unique<int[]>(m.natoms());

  int queries_matching = 0;

  int embeddings_processed = 0;

  std::unique_ptr<int[]> dtb;
  if (!query_is_down_the_bond.empty()) {
    dtb = std::make_unique<int[]>(m.natoms());
  }

  int inc = 1;
  if (extend_shell) {
    inc = extend_shell + 1;
  }

  Molecule_to_Match target(&m);

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    if (verbose > 2) {
      cerr << m.name() << " " << nhits << " hits to query " << i << '\n';
    }

    queries_matching++;

    int jstop = nhits;
    if (break_at_first_match) {
      jstop = 1;
    }

    for (int j = 0; j < jstop; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);
      // cerr << "Doing embedding " << (*e) << '\n';

      embeddings_processed++;
      if (i < static_cast<int>(query_is_down_the_bond.size()) &&
          query_is_down_the_bond[i]) {
        if (!MarkDownTheBond(m, *e, *query_is_down_the_bond[i], dtb.get(),
                             include_these_atoms.get(), inc)) {
          return 0;
        }
      } else {
        e->set_vector(include_these_atoms.get(), inc);
      }
    }

    if (break_at_first_match) {
      return fingerprint_substructure(m, include_these_atoms.get(), embeddings_processed,
                                      output);
    }
  }

  if (0 == queries_matching) {
    return WriteNoQueryMatches(m, output);
  }

  return fingerprint_substructure(m, include_these_atoms.get(), embeddings_processed,
                                  output);
}

int
Options::fingerprint_substructure_filter(const const_IWSubstring& smiles,
                                         IWString_and_File_Descriptor& output) {
  Molecule m;

  if (!m.build_from_smiles(smiles)) {
    cerr << "fingerprint_substructure_filter:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return fingerprint_substructure(m, output);
}

int
Options::fingerprint_substructure_filter(iwstring_data_source& input,
                                         IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(input_smiles_tag)) {
      continue;
    }

    buffer.remove_leading_chars(input_smiles_tag.length());
    buffer.chop();

    if (!fingerprint_substructure_filter(buffer, output)) {
      cerr << "Fatal error on line " << input.lines_read() << '\n';
      return 0;
    }
  }

  return output.good();
}

int
Options::fingerprint_substructure_filter(IWString_and_File_Descriptor& output) {
  iwstring_data_source input("-");

  if (!input.good()) {
    cerr << "fingerprint_substructure_filter:cannot open stdin?\n";
    return 0;
  }

  return fingerprint_substructure_filter(input, output);
}

int
Options::fingerprint_substructure(data_source_and_type<Molecule>& input,
                                  IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    //  output << input_smiles_tag << m->smiles() << ">\n";
    //  output << identifier_tag << m->name() << ">\n";

    preprocess(*m);

    if (!fingerprint_substructure(*m, output)) {
      return 0;
    }

    //  output << "|\n";
  }

  return output.good();
}

int
Options::fingerprint_substructure(const char* fname, FileType input_type,
                                  IWString_and_File_Descriptor& output) {
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

  return fingerprint_substructure(input, output);
}

int
Options::Main(int argc, char** argv) {
  prog_name = argv[0];

  Command_Line cl(argc, argv, "vA:E:i:g:lq:s:Y:fx:z:T:I:o:P:MJ:d:O:ej:C:");

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

  // A quick scan through the -z option to see if help is requested.
  // Otherwise execution will stop if just involed as 'fingerprint_subset -z help'.
  // the full option processing happens below.
  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == "help") {
        DisplayDashzOptions(0);
      }
    }
  }

  if (cl.option_present('O')) {
    IWString o;
    for (int i = 0; cl.value('O', o, i); ++i) {
      if (o == "trunc01") {
        truncate_counted_fingerprints_to_one = 1;
        if (verbose) {
          cerr << "Will truncate counted fingerprints to 0/1\n";
        }
      } else if (o.starts_with("NATAG=")) {
        o.remove_leading_chars(6);
        natoms_tag = o;
        natoms_tag.EnsureEndsWith('<');
        if (verbose) {
          cerr << "The number of atoms in the subset written to " << natoms_tag << '\n';
        }
      } else if (o.starts_with("join=")) {
        o.remove_leading_chars(5);
        if (!o.numeric_value(join_disconnected_sections) ||
            join_disconnected_sections < 1) {
          cerr << "Invalid join disconnected regions directive '" << o
               << "' must be +ve whole number\n";
          DisplayDashOOptions(1);
        }
        if (verbose) {
          cerr << "Will join disconnected sections closer than "
               << join_disconnected_sections << " bonds apart\n";
        }
      } else if (o == "isosmi") {
        if (cl.option_present('f')) {
          cerr << "The '-O isosmi' directive cannot be used with the -f option\n";
          return 1;
        }
        write_smiles_with_included_atoms_marked = 1;
        if (verbose) {
          cerr << "Will write the subset with isotope labels\n";
        }
      } else if (o.starts_with("INTAG=")) {
        o.remove_leading_chars(6);
        input_smiles_tag = o;
        input_smiles_tag.EnsureEndsWith('<');
      } else if (o.starts_with("OUTTAG=")) {
        o.remove_leading_chars(7);
        output_smiles_tag = o;
        if (output_smiles_tag.empty()) {
          cerr << "Empty output smiles tag\n";  // die?
        } else {
          output_smiles_tag.EnsureEndsWith('<');
        }
      } else if (o.starts_with("ISO=")) {
        o.remove_leading_chars(4);
        tag_for_isotopically_labelled_parent = o;
        tag_for_isotopically_labelled_parent.EnsureEndsWith('<');
      } else if (o == "help") {
        DisplayDashOOptions(0);
      } else {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        DisplayDashOOptions(1);
      }
    }
  }

  if (cl.option_present('T')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T')) {
      cerr << "Cannot initialise element transformations (-T)\n";
      return 0;
    }
  }

  if (cl.option_present('e')) {
    truncate_counted_fingerprints_to_one = 1;

    if (verbose) {
      cerr << "Will truncate counted fingerprints to 0/1\n";
    }
  }

  if (cl.option_present('d')) {
    IWString d;
    for (int i = 0; cl.value('d', d, i); ++i) {
      if (!AddQueryIsDownTheBond(d, query_is_down_the_bond)) {
        cerr << "Invalid query is down the bond (-d) specification '" << d << "'\n";
        return 1;
      }
    }
    if (verbose) {
      cerr << "The query defines a down the bond type query\n";
    }
  }

  if (cl.option_present('x')) {
    if (!cl.value('x', extend_shell) || extend_shell < 1) {
      cerr << "The extend shell option (-x) must be a whole positive number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will extend the shell " << extend_shell << " atoms\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', fingerprint_tag);

    fingerprint_tag.EnsureEndsWith('<');
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', ec_fingerprint_radius) || ec_fingerprint_radius < 0) {
      cerr << "The EC fingerprint radius (-C) must be a whole non-negative number\n";
      return 1;
    }
    if (verbose) {
      cerr << "Will generate EC fingerprints with radius " << ec_fingerprint_radius << '\n';
    }
  }

  const int fingerprint_types_requested = (cl.option_present('Y') ? 1 : 0) +
                                          (cl.option_present('M') ? 1 : 0) +
                                          (cl.option_present('C') ? 1 : 0);
  if (fingerprint_types_requested > 1) {
    cerr << "Only one of -C, -M and -Y may be specified\n";
    return 1;
  }

  if (cl.option_present('Y')) {
    set_include_bits_for_rings(3);  // set here so it can be changed on the command line

    if (!parse_misc_fingerprint_options(cl, 'Y', verbose)) {
      cerr << "Cannot parse fingerprint options (-Y)\n";
      return 4;
    }

    if (0 == fingerprint_tag.length()) {
      fingerprint_tag = "FPSUB<";
    }
  } else if (cl.option_present('M')) {
    produce_atom_pair_fingerprint = 1;

    if (verbose) {
      cerr << "Will produce atom pair fingerprints\n";
    }

    if (fingerprint_tag.empty()) {
      fingerprint_tag = "NCAPSUB<";
    }

    output_smiles_tag = "";  // the subset molecule is never formed
  } else if (cl.option_present('C')) {
    ec_fingerprint_generator.set_max_radius(ec_fingerprint_radius);

    if (fingerprint_tag.empty()) {
      fingerprint_tag = "NCESUB<";
    }
  } else {
    set_iwmfingerprint_nbits(1024);

    if (fingerprint_tag.empty()) {
      fingerprint_tag = "FPSUB<";
    }
  }


  if (cl.option_present('P')) {
    const const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Cannot initialise atom typing specification '" << p << "'\n";
      return 2;
    }
  } else {
    atom_typing_specification.set_user_specified_type(IWATTYPE_USP_Z | IWATTYPE_USP_A);
  }

  if (!cl.option_present('q') && !cl.option_present('s')) {
    cerr << "Must specify a substructure query via the -q or -s option\n";
    usage(3);
  }

  queries.reserve(cl.option_count('q') + cl.option_count('s') + 100);

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose)) {
      cerr << prog_name << ": cannot process queries from -q option(s)\n";
      return 1;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    int i = 0;
    while (cl.value('s', smarts, i++)) {
      std::unique_ptr<Substructure_Hit_Statistics> q =
          std::make_unique<Substructure_Hit_Statistics>();
      if (!q->create_from_smarts(smarts)) {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 1;
      }

      queries.add(q.release());
    }
  }

  if (queries.empty()) {
    cerr << "No queries read, use the -q and/or -s options to specify queries\n";
    usage(4);
  }

  if (query_is_down_the_bond.size() > queries.size()) {
    cerr << "More -d specifications (" << query_is_down_the_bond.size()
         << ") than queries (" << queries.size() << ")\n";
    return 1;
  }
  query_is_down_the_bond.resize(queries.number_elements());

  for (Substructure_Hit_Statistics* q : queries) {
    q->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if ("i" == z) {
        ignore_molecules_not_matching_any_queries = 1;
        if (verbose) {
          cerr << "Will ignore molecules not hitting any queries\n";
        }
      } else if ('f' == z) {
        break_at_first_match = 1;

        if (verbose) {
          cerr << "Will only process the first of multiple matches\n";
        }
      } else if ('e' == z) {
        fingerprint_each_substructure_match = 1;

        if (verbose) {
          cerr << "Will fingerprint all substructure matches\n";
        }
      } else if (z == "help") {
        DisplayDashzOptions(0);
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        DisplayDashzOptions(1);
      }
    }

    if (fingerprint_each_substructure_match && break_at_first_match) {
      cerr << "The -z e and -z f options are mutually incompatible\n";
      return 1;
    }

    if (fingerprint_each_substructure_match && ec_fingerprint_radius >= 0) {
      cerr << "The -z e and -C options are mutually incompatible\n";
      return 1;
    }

    if (fingerprint_each_substructure_match &&
        (0 == fingerprint_tag.length() || fingerprint_tag.starts_with("FP"))) {
      fingerprint_tag = "NCSUB<";
    }

    if (fingerprint_each_substructure_match && !atom_typing_specification.active()) {
      atom_typing_specification.set_user_specified_type(IWATTYPE_USP_Y);
    }
  }

  if (produce_atom_pair_fingerprint || fingerprint_each_substructure_match ||
      ec_fingerprint_radius >= 0) {
    empty_fingerprint_with_closing_angle_bracket_and_newline << fingerprint_tag << ">\n";
  } else {
    int produce_fingerprints = iwmfingerprint_nbits();
    assert(produce_fingerprints > 0);

    IW_Bits_Base fp(produce_fingerprints);
    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    empty_fingerprint_with_closing_angle_bracket_and_newline << fingerprint_tag << tmp
                                                             << ">\n";
  }

  if (cl.option_present('j')) {
    if (!cl.value('j', join_disconnected_sections) || join_disconnected_sections < 1) {
      cerr << "The join disconnected sections option (-j) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will join disconnected sections " << join_disconnected_sections
           << " or closer\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('I')) {
    const_IWSubstring i = cl.string_value('I');

    if (!cl.option_present('o')) {
      stream_for_labelled_molecules.add_output_type(FILE_TYPE_SMI);
    } else if (!stream_for_labelled_molecules.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      return 4;
    }

    if (stream_for_labelled_molecules.would_overwrite_input_files(cl, i)) {
      cerr << "-I option '" << i << "' cannot overwrite input file(s)\n";
      return 4;
    }

    if (!stream_for_labelled_molecules.new_stem(i)) {
      cerr << "Cannot initialise stream for labelled molecules '" << i << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Labelled molecules written to '" << i << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('f')) {
    working_as_filter = 1;

    tag_for_isotopically_labelled_parent.resize(0);

    rc = fingerprint_substructure_filter(output) ? 0 : 1;
  } else {
    FileType input_type = FILE_TYPE_INVALID;

    if (cl.option_present('i')) {
      if (!process_input_type(cl, input_type)) {
        cerr << "Cannot determine input type\n";
        usage(6);
      }
    } else if (!all_files_recognised_by_suffix(cl)) {
      return 4;
    }

    for (const char* fname : cl) {
      if (!fingerprint_substructure(fname, input_type, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        rc = 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "processed " << molecules_processed << " molecules\n";
    if (molecules_not_matching_queries) {
      cerr << molecules_not_matching_queries << " did not match any queries\n";
    }

    if (join_disconnected_sections) {
      for (int i = 0; i < disconnected_groups.number_elements(); i++) {
        if (disconnected_groups[i]) {
          cerr << disconnected_groups[i] << " molecules had " << i
               << " disconnected groups\n";
        }
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  Options options;
  int rc = options.Main(argc, argv);

  return rc;
}
