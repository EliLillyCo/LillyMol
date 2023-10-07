/*
  Rearrange a smiles to make a matched atom the first
*/

#include <stdlib.h>

#include <iostream>
#include <memory>
#include <vector>

#include "absl/algorithm/container.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_written = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

/*

*/

static int query_atom_number_to_match = 0;

static int ignore_molecules_not_hitting = 0;

static int write_molecules_not_hitting = 0;

static int take_first_of_multiple_hits = 0;

static int skip_queries_with_multiple_hits = 0;

static int do_each_of_multiple_hits = 0;

// Because we write the output as text, we can add a prefix

static IWString prefix;

static Chemical_Standardisation chemical_standardisation;

// We want to be able to write a long chain molecule starting
// at one end of the chain.
static int first_atom_is_highest_eccentricity = 0;

/*
  Useless option (seemingly). Just used for testing.
  Write a molecular subset - just those atoms matched
*/

static int write_subset = 0;

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
  cerr << "Creates smiles starting with an atom hit by a query\n";
  cerr << "  -q <query>     specify query or queries\n";
  cerr << "  -s <smarts>    specify one or more smarts queries\n";
  cerr << "  -m <number>    query atom number to use (default 0)\n";
  cerr << "  -z i           ignore molecules not hitting any queries\n";
  cerr << "  -z w           write molecules that don't hit any queries\n";
  cerr << "  -z f           take the first of multiple hits\n";
  cerr << "  -z m           skip queries with multiple hits\n";
  cerr << "  -z each        when a query matches multiple times, do each match\n";
  cerr << "  -k             do NOT perceive symmetry related matches\n";
  cerr << "  -u             unique embeddings only\n";
  cerr << "  -P <prefix>    prepend <  prefix> to all smiles (useful for adding special atoms)\n";
  cerr << "  -e             first atom has longest distance to another atom - good for long chains\n";
  cerr << "  -b             write molecular subset (matched atoms only)\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  (void)display_standard_chemical_standardisation_options(cerr, 'g');
  (void)display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
firstatom(Molecule& m, atom_number_t a, const int* process_these_atoms,
          std::ostream& output) {
  Smiles_Information sminfo;

  IWString smiles = m.smiles_starting_with_atom(a, sminfo, process_these_atoms);

  if (prefix.length()) {
    output << prefix;
  }

  output << smiles << ' ' << m.name() << '\n';

  molecules_written++;

  return output.good();
}

// For now, we just compute the longest distance to another atom.
int
Eccentricity(Molecule& m,
             atom_number_t zatom) {
  const int matoms = m.natoms();

  const int zfrag = m.fragment_membership(zatom);

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (i == zatom) {
      continue;
    }
    if (m.fragment_membership(i) != zfrag) {
      continue;
    }

    if (int d = m.bonds_between(zatom, i); d > rc) {
      rc = d;
    }
  }

  return rc;
}

static int
FirstAtomIsHighestEccentricity(Molecule& m, std::ostream& output) {
  Molecule_to_Match target(&m);

  struct AtomEccentricity {
    atom_number_t atom;
    int eccentricity;
  };

  m.recompute_distance_matrix();

  std::vector<AtomEccentricity> candidates;
  candidates.reserve(m.natoms());

  for (Substructure_Hit_Statistics* q : queries) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }
    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      atom_number_t a = (*embedding)[query_atom_number_to_match];
      int eccentricity = Eccentricity(m, a);
      candidates.push_back(AtomEccentricity(a, eccentricity));
    }
  }

  absl::c_sort(candidates, [](const AtomEccentricity& s1, const AtomEccentricity& s2) {
    return s1.eccentricity > s2.eccentricity;
  });

#ifdef DEBUG_CANDIDATE_SORT
  cerr << "Sorging complete\n";
  for (const auto& x : candidates) {
    cerr << "atom " << x.atom << " value " << x.eccentricity << '\n';
  }
#endif

  return firstatom(m, candidates[0].atom, nullptr, output);
}

static int
firstatom(Molecule& m, const Set_of_Atoms* embedding, std::ostream& output) {
  if (!embedding->ok_index(query_atom_number_to_match)) {
    cerr << "Yipes, you asked for matched atom " << query_atom_number_to_match
         << " but embedding contains only " << embedding->number_elements() << " atoms\n";
    return 0;
  }

  atom_number_t a = embedding->item(query_atom_number_to_match);

  if (nullptr != m.chiral_centre_at_atom(a)) {
    if (verbose > 1) {
      cerr << "Cannot start smiles with chiral centre, skipping...\n";
    }
    return 1;
  }

  int* process_these_atoms;
  if (write_subset) {
    process_these_atoms = new_int(m.natoms());
    embedding->set_vector(process_these_atoms, 1);
  } else {
    process_these_atoms = nullptr;
  }

  int rc = firstatom(m, a, process_these_atoms, output);

  if (nullptr != process_these_atoms) {
    delete process_these_atoms;
  }

  return rc;
}

static int
firstatom(Molecule& m, std::ostream& output) {
  if (first_atom_is_highest_eccentricity) {
    return FirstAtomIsHighestEccentricity(m, output);
  }

  Molecule_to_Match target(&m);

  int nq = queries.number_elements();
  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (verbose > 1) {
      cerr << nhits << " hits to query " << i << '\n';
    }

    if (1 == nhits) {  // hopefully the most common case
      return firstatom(m, sresults.embedding(0), output);
    }

    if (0 == nhits) {  // maybe one of the other queries will hit
      continue;
    }

    //  must have > 1 hits

    if (take_first_of_multiple_hits) {
      return firstatom(m, sresults.embedding(0), output);
    }

    if (skip_queries_with_multiple_hits) {
      continue;
    }

    if (do_each_of_multiple_hits) {
      for (int j = 0; j < nhits; j++) {
        if (!firstatom(m, sresults.embedding(j), output)) {
          return 0;
        }
      }

      return 1;
    }

    cerr << nhits << " hits to query " << i << " Skipping to next query\n";
  }

  if (verbose) {
    cerr << "None of " << nq << " queries hit '" << m.name()
         << "' or multiple matches not handled\n";
  }

  // If we come to here, none of the queries hit

  if (!ignore_molecules_not_hitting) {
    cerr << "No queries hit '" << m.name() << "' fatal\n";
    return 0;
  }

  if (write_molecules_not_hitting) {
    output << m.smiles() << ' ' << m.name() << '\n';
    return output.good();
  }

  return 1;
}

static int
firstatom(data_source_and_type<Molecule>& input, std::ostream& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1) {
      cerr << molecules_read << " read '" << m->name() << "'\n";
    }

    if (chemical_standardisation.active()) {
      chemical_standardisation.process(*m);
    }

    if (!firstatom(*m, output)) {
      return 0;
    }
  }

  return 1;
}

static int
firstatom(const char* fname, FileType input_type, std::ostream& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 1;
  }

  return firstatom(input, output);
}

static int
firstatom(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:q:s:i:z:ubP:m:g:ke");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    cerr << "Cannot process -E option\n";
    usage(11);
  }

  if (!process_standard_aromaticity_options(cl)) {
    cerr << "Cannot process aromaticity -A options\n";
    usage(13);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      cerr << "Cannot parse chemical sstandardisation specification (-g option)\n";
      return 17;
    }
  }

  if (!cl.option_present('q') && !cl.option_present('s')) {
    cerr << "Must specify queries via the -q or -s options\n";
    usage(2);
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose, 'q')) {
      cerr << "Cannot process queries (-q option)\n";
      return 4;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    int i = 0;
    while (cl.value('s', s, i++)) {
      Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;
      if (!q->create_from_smarts(s)) {
        cerr << "Cannot process smarts '" << s << "'\n";
        return 12;
      }

      queries.add(q);
    }
  }

  if (verbose) {
    cerr << "Read " << queries.number_elements() << " queries\n";
  }

  if (cl.option_present('k')) {
    for (int i = 0; i < queries.number_elements(); ++i) {
      queries[i]->set_do_not_perceive_symmetry_equivalent_matches(1);
    }

    if (verbose) {
      cerr << "Will not perceive symmetry related matches\n";
    }
  }

  if (cl.option_present('e')) {
    first_atom_is_highest_eccentricity = 1;
    if (verbose) {
      cerr << "The first atom will be the one with longest distance to another atom\n";
    }
  }

  if (cl.option_present('P')) {
    cl.value('P', prefix);
    if (verbose) {
      cerr << "'" << prefix << "' will be prepended to all output\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', query_atom_number_to_match) || query_atom_number_to_match < 0) {
      cerr << "The query atom number (-m option) must be a non-negative number\n";
      usage(21);
    }

    if (verbose) {
      cerr << "Will use matched query atom " << query_atom_number_to_match << '\n';
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
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

  if (cl.option_present('z')) {
    const_IWSubstring z;
    int i = 0;
    while (cl.value('z', z, i++)) {
      if ('i' == z) {
        ignore_molecules_not_hitting = 1;
        if (verbose) {
          cerr << "Will ignore molecules not hitting any queries\n";
        }
      } else if ('w' == z) {
        write_molecules_not_hitting = 1;
        if (verbose) {
          cerr << "Will write molecules not hitting any queries\n";
        }
      } else if ('f' == z) {
        take_first_of_multiple_hits = 1;
        if (verbose) {
          cerr << "Will process the first of any multiple hits\n";
        }
      } else if ('m' == z) {
        skip_queries_with_multiple_hits = 1;
        if (verbose) {
          cerr << "Will skip any queries with multiple hits\n";
        }
      } else if ("each" == z) {
        do_each_of_multiple_hits = 1;
        if (verbose) {
          cerr << "Will process each instance of multiple hits\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(19);
      }
    }

    if (take_first_of_multiple_hits && skip_queries_with_multiple_hits) {
      cerr << "The combinations '-z f' and '-z m' are mutually exclusive\n";
      usage(21);
    }
  }

  if (cl.option_present('u')) {
    int nq = queries.number_elements();
    for (int i = 0; i < nq; i++) {
      queries[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('b')) {
    write_subset = 1;
    if (verbose) {
      cerr << "Only the matched atoms will be written - why do you want this?\n";
    }
  }

  std::ofstream output;

  std::ostream* optr;  // will point to output or cout

  if (cl.option_present('S')) {
    output.open(cl.option_value('S'), std::ios::out);

    if (!output.good()) {
      cerr << "Cannot open '" << cl.string_value('S') << "'\n";
      return 17;
    }

    if (verbose) {
      cerr << "Output written to '" << cl.string_value('S') << "'\n";
    }

    optr = &output;
  } else {
    optr = &std::cout;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!firstatom(cl[i], input_type, *optr)) {
      cerr << "Cannot process input file '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    int nq = queries.number_elements();
    for (int i = 0; i < nq; i++) {
      cerr << "Report on query " << i << '\n';
      queries[i]->report(cerr, verbose);
    }
    cerr << molecules_written << " molecules written\n";
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = firstatom(argc, argv);

  return rc;
}
