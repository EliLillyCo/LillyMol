// Implementation of Positional Analogue Scanning

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "absl/container/flat_hash_set.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace internal {

using std::cerr;

class SetOfMolecules {
  private:
    resizable_array_p<Molecule> _molecule;

  public:
    int Build(IWString& fname);

    uint64_t size() const {
      return _molecule.size();
    }

    Molecule** begin() const {
      return _molecule.begin();
    }

    Molecule** end() const {
      return _molecule.end();
    }

    Molecule* operator[](int ndx) const {
      return _molecule[ndx];
    }
    Molecule* at(int ndx) const {
      return _molecule[ndx];
    }
};

int
SetOfMolecules::Build(IWString& fname) {
  data_source_and_type<Molecule> input(FileType::FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "SetOfMolecules::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    cerr << m->smiles() << ' ' << m->name() << " read\n";
    _molecule << m;
  }

  return _molecule.number_elements();
}

}  // namespace internal

// We have two largely separate sets of functionality in the same file.
// Not sure this is a great idea, but nor is a proliferation of executables.
// This functionality is just a translation of the python impelementation.
namespace from_python {

using std::cerr;
using internal::SetOfMolecules;

void
Usage() {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Positional Analogue Scanning
 -s <smarts>             Smarts specifying the sites at which atoms are to be added.
 -q <query>              Query  specifying the sites at which atoms are to be added.
 -z i                    Ignore molecules not matching the query.
 -F <fname>              File(s) containing fragment(s) to be added.
 -p                      Write the parent molecule to the output stream.
 -x <max>                Max number of products per starting molecule.
 -v                      verbose output.
)";
  // clang-format on

  ::exit(0);
}

class SimpleEnumerator {
  private:
    int _verbose;

    uint64_t _molecules_processed;
    uint64_t _molecules_written;

    resizable_array_p<Substructure_Query> _query;

    // The k value from n choose k.
    int _combo;
    // The user will supply sets of molecules. The number of sets defined
    // is k.
    resizable_array_p<SetOfMolecules> _molecules;
    // If each SetOfMolecules has just one molecule, processing is handled
    // separately.
    int _just_one_fragment_everywhere;

    // If all the fragments being added are the same, we can cut in half
    // the number of combinations to consider.
    int _fragments_identical;

    uint32_t _made_this_molecule;

    int _write_starting_molecule;

    int _ignore_molecules_not_matchng_queries;

    // We may have matches, but not enough to attach all the fragments.
    int _ignore_not_enough_matched_atoms;

    uint32_t _max_products_per_starting_molecule;

    absl::flat_hash_set<IWString> _seen;

    uint32_t _duplicates_suppressed;
    uint32_t _invalid_valence_skipped;

  // private functions;
    int Process(Molecule& m, const Set_of_Atoms& d,
                IWString_and_File_Descriptor& output);
    int ProcessNewMolecule(Molecule& m, IWString_and_File_Descriptor& output);
    int ProcessSingleFragments(Molecule& m, const Set_of_Atoms& d,
                        IWString_and_File_Descriptor& output);
    int ProcessMultipleFragments(Molecule& m, const Set_of_Atoms& matched_atoms,
                        IWString_and_File_Descriptor& output);
    int MakeMolecule(Molecule& m, const Set_of_Atoms& matched_atom,
                const std::vector<uint32_t>& frag,
                IWString_and_File_Descriptor& output);
    int AddSingleFragments(Molecule& m,
                                 const Set_of_Atoms& matched_atoms,
                                 IWString_and_File_Descriptor& output);

  public:
    SimpleEnumerator();

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

SimpleEnumerator::SimpleEnumerator() {
  _molecules_processed = 0;
  _molecules_written = 0;
  _combo = 0;
  _ignore_molecules_not_matchng_queries = 0;
  _ignore_not_enough_matched_atoms = 0;
  _write_starting_molecule = 0;
  _made_this_molecule = 0;
  _just_one_fragment_everywhere = 0;
  _max_products_per_starting_molecule = std::numeric_limits<uint32_t>::max();
  _duplicates_suppressed = 0;
  _invalid_valence_skipped = 0;
  _fragments_identical = 0;
}

int
SimpleEnumerator::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');
  if (! cl.option_present('q')) {
    cerr << "SimpleEnumerator::Initialise:must specify query(s) via the -q option\n";
    return 0;
  } else {
    if (! process_queries(cl, _query, _verbose, 'q')) {
      cerr << "SimpleEnumerator::Initialise:cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matchng_queries = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching all the queries\n";
        }
      } else if (z == "nqm") {
        _ignore_not_enough_matched_atoms = 1;
        if (_verbose) {
          cerr << "Will ignore cases where there are not enough query matches\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (! cl.option_present('F')) {
    cerr << "Must specify fragments to add via the -F option\n";
    Usage();
  }

  if (cl.option_present('F')) {
    IWString fname;
    for (int i = 0; cl.value('F', fname, i); ++i) {
      std::unique_ptr<SetOfMolecules> m = std::make_unique<SetOfMolecules>();
      if (! m->Build(fname)) {
        cerr << "SimpleEnumerator::Initialise:cannot read " << fname << '\n';
        return 0;
      }

      _molecules << m.release();
    }

    _combo = _molecules.number_elements();
  }

  _just_one_fragment_everywhere = 1;
  for (const SetOfMolecules* mols : _molecules) {
    if (mols->size() > 1) {
      _just_one_fragment_everywhere = 0;
      break;
    }
  }

  if (_combo == _query.number_elements()) {
  } else if (_query.number_elements() > _molecules.number_elements()) {
    cerr << "SimpleEnumerator::Initialise:more queries " << _query.size() <<
            " than molecules " << _molecules.size() << " impossible\n";
  } else if (_query.size() == 1) {
    cerr << "Single query will be applied to all sidechains\n";
  } else {
    cerr << "SimpleEnumerator::Initialise:Have " << _query.size() << " queries and " <<
                _molecules.size() << " fragments. Impossible, must be fewer fragments than queries\n";
    return 0;
  }

  if (cl.option_present('p')) {
    _write_starting_molecule = 1;
    if (_verbose) {
      cerr << "Will write the starting molecule\n";
    }
  }

  if (cl.option_present('z')) {
    _ignore_molecules_not_matchng_queries = 1;
    if (_verbose) {
      cerr << "Will ignore molecules not matching any query\n";
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_products_per_starting_molecule) || _max_products_per_starting_molecule < 1) {
      cerr << "The max products per starting molecule (-x) option must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will generate a max of " << _max_products_per_starting_molecule <<
              " products per starting molecule\n";
    }
  }

  if (cl.option_present('h')) {
    _fragments_identical = 1;
  }

  return 1;
}

/* Idea from
https://stackoverflow.com/questions/28711797/generating-n-choose-k-permutations-in-c

void PermGenerator(int n, int k)
{
    std::vector<int> d(n);
    std::iota(d.begin(),d.end(),1);
    cout << "These are the Possible Permutations: " << endl;
    do
    {
        for (int i = 0; i < k; i++)
        {
            cout << d[i] << " ";
        }
        cout << endl;
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
}
*/

int
SimpleEnumerator::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  if (_write_starting_molecule) {
    output << m.smiles() << ' ' << m.name() << '\n';
  }

  Molecule_to_Match target(&m);

  Substructure_Results sresults;

  Set_of_Atoms e0;
  e0.reserve(_query.size());
  for (Substructure_Query* q : _query) {
    if (q->substructure_search(target, sresults) == 0) {
      if (_verbose) {
        cerr << m.smiles() << ' ' << m.name() << " no query match " << q->comment() << '\n';
      }
      return _ignore_molecules_not_matchng_queries;
    }

    // First atom of each embedding
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      e0 << e->front();
    }
  }

  cerr << "Find " << e0.size() << " matched atoms\n";

  if (e0.size() < _molecules.size()) {
    cerr << "SimpleEnumerator::Process:not enough query matches " << e0.size() << 
            " for " << _molecules.size() << " molecules\n";
    return _ignore_not_enough_matched_atoms;
  }

  _made_this_molecule = 0;

  AddSingleFragments(m, e0, output);

#ifdef DEBUG_COMINATORICS
  std::vector<int> foo(10);
  std::iota(foo.begin(), foo.end(), 0);
  do {
    for (int i = 0; i < 3; ++i) {
      cerr << ' ' << foo[i];
    }
    cerr << '\n';
    std::reverse(foo.begin() + 3, foo.end());
  } while (std::next_permutation(foo.begin(), foo.end()));
#endif

  std::vector<int> d(e0.size(), 0);
  std::iota(d.begin(), d.end(), 0);
  const int k = _combo;
  do {
    Set_of_Atoms x;
    for (int i = 0;i < k; ++i) {
      x << e0[d[i]];
    }
    Process(m, x, output);
    if (_made_this_molecule > _max_products_per_starting_molecule) {
      break;
    }
    std::reverse(d.begin() + k, d.end());
  }  while (std::next_permutation(d.begin(), d.end()));
  
  return 1;
}

// Make a single bond btw `a1` and `a2` and unset any implicit
// Hydrogen information on each one.
int
AddBond(Molecule& m, atom_number_t a1, atom_number_t a2) {
  m.add_bond(a1, a2, SINGLE_BOND);
  m.unset_all_implicit_hydrogen_information(a1);
  m.unset_all_implicit_hydrogen_information(a2);

  return 1;
}

int
SimpleEnumerator::AddSingleFragments(Molecule& m,
                                 const Set_of_Atoms& matched_atoms,
                                 IWString_and_File_Descriptor& output) {
  if (matched_atoms.size() != 200005) {
    return 1;
  }

  for (const SetOfMolecules* mols : _molecules) {
    for (const Molecule * frag : *mols) {
      for (atom_number_t a : matched_atoms) {
        const int save_natoms = m.natoms();
        m.add_molecule(frag);
        AddBond(m, a, save_natoms);
        ProcessNewMolecule(m, output);
        m.resize(save_natoms);
      }
    }
  }

  return 1;
}

bool
AscendingOrder(const Set_of_Atoms& matched_atoms) {
  int n = matched_atoms.number_elements();
  if (n == 2 && matched_atoms[0] > matched_atoms[1]) {
    return false;
  }

  for (int i = 1; i < n; ++i) {
    if (matched_atoms[i-1] > matched_atoms[i]) {
      return 0;
    }
  }

  return 1;
}
//#define DEBUG_SIMPLE_ENUMERATOR_PROCESS
// `matched_atoms` contains atom numbers for the substitution points
int
SimpleEnumerator::Process(Molecule& m, const Set_of_Atoms& matched_atoms,
                 IWString_and_File_Descriptor& output) {

#ifdef DEBUG_SIMPLE_ENUMERATOR_PROCESS
  cerr << m.name() <<  " matched atoms";
  for (int i : matched_atoms) {
    cerr << ' ' << i;
  }
  cerr << '\n';
#endif
  if (_fragments_identical && ! AscendingOrder(matched_atoms)) {
    return 0;
  }

  if (_just_one_fragment_everywhere) {
    return ProcessSingleFragments(m, matched_atoms, output);
  }

  return ProcessMultipleFragments(m, matched_atoms, output);
}

int
SimpleEnumerator::ProcessSingleFragments(Molecule& m, const Set_of_Atoms& matched_atoms,
                        IWString_and_File_Descriptor& output) {
  const int initial_matoms = m.natoms();
  const int n = _molecules.number_elements();

  for (int i = 0; i < n; ++i) {
    int save_natoms = m.natoms();
    m.add_molecule(_molecules[i]->at(0));
    AddBond(m, matched_atoms[i], save_natoms);
  }

  int rc = ProcessNewMolecule(m, output);

  m.resize(initial_matoms);

  return rc;
}

int
SimpleEnumerator::ProcessMultipleFragments(Molecule& m, const Set_of_Atoms& matched_atoms,
                        IWString_and_File_Descriptor& output) {
  std::vector<uint32_t> count(_combo);
  for (int i = 0; i < _combo; ++i) {
    count[i] = _molecules[i]->size();
  }

  std::vector<uint32_t> state(_combo, 0);

  combinations::Combinations<uint32_t> c(count);

  do {
    MakeMolecule(m, matched_atoms, state, output);
    if (_made_this_molecule > _max_products_per_starting_molecule) {
      return 1;
    }
  } while (c.Next(state));

  return 1;
}

int
SimpleEnumerator::MakeMolecule(Molecule& m, const Set_of_Atoms& matched_atom,
                const std::vector<uint32_t>& frag,
                IWString_and_File_Descriptor& output) {
  const int initial_matoms = m.natoms();

  for (int i = 0; i < _combo; ++i) {
    const int save_natoms = m.natoms();

    atom_number_t ma = matched_atom[i];

    if (m.hcount(ma) == 0) {
      cerr << "SimpleEnumerator::MakeMolecule:no H on " << m.smarts_equivalent_for_atom(ma) <<
           ' ' << m.name() << " atom " << ma << '\n';
      return 0;
    }

    m.add_molecule(_molecules[i]->at(frag[i]));
    from_python::AddBond(m, ma, save_natoms);
  }

  int rc = ProcessMultipleFragments(m, matched_atom, output);

  m.resize(initial_matoms);

  return rc;
}

int
SimpleEnumerator::ProcessNewMolecule(Molecule& m, IWString_and_File_Descriptor& output) {
  if (! m.valence_ok()) {
    cerr << "SimpleEnumerator::ProcessNewMolecule:invalid valence\n";
    cerr << m.smiles() << ' ' << m.name() << '\n';

    ++_invalid_valence_skipped;

    return 0;
  }

  if (const auto iter = _seen.find(m.unique_smiles()); iter != _seen.end()) {
    ++_duplicates_suppressed;
    return 0;
  }

  _seen.insert(m.unique_smiles());

  ++_made_this_molecule;

  output << m.smiles() << ' ' << m.name() << '.' << _made_this_molecule << '\n';

  ++_molecules_written;

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
SimpleEnumerator::Report(std::ostream& output) const {
  output << "Processed " << _molecules_processed << " molecules, wrote " <<
            _molecules_written << '\n';

  output << _duplicates_suppressed << " duplicates suppressed\n";
  output << _invalid_valence_skipped << " invalid valence suppressed\n";
  return 1;
}

int
PositionalAnalogueScanning(Molecule& m,
                SimpleEnumerator& options,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int PositionalAnalogueScanning(data_source_and_type<Molecule>& input,
                SimpleEnumerator& options,
                IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    if (! PositionalAnalogueScanning(*m, options, output)) {
      cerr << m->smiles() << ' ' << m->name() << " failed\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
PositionalAnalogueScanning(const char* fname, SimpleEnumerator& options,
                        IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(fname);
  if (! input.good()) {
    cerr << "PositionalAnalogueScanning:cannot open '" << fname << "'\n";
    return 0;
  }

  return PositionalAnalogueScanning(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:Eq:c:z:x:pF:h");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage();
  }

  const int verbose = cl.option_present('v');

  SimpleEnumerator options;

  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage();
  }

  if (cl.empty()) {
    cerr << "Must provile name of input file\n";
    Usage();
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    if (! PositionalAnalogueScanning(fname, options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

} // namespace from_python

namespace positional_analogue_scanning {

using std::cerr;
using internal::SetOfMolecules;

class SetOfQueries {
  private:
    resizable_array_p<Substructure_Query> _query;

  public:
    // `fname` is passed to queries_from_file to populate `_query`.
    int Build(IWString& fname);

    int number_queries() const {
      return _query.number_elements();
    }

    int LabelMatchedAtoms(Molecule_to_Match& target, int* matched,
                          int flag, int first_match_only,
                          Set_of_Atoms& destination);
};

int
SetOfQueries::Build(IWString& fname) {
  static constexpr int kInheritDirectoryPath = 1;
  static constexpr int kVerbose = 0;
  if (!  queries_from_file(fname, _query, kInheritDirectoryPath, kVerbose)) {
    cerr << "SetOfQueries::Build:cannot read queries from '" << fname << "'\n";
    return 0;
  }

  return _query.number_elements();
}

// For each query match, fetch the first matched atom, and set the corresponding value
// in `matched` to `flag`.
// Make sure that only one set of queries matches any particular atom.
int
SetOfQueries::LabelMatchedAtoms(Molecule_to_Match& target, int* matched, int flag, int first_match_only,
                Set_of_Atoms& destination) {
  int rc = 0;
  for (Substructure_Query* q : _query) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      atom_number_t a = e->first();
      destination << a;

//    if (matched[a] == 0) {
//      matched[a] = flag;
//      destination << a;
//    } 
    }

    if (first_match_only) {
      return 1;
    }
    ++rc;
  }

  return rc;
}

void
Usage() {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(Positional analogue implementation.
 -c <n>         number of combinations to generate for each set of substituents, usually 1, 2 or 3.
 -q <file>      file with queries specifing matched atoms. 
 -z i           ignore molecules not matching any of the queries.
 -p             write the parent molecule to the output.
 -x <max>       max number of products per starting molecule.
 -v             verbose output.
)";
  // clang-format off

  ::exit(0);
}

class Options {
  private:
    int _verbose;

    int _combo;

    int _nq;
    std::unique_ptr<SetOfQueries[]> _query;

    std::unique_ptr<SetOfMolecules[]> _molecule;

    // For each query, the atoms matched;
    std::unique_ptr<Set_of_Atoms[]> _atoms;

    uint64_t _molecules_processed;
    uint64_t _molecules_written;

    absl::flat_hash_set<IWString> _seen;

    // By default, we do NOT write the starting molecule.
    int _write_starting_molecule;

    // We need to assign a unique name to each product molecule.
    // This is sequential with each starting molecule.
    uint32_t _made_this_molecule;

    extending_resizable_array<uint32_t> _generated;

    int _ignore_no_query_matches;

    uint32_t _max_products_per_starting_molecule;

  // Private functions
    int MaybeOutput(Molecule& m, IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                 const Set_of_Atoms& matched_atoms,
                 int istart,
                 IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                 const std::vector<uint32_t>& indices,
                 IWString_and_File_Descriptor& output);

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _combo = 2;
  _nq = 0;
  _molecules_processed = 0;
  _molecules_written = 0;
  _made_this_molecule = 0;
  _ignore_no_query_matches = 0;
  _max_products_per_starting_molecule = std::numeric_limits<uint32_t>::max();
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('c')) {
    if (! cl.value('c', _combo) || _combo < 1) {
      cerr << "The number of simultaneous attachments (-c) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will generate " << _combo << " way changes\n";
    }
  }

  if (! cl.option_present('q')) {
    cerr << "Options::Initialise:must specify one or more queries for the attachment points via the -q option\n";
    return 0;
  }

  _nq = cl.option_count('q');
  _query.reset(new SetOfQueries[_nq]);
  _atoms.reset(new Set_of_Atoms[_nq]);

  if (cl.option_present('q')) {
    IWString fname;
    for (int i = 0; cl.value('q', fname, i); ++i) {
      if (! _query[i].Build(fname)) {
        cerr << "Options::Initialise:cannot read queries from '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      int totalq = 0;
      for (int i = 0; i < _nq; ++i) {
        totalq += _query[i].number_queries();
      }

      cerr << "Read " << _nq << " query sets containing " << totalq << " substructure queries\n";
    }
  }

  if ((cl.number_elements() - 1) != _nq) {
    cerr << "Options::Initialise:read " << _nq << " sets of queries but " << (cl.size() - 1) << 
                " input files on command line, impossible, must be the same\n";
//  return 0;
  }

  cerr << _nq << " nq\n";
  _molecule.reset(new SetOfMolecules[_nq]);

  for (int i = 1; i < cl.number_elements(); ++i) {
    IWString fname(cl[i]);
    cerr << "i = " << i << " file " << fname << "\n";
    if (! _molecule[i - 1].Build(fname)) {
      cerr << "Cannot read molecules from '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('z')) {
    _ignore_no_query_matches = 1;
    if (_verbose) {
      cerr << "Will ignore molecules not matching any query\n";
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_products_per_starting_molecule) || _max_products_per_starting_molecule < 1) {
      cerr << "The max products per starting molecule (-x) option must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will generate a max of " << _max_products_per_starting_molecule <<
              " products per starting molecule\n";
    }
  }

  if (cl.option_present('p')) {
    _write_starting_molecule = 1;
    if (_verbose) {
      cerr << "Will write the starting molecule\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Options:Read " << _molecules_processed << " wrote " << _molecules_written << '\n';

  for (int i = 0; i < _generated.number_elements() ; ++i) {
    if (_generated[i]) {
      output << _generated[i] << " starting molecules made " << i << " products\n";
    }
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  _made_this_molecule = 0;

  if (_write_starting_molecule) {
    output << m.smiles() << ' ' << m.name() << '\n';
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> label = std::make_unique<int[]>(matoms);
  std::fill_n(label.get(), matoms, 0);

  Molecule_to_Match target(&m);

  static constexpr int kFirstMatchOnly = 1;
  for (int i = 0; i < _nq; ++i) {
    if (! _query[i].LabelMatchedAtoms(target, label.get(), i + 1, kFirstMatchOnly, _atoms[i])) {
      // cerr << "Options::Process:no query matches " << m.name() << " query set " << i << '\n';
      return _ignore_no_query_matches;
    }
    // cerr << "query " << i << " " << _atoms[i] << '\n';
  }

  std::vector<uint32_t> count(_nq);
  for (int i = 0; i < _nq; ++i) {
    count[i] = _atoms[i].size();
  }
  combinations::Combinations<uint32_t> indices(count);

  std::vector<uint32_t> state;
  state.resize(_nq, 0);

  do {
    Process(m, state, output);
    if (_made_this_molecule > _max_products_per_starting_molecule) {
      break;
    }
  } while (indices.Next(state));

  ++_generated[_made_this_molecule];

  return 1;
}
#ifdef N_CHOOSE_K
from https://stackoverflow.com/questions/28711797/generating-n-choose-k-permutations-in-c
void PermGenerator(int n, int k)
{
    std::vector<int> d(n);
    std::iota(d.begin(),d.end(),1);
    cout << "These are the Possible Permutations: " << endl;
    do
    {
        for (int i = 0; i < k; i++)
        {
            cout << d[i] << " ";
        }
        cout << endl;
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
}

#endif

int
Options::Process(Molecule& m,
                 const std::vector<uint32_t>& indices,
                 IWString_and_File_Descriptor& output) {
  Set_of_Atoms matched_atoms;
  for (int i = 0; i < _nq; ++i) {
    int j = indices[i];
    // cerr << "Query " << i << " atom " << j << '\n';
    if (! matched_atoms.add_if_not_already_present(_atoms[i][j])) {
      return 0;
    }
  }

  return Process(m, matched_atoms, 0, output);
}

int
Options::Process(Molecule& m,
                 const Set_of_Atoms& matched_atoms,
                 int istart,
                 IWString_and_File_Descriptor& output) {
  atom_number_t zatom = matched_atoms[istart];

  if (m.hcount(zatom) == 0) {
    cerr << m.smiles() << " no available H on atom " << zatom <<
            m.smarts_equivalent_for_atom(zatom) << '\n';
    return 0;
  }

  const int initial_natoms = m.natoms();
  // cerr << "Starting moleculehas " << initial_natoms << " atoms\n";

  for (int i = istart; i < _nq; ++i) {
    for (Molecule* frag : _molecule[i]) {
      m.add_molecule(frag);
      from_python::AddBond(m, zatom, m.natoms() - 1);
      MaybeOutput(m, output);
      if (istart < (_nq - 1)) {
        Process(m, matched_atoms, istart + 1, output);
      }
      m.resize(initial_natoms);
    }
    if (_made_this_molecule > _max_products_per_starting_molecule) {
      return 1;
    }
  }
  // cerr << "Finished iteration\n";

  return 1;
}

int
Options::MaybeOutput(Molecule& m, IWString_and_File_Descriptor& output) {
  if (! m.valence_ok()) {
    return 0;
  }

  if (const auto iter = _seen.find(m.unique_smiles()); iter != _seen.end()) {
    return 0;
  }

  _seen.insert(m.unique_smiles());

  ++_made_this_molecule;

  output << m.smiles() << ' ' << m.name() << '.' << _made_this_molecule << '\n';

  ++_molecules_written;

  return 1;
}

int
PositionalAnalogueScanning(Molecule& m,
                Options& options,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
PositionalAnalogueScanning(data_source_and_type<Molecule>& input,
                Options& options,
                IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    if (! PositionalAnalogueScanning(*m, options, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
PositionalAnalogueScanning(const char* fname, Options& options,
                        IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(fname);
  if (! input.good()) {
    cerr << "PositionalAnalogueScanning:cannot open '" << fname << "'\n";
    return 0;
  }

  return PositionalAnalogueScanning(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:Eq:c:z:x:ph");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage();
  }

  const int verbose = cl.option_present('v');

  Options options;

  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage();
  }

  if (cl.empty()) {
    cerr << "Insufficent arguments\n";
    Usage();
  }

  IWString_and_File_Descriptor output(1);

  if (! PositionalAnalogueScanning(cl[0], options, output)) {
    return 1;
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace pos

int 
main(int argc, char ** argv) {
  IWString tmp(argv[0]);

  int rc;
  if (tmp == "positional_analogue_scanning") {
    rc = positional_analogue_scanning::Main(argc, argv);
  } else {
    rc = from_python::Main(argc, argv);
  }


  return rc;
}
